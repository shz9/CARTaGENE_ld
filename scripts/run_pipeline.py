#!/usr/bin/env python3
"""Prepare local dependencies, validate requested chromosomes, then submit LD jobs."""

import argparse
import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import urllib.request
from pathlib import Path

import pandas as pd
from annotate_bim import load_genetic_map
from magenpy import __version__ as magenpy_version
from magenpy.stats.ld import c_utils  # Verify the compiled native LD extension imports.
from prepare_blocks import prepare as prepare_blocks
from prepare_panel import chromosome, identity, sha256
from run_chromosome import (
    ROOT,
    load_keep,
    local_references,
    path,
    read_config,
    sample_ids,
    samples,
    validate_config,
)

DEFAULT_DBSNP_URL = (
    "https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
)


def chromosomes(value):
    result = set()
    for part in value.split(","):
        if not re.fullmatch(r"\d+(?:-\d+)?", part):
            raise argparse.ArgumentTypeError("Use autosome lists/ranges such as 1-3,22")
        bounds = list(map(int, part.split("-")))
        start, stop = bounds[0], bounds[-1]
        if not 1 <= start <= stop <= 22:
            raise argparse.ArgumentTypeError(
                "Chromosomes must be in 1-22 with ascending ranges"
            )
        result.update(range(start, stop + 1))
    return sorted(result)


def command(*args, capture_output=False):
    print("RUN", " ".join(map(str, args)), flush=True)
    return subprocess.run(
        list(map(str, args)),
        check=True,
        cwd=ROOT,
        text=True,
        capture_output=capture_output,
    )


def nonempty(filename):
    if not filename.is_file() or filename.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty file: {filename}")
    with filename.open("rb") as handle:
        handle.read(1)


def validate_pvar(filename, chrom, info_field=None):
    """Inspect the header and first variant without scanning the full chromosome."""
    process = None
    if filename.suffix == ".zst":
        process = subprocess.Popen(
            ["zstd", "-dc", str(filename)], stdout=subprocess.PIPE, text=True
        )
        handle = process.stdout
    else:
        handle = (
            gzip.open(filename, "rt") if filename.suffix == ".gz" else filename.open()
        )
    try:
        line = next((line for line in handle if not line.startswith("##")), "")
        header = line.split()
        if header[:5] != ["#CHROM", "POS", "ID", "REF", "ALT"]:
            raise ValueError(f"Invalid PVAR header: {filename}")
        if info_field and "INFO" not in header:
            raise ValueError(
                f"Imputation quality requested but PVAR has no INFO column: {filename}"
            )
        first = next(handle, "").split()
        if len(first) < 5 or chromosome(first[0]) != chrom or int(first[1]) < 1:
            raise ValueError(f"Invalid or empty chromosome {chrom} PVAR: {filename}")
    finally:
        if process and process.poll() is None:
            process.terminate()  # Only the header was requested; do not decompress the entire file.
        handle.close()
        if process:
            process.wait()


def validate_genotypes(cfg, chroms):
    metadata = path(cfg["metadata"])
    nonempty(metadata)
    keep = load_keep(path(cfg["sample_keep"])) if cfg["sample_keep"] else None
    # Validate metadata even if no chromosome files are accessible.
    samples(metadata, set(), cfg["ancestries"], keep)
    report, errors = {}, []
    for chrom in chroms:
        try:
            prefix = path(cfg["input_prefix"], chrom)
            pgen, psam = (Path(str(prefix) + ext) for ext in (".pgen", ".psam"))
            pvars = [
                Path(str(prefix) + ext) for ext in (".pvar", ".pvar.gz", ".pvar.zst")
            ]
            pvars = [p for p in pvars if p.is_file()]
            if len(pvars) != 1:
                raise ValueError(f"Expected exactly one PVAR, found {pvars}")
            pvar = pvars[0]
            for filename in (pgen, pvar, psam):
                nonempty(filename)
            with pgen.open("rb") as handle:
                header = handle.read(3)
            if len(header) < 3 or header[:2] != b"\x6c\x1b":
                raise ValueError(f"Invalid PGEN magic header: {pgen}")
            if pvar.suffix == ".zst" and not shutil.which("zstd"):
                raise ValueError("zstd is required for .pvar.zst files")
            validate_pvar(pvar, chrom, cfg["info_field"])
            available = sample_ids(psam)
            selected, unmatched = samples(metadata, available, cfg["ancestries"], keep)
            counts = {g: len(pairs) for g, pairs in selected.items()}
            if not any(n >= cfg["min_samples"] for n in counts.values()):
                raise ValueError(
                    f"No ancestry has {cfg['min_samples']} eligible samples: {counts}"
                )
            report[str(chrom)] = {
                "files": [identity(p) for p in (pgen, pvar, psam)],
                "samples": counts,
                "unmatched_metadata": unmatched,
            }
            print(
                f"chr{chrom}: genotype headers valid; matched samples {counts}",
                flush=True,
            )
        except (OSError, ValueError, EOFError, pd.errors.ParserError) as exc:
            errors.append(f"chr{chrom}: {exc}")
    if errors:
        raise ValueError(
            "Genotype validation failed; no jobs submitted:\n" + "\n".join(errors)
        )
    return report


def panel_current(panel, rsids, vcf):
    """Lightweight receipt/stat check; full scans belong on compute nodes."""
    if not panel.is_file() or not panel.stat().st_size:
        return False
    if not rsids:
        return True  # An independently prepared panel may omit its source lists.
    receipt = panel.with_suffix('.json')
    if not receipt.is_file():
        return False
    try:
        report = json.loads(receipt.read_text())
    except (ValueError, OSError):
        return False
    recorded = [{k: row.get(k) for k in ('path', 'size', 'mtime_ns')}
                for row in report.get('lists', [])]
    return (report.get('assembly') == 'GRCh38'
            and recorded == [identity(p) for p in rsids]
            and report.get('panel') == identity(panel)
            and (not vcf.exists() or report.get('dbsnp') == identity(vcf)))


def panel_job_active(job_id):
    queued = subprocess.run(['squeue', '--noheader', '--jobs', job_id, '--format=%T'],
                            text=True, capture_output=True, check=True)
    return bool(queued.stdout.strip())


def submit_panel(panel, rsids, vcf):
    inputs = {'panel': str(panel.resolve()), 'vcf': identity(vcf),
              'lists': [identity(p) for p in rsids],
              'code': sha256(ROOT / 'scripts/prepare_panel.py')}
    token = hashlib.sha256(str(panel.resolve()).encode()).hexdigest()[:16]
    logs = ROOT / 'logs'
    logs.mkdir(parents=True, exist_ok=True)
    receipt = logs / f'panel-{token}.json'
    if receipt.is_file():
        state = json.loads(receipt.read_text())
        if panel_job_active(state['job_id']):
            print(f'Panel preparation job {state["job_id"]} is still queued/running. '
                  'Please wait until it completes, then rerun the same pipeline command.', flush=True)
            return
    submitted = command('sbatch', '--parsable', ROOT / 'scripts/prepare_panel.sbatch',
                        vcf, panel, *rsids, capture_output=True)
    job_id = submitted.stdout.strip().split(';')[0]
    if not job_id.isdigit():
        raise ValueError(f'Unexpected sbatch response: {submitted.stdout!r}')
    state = dict(inputs, job_id=job_id, status='awaiting_panel')
    receipt.write_text(json.dumps(state, indent=2) + '\n')
    (logs / 'preflight.json').write_text(json.dumps(state, indent=2) + '\n')
    print(f'Submitted panel preparation job {job_id}. Please wait until it completes, '
          f'check logs/panel-{job_id}.out, then rerun the same pipeline command. '
          'No LD jobs have been submitted.', flush=True)


def download_dbsnp(url, destination):
    """Download on login; save the expected MD5 for verification in the panel job."""
    if not shutil.which("curl"):
        raise FileNotFoundError(
            "curl is required to download dbSNP; or supply a local --dbsnp-vcf"
        )
    with urllib.request.urlopen(url + ".md5", timeout=60) as response:
        checksum = re.search(r"\b[a-fA-F0-9]{32}\b", response.read().decode())
    if checksum is None:
        raise ValueError("Invalid dbSNP MD5 checksum response")
    checksum = checksum.group()
    destination.parent.mkdir(parents=True, exist_ok=True)
    partial = destination.with_name(destination.name + ".part")
    command(
        "curl",
        "--fail",
        "--location",
        "--retry",
        "3",
        "--continue-at",
        "-",
        "--output",
        partial,
        url,
    )
    partial.replace(destination)
    Path(str(destination) + ".md5").write_text(f"{checksum}  {destination.name}\n")


def prepare_references(cfg, chroms, args):
    panel, vcf = path(cfg["panel"]), path(args.dbsnp_vcf)
    rsids = [path(p) for p in (args.rsids if args.rsids is not None else cfg["rsids"])]
    for filename in rsids:
        nonempty(filename)
    if not panel_current(panel, rsids, vcf):
        if args.check_only:
            raise ValueError(
                "Prepared panel is missing or stale; rerun without --check-only"
            )
        if not rsids:
            raise ValueError("Supply --rsids FILE [FILE...] to build the missing panel")
        if not vcf.is_file():
            download_dbsnp(args.dbsnp_url, vcf)
        submit_panel(panel, rsids, vcf)
        return False  # Stop after scheduling preparation; never submit LD jobs yet.
    maps = {c: path(cfg["genetic_map"], c) for c in chroms}
    if not all(p.is_file() for p in maps.values()):
        if args.check_only:
            raise FileNotFoundError("Genetic maps missing; rerun without --check-only")
        cache = ROOT / "data/genetic_map"
        command(
            sys.executable,
            ROOT / "scripts/prepare_genetic_map.py",
            "--output-dir",
            cache,
        )
        for chrom, destination in maps.items():
            source = cache / f"chr_{chrom}.tsv"
            destination.parent.mkdir(parents=True, exist_ok=True)
            if source.resolve() != destination.resolve():
                shutil.copyfile(source, destination)
    if cfg["estimator"] == "block" and not args.check_only:
        # Download each source once, retaining the existing provenance/cache behavior.
        cache = ROOT / "data/blocks"
        prepare_blocks(cache, cfg["ancestries"], cfg["block_ancestry_map"])
        for group in cfg["ancestries"]:
            for chrom in chroms:
                destination = path(cfg["blocks"].replace("{ancestry}", group), chrom)
                source = cache / f"{group}.bed"
                destination.parent.mkdir(parents=True, exist_ok=True)
                if source.resolve() != destination.resolve():
                    shutil.copyfile(source, destination)
    for chrom in chroms:
        local_references(cfg, chrom)
        load_genetic_map(maps[chrom])
        with panel.open() as handle:
            if handle.readline().split() != ['CHROM', 'POS', 'REF', 'ALT', 'RSID']:
                raise ValueError(f'Invalid prepared panel header: {panel}')
        print(f'chr{chrom}: local reference headers valid; full panel validation runs in the LD job', flush=True)
    if cfg["estimator"] == "block":
        # Prevent a config edit from silently using blocks prepared under an old mapping.
        manifest = ROOT / "data/blocks/summary.json"
        if not manifest.is_file():
            raise ValueError(
                "Block provenance missing; run preparation before submitting"
            )
        records = json.loads(manifest.read_text())["ancestries"]
        for group in cfg["ancestries"]:
            record = records[group]
            for chrom in chroms:
                blockfile = path(cfg["blocks"].replace("{ancestry}", group), chrom)
                if record["source_ancestry"] != cfg["block_ancestry_map"][
                    group
                ] or record["sha256"] != sha256(blockfile):
                    raise ValueError(
                        f"Stale/mismatched blocks for {group}; run preparation again"
                    )

    return True


def writable_directory(directory):
    directory.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryFile(dir=directory):
        pass


def main():
    if os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Run this script on the login node, outside a SLURM job.")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/pipeline.json")
    parser.add_argument("--chromosomes", type=chromosomes, default=chromosomes("1-22"))
    parser.add_argument("--max-concurrent", type=int, default=4)
    parser.add_argument("--rsids", nargs="+", help="Override config rsids (default: metadata/rsids.txt)")
    parser.add_argument("--dbsnp-vcf", default="data/reference/GCF_000001405.40.gz")
    parser.add_argument("--dbsnp-url", default=DEFAULT_DBSNP_URL)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check-only", action="store_true")
    mode.add_argument("--prepare-only", action="store_true", help="Prepare references (may submit a panel job); do not submit LD jobs")
    args = parser.parse_args()
    if args.max_concurrent < 1:
        parser.error("--max-concurrent must be positive")
    os.chdir(ROOT)
    config_file = path(args.config).resolve()
    cfg = read_config(config_file)
    validate_config(cfg)
    if magenpy_version != "0.2.5":
        raise ValueError("Expected magenpy 0.2.5")
    plink = shutil.which("plink2") or shutil.which("plink")
    if not plink or "PLINK v2" not in subprocess.check_output(
        [plink, "--version"], text=True
    ):
        raise ValueError("Load a PLINK2 module before running the pipeline")
    if not args.check_only and not shutil.which("sbatch"):
        raise FileNotFoundError("sbatch is unavailable")
    genotypes = validate_genotypes(cfg, args.chromosomes)
    if not prepare_references(cfg, args.chromosomes, args):
        return
    writable_directory(path(cfg["output_dir"]))
    writable_directory(ROOT / "logs")
    if cfg["scratch_dir"]:
        writable_directory(path(cfg["scratch_dir"]))
    report = {
        "config": cfg,
        "chromosomes": args.chromosomes,
        "genotypes": genotypes,
        "python": sys.version,
        "python_module": os.environ.get("PYTHON_MODULE", "python/3.13"),
        "magenpy": magenpy_version,
        "status": "validated",
    }
    receipt = ROOT / "logs/preflight.json"
    receipt.write_text(json.dumps(report, indent=2) + "\n")
    print(f"All requested chromosomes validated. Report: {receipt}", flush=True)
    if not args.prepare_only and not args.check_only:
        # Freeze the validated settings rather than rereading a mutable config on workers.
        settings = json.dumps(cfg, indent=2) + "\n"
        token = hashlib.sha256(settings.encode()).hexdigest()[:16]
        snapshot = ROOT / f"logs/config-{token}.json"
        snapshot.write_text(settings)
        array = ",".join(map(str, args.chromosomes)) + f"%{args.max_concurrent}"
        submitted = command(
            "sbatch",
            "--parsable",
            f"--array={array}",
            ROOT / "scripts/chromosome.sbatch",
            snapshot,
            capture_output=True,
        )
        report.update(
            status="submitted",
            job_id=submitted.stdout.strip(),
            config_snapshot=str(snapshot),
        )
        receipt.write_text(json.dumps(report, indent=2) + "\n")
        print(f"Submitted SLURM job {report['job_id']}", flush=True)
    else:
        print("No jobs submitted.", flush=True)


if __name__ == "__main__":
    main()
