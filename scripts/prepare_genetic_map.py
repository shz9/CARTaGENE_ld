#!/usr/bin/env python3
"""Download and normalize the Beagle GRCh38 genetic map."""

import argparse
import json
import os
import pathlib
import urllib.request
import zipfile

import pandas as pd

from annotate_bim import validate_genetic_map


ROOT = pathlib.Path(__file__).resolve().parent.parent
DEFAULT_OUTPUT_DIR = ROOT / "data" / "genetic_map"
DEFAULT_URL = (
    "https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/"
    "plink.GRCh38.map.zip"
)
AUTOSOMES = tuple(range(1, 23))


def fetch(url, path):
    if path.is_file() and path.stat().st_size > 0:
        print(f"Reusing {path}")
        return

    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".part")
    print(f"Downloading {url}", flush=True)
    urllib.request.urlretrieve(url, temporary)
    temporary.replace(path)


def find_members(names):
    """Resolve every autosomal map member in one archive-name pass."""
    members = {}
    for name in names:
        if "no_chr_in_chrom_field/" not in name:
            continue
        for chrom in AUTOSOMES:
            if name.endswith(f"/plink.chr{chrom}.GRCh38.map"):
                if chrom in members:
                    raise RuntimeError(
                        f"Duplicate GRCh38 map for chromosome {chrom}"
                    )
                members[chrom] = name
                break

    missing = [chrom for chrom in AUTOSOMES if chrom not in members]
    if missing:
        raise RuntimeError(f"Missing GRCh38 maps for chromosomes: {missing}")
    return members


def statistics(table):
    return {'markers': len(table), 'first': table.head(1).to_dict('records')[0], 'last': table.tail(1).to_dict('records')[0]}


def normalize_map(source, destination, chrom):
    with source:
        table = pd.read_csv(source, sep=r'\s+', header=None)
    table.columns = ['chrom', 'id', 'cm', 'position']
    if not table.chrom.eq(chrom).all():
        raise RuntimeError(f'Unexpected chromosome in chromosome {chrom} map')
    table = validate_genetic_map(table[['position', 'cm']])
    temporary = destination.with_name(destination.name + '.tmp')
    table.to_csv(temporary, sep='\t', index=False)
    temporary.replace(destination)
    return statistics(table)


def map_statistics(filename):
    return statistics(validate_genetic_map(pd.read_csv(filename, sep='\t')))


def main():
    if os.environ.get("SLURM_JOB_ID"):
        raise SystemExit("Prepare genetic maps directly on the login node, outside a SLURM job.")
    parser = argparse.ArgumentParser(
        description="Prepare the Beagle/HapMap GRCh38 genetic map"
    )
    parser.add_argument("--archive", type=pathlib.Path)
    parser.add_argument("--url", default=DEFAULT_URL)
    parser.add_argument(
        "--output-dir",
        type=pathlib.Path,
        default=DEFAULT_OUTPUT_DIR,
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.archive:
        archive = args.archive.expanduser().resolve()
        if not archive.is_file():
            raise RuntimeError(f"Missing {archive}")
    else:
        archive = output_dir / "plink.GRCh38.map.zip"
        fetch(args.url, archive)

    summary = {
        "assembly": "GRCh38",
        "source": args.url,
        "archive": str(archive),
        "interpolation": "linear; nearest endpoint outside map range",
        "chromosomes": {},
    }

    with zipfile.ZipFile(archive) as zipped:
        members = find_members(zipped.namelist())

        for index, chrom in enumerate(AUTOSOMES, start=1):
            destination = output_dir / f"chr_{chrom}.tsv"
            if args.force or not destination.is_file():
                statistics = normalize_map(
                    zipped.open(members[chrom]), destination, chrom
                )
                action = "normalized"
            else:
                statistics = map_statistics(destination)
                action = "validated existing"
            summary["chromosomes"][str(chrom)] = statistics
            print(
                f"[{index}/{len(AUTOSOMES)}] chr{chrom}: "
                f"{statistics['markers']:,} markers ({action})",
                flush=True,
            )

    temporary = output_dir / "summary.json.tmp"
    temporary.write_text(json.dumps(summary, indent=2) + "\n")
    temporary.replace(output_dir / "summary.json")
    print(f"GRCh38 genetic maps: {output_dir}")


if __name__ == "__main__":
    main()
