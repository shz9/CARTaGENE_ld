#!/usr/bin/env python3
"""One chromosome per SLURM task; ancestries processed sequentially."""
import argparse
from collections import Counter
import contextlib
import csv
import gzip
import importlib.metadata
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import zipfile

from annotate_bim import load_genetic_map, interpolate_sorted, format_cm
from prepare_panel import chromosome, identity, sha256, RSID

ROOT = Path(__file__).resolve().parent.parent
GROUPS = {'AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID'}


def path(value, chrom=None):
    value = os.path.expandvars(os.path.expanduser(value))
    if chrom is not None:
        value = value.format(chrom=chrom)
    result = Path(value)
    return result if result.is_absolute() else ROOT / result


@contextlib.contextmanager
def text_file(filename):
    if str(filename).endswith('.zst'):
        process = subprocess.Popen(['zstd', '-dc', str(filename)], stdout=subprocess.PIPE, text=True)
        try:
            yield process.stdout
        finally:
            process.stdout.close()
            code = process.wait()
            if code:
                raise RuntimeError(f'zstd failed for {filename}: {code}')
    elif str(filename).endswith('.gz'):
        with gzip.open(filename, 'rt') as handle:
            yield handle
    else:
        with filename.open() as handle:
            yield handle


def sample_ids(psam):
    with psam.open() as handle:
        header = next(handle).split()
        header[0] = header[0].lstrip('#')
        if 'IID' not in header:
            raise ValueError('PSAM missing IID')
        seen = set()
        for line in handle:
            f = dict(zip(header, line.split()))
            pair = (f.get('FID', '0'), f['IID'])
            if pair in seen:
                raise ValueError(f'Duplicate PSAM sample {pair}')
            seen.add(pair)
    return seen


def samples(metadata, available, groups, keep=None):
    selected = {g: [] for g in groups}
    seen, iids = set(), set()
    unmatched = 0
    with metadata.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'FID', 'IID', 'Ancestry'} <= set(reader.fieldnames or []):
            raise ValueError('Metadata requires FID, IID, Ancestry headers')
        for row in reader:
            pair = (row['FID'], row['IID'])
            if any(not v or any(c.isspace() for c in v) for v in pair):
                raise ValueError(f'Invalid sample ID {pair}')
            if pair in seen or pair[1] in iids:
                raise ValueError(f'Duplicate metadata ID {pair}')
            seen.add(pair)
            iids.add(pair[1])
            group = row['Ancestry']
            if group not in GROUPS:
                raise ValueError(f'Unknown ancestry {group!r}')
            if pair not in available:
                unmatched += 1
                continue
            if group in selected and (keep is None or pair in keep):
                selected[group].append(pair)
    return selected, unmatched


def load_keep(filename):
    with filename.open() as handle:
        if next(handle).split() not in (['FID', 'IID'], ['#FID', 'IID']):
            raise ValueError('Sample keep file requires FID IID header')
        result = set()
        for line in handle:
            fields = line.split()
            if len(fields) != 2:
                raise ValueError('Expected FID IID keep rows')
            result.add(tuple(fields))
    return result


def load_panel(filename, chrom):
    panel = {}
    rsids = set()
    with filename.open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if int(row['CHROM']) != chrom:
                continue
            key = (int(row['POS']), row['REF'], row['ALT'])
            if key[0] < 1 or len(key[1]) != 1 or len(key[2]) != 1 or key[1] not in 'ACGT' or key[2] not in 'ACGT' or key[1] == key[2] or not RSID.fullmatch(row['RSID']):
                raise ValueError(f'Invalid panel row {row}')
            if key in panel or row['RSID'] in rsids:
                raise ValueError('Ambiguous panel; use prepare_panel.py')
            panel[key] = row['RSID']
            rsids.add(row['RSID'])
    if not panel:
        raise ValueError(f'No panel variants on chromosome {chrom}')
    return panel


def match_pvar(filename, panel, chrom, info_field=None, min_info=0.8):
    found, counts, id_counts, key_counts = {}, Counter(), Counter(), Counter()
    with text_file(filename) as handle:
        header = None
        previous = 0
        for line in handle:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.split()
                continue
            if header is None:
                raise ValueError('PVAR requires #CHROM header')
            f = line.split()
            c, pos, old, ref, alt = f[:5]
            if chromosome(c) != chrom:
                raise ValueError('Unexpected chromosome in PVAR')
            pos = int(pos)
            if pos < previous:
                raise ValueError('Unsorted PVAR')
            previous = pos
            key = (pos, ref.upper(), alt.upper())
            rsid = panel.get(key)
            if rsid is None:
                continue
            counts['panel_records'] += 1
            id_counts[old] += 1
            key_counts[key] += 1
            if old == '.':
                continue
            if info_field:
                if 'INFO' not in header:
                    raise ValueError('Imputation quality requested but PVAR has no INFO column')
                info = dict(item.split('=', 1) for item in f[header.index('INFO')].split(';') if '=' in item)
                try:
                    value = float(info.get(info_field, 'nan'))
                except ValueError:
                    value = float('nan')
                if not math.isfinite(value) or value < min_info:
                    counts['info_failed'] += 1
                    continue
            found[old] = (key, rsid)
    result = {old: row for old, row in found.items() if id_counts[old] == 1 and key_counts[row[0]] == 1}
    counts.update(selected=len(result), duplicate_or_missing_records=counts['panel_records'] - len(result) - counts['info_failed'])
    if not result:
        raise ValueError('No panel matches after quality and uniqueness filters')
    return result, dict(counts)


def write_keep(filename, pairs):
    filename.write_text('#FID\tIID\n' + ''.join(f'{fid}\t{iid}\n' for fid, iid in pairs))


def line_count(filename):
    with filename.open() as handle:
        return sum(1 for _ in handle)


def annotate(bim, genetic_map, mapping):
    positions, cms = load_genetic_map(genetic_map)
    counts, seen = Counter(), set()
    temporary = bim.with_suffix('.bim.tmp')
    index, previous = 0, 0
    with bim.open() as source, temporary.open('w') as output:
        for line in source:
            f = line.split()
            key, rsid = mapping[f[1]]
            if int(f[3]) != key[0] or {f[4], f[5]} != {key[1], key[2]}:
                raise ValueError('BIM allele/position mismatch')
            if key[0] < previous or rsid in seen:
                raise ValueError('Unsorted/duplicate BIM')
            previous = key[0]
            seen.add(rsid)
            cm, kind, index = interpolate_sorted(key[0], positions, cms, index)
            counts[kind] += 1
            f[1], f[2] = rsid, format_cm(cm)
            output.write('\t'.join(f) + '\n')
    temporary.replace(bim)
    return dict(counts)


def validate_zip(filename):
    with zipfile.ZipFile(filename) as zipped:
        if '.zmetadata' not in zipped.namelist() or zipped.testzip() is not None:
            raise ValueError(f'Invalid LD ZIP {filename}')


def compute_ld(bfile, output, temporary, cfg, chrom, group, threads):
    import magenpy as mgp
    import zarr
    kwargs = {'cm_window_size': cfg['window_cm']} if cfg['estimator'] == 'windowed' else {'ld_blocks_file': str(path(cfg['blocks'].replace('{ancestry}', group), chrom))}
    loader = mgp.GWADataLoader(bed_files=str(bfile), backend='magenpy', threads=threads,
                               temp_dir=str(temporary), genome_build='GRCh38')
    try:
        loader.compute_ld(cfg['estimator'], str(output), dtype=cfg['storage_dtype'],
                          compressor_name='zstd', compression_level=7, **kwargs)
        matrices = loader.get_ld_matrices()
        if len(matrices) != 1:
            raise ValueError('Expected exactly one LD chromosome')
        for matrix in matrices.values():
            matrix.set_store_attr('Biobank', 'CARTaGENE')
            matrix.set_store_attr('Ancestry', group)
            zarr.consolidate_metadata(matrix._zg.store)
            store = Path(matrix._zarr_store_path(matrix._zg.store))
            archive = output / f'chr_{chrom}.zip'
            with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_STORED, allowZip64=True) as zipped:
                for f in sorted(store.rglob('*')):
                    if f.is_file() and f.name != '.zmetadata':
                        zipped.write(f, f.relative_to(store).as_posix())
                zipped.write(store / '.zmetadata', '.zmetadata')
            validate_zip(archive)
    finally:
        loader.cleanup()
    # Only retain the verified ZIP, not a second copy of the Zarr store.
    for entry in output.iterdir():
        if entry.is_dir():
            shutil.rmtree(entry)


def run(args):
    cfg = json.loads(path(args.config).read_text())
    defaults = json.loads((ROOT / 'config/pipeline.json').read_text())
    unknown = set(cfg) - set(defaults)
    if unknown:
        raise ValueError(f'Unknown configuration keys: {unknown}')
    defaults.update(cfg)
    cfg = defaults
    if not cfg['ancestries'] or len(set(cfg['ancestries'])) != len(cfg['ancestries']) or not set(cfg['ancestries']) <= GROUPS:
        raise ValueError('Invalid ancestries')
    if cfg['estimator'] not in ('windowed', 'block') or (cfg['estimator'] == 'block' and not cfg['blocks']):
        raise ValueError('Use windowed LD or provide block file template with {ancestry}')
    for name in ('maf', 'geno', 'mind', 'hard_call_threshold'):
        if not 0 <= cfg[name] < 0.5:
            raise ValueError(f'Invalid {name}')
    if not 0 < cfg['hwe'] <= 1 or cfg['mac'] < 1 or cfg['min_samples'] < 2 or cfg['window_cm'] <= 0 or args.threads < 1:
        raise ValueError('Invalid QC/LD parameters')
    if cfg['storage_dtype'] not in ('int8', 'int16', 'float32', 'float64'):
        raise ValueError('Invalid storage dtype')
    prefix = path(cfg['input_prefix'], args.chrom)
    pgen, psam = Path(str(prefix) + '.pgen'), Path(str(prefix) + '.psam')
    pvars = [Path(str(prefix) + ext) for ext in ('.pvar', '.pvar.gz', '.pvar.zst')]
    pvars = [p for p in pvars if p.is_file()]
    if len(pvars) != 1:
        raise ValueError(f'Expected exactly one PVAR: {pvars}')
    pvar = pvars[0]
    references = [path(cfg['metadata']), path(cfg['panel']), path(cfg['genetic_map'], args.chrom)]
    if cfg['sample_keep']:
        references.append(path(cfg['sample_keep']))
    if cfg['estimator'] == 'block':
        references += [path(cfg['blocks'].replace('{ancestry}', g), args.chrom) for g in cfg['ancestries']]
    version = subprocess.check_output([args.plink, '--version'], text=True).strip()
    if 'PLINK v2' not in version:
        raise ValueError(f'Expected PLINK2: {version}')
    mgp_version = importlib.metadata.version('magenpy')
    if mgp_version != '0.2.5':
        raise ValueError(f'Expected magenpy 0.2.5, found {mgp_version}')
    fingerprint = {'config': cfg, 'chrom': args.chrom, 'plink': version, 'magenpy': mgp_version,
                   'packages': sorted((d.metadata['Name'], d.version) for d in importlib.metadata.distributions() if d.metadata['Name']),
                   'inputs': [identity(f) for f in (pgen, pvar, psam)],
                   'references': {str(f): sha256(f) for f in references},
                   'code': {f.name: sha256(f) for f in sorted((ROOT / 'scripts').glob('*')) if f.is_file()}}
    import hashlib
    token = hashlib.sha256(json.dumps(fingerprint, sort_keys=True).encode()).hexdigest()[:20]
    parent = path(cfg['output_dir']) / token
    destination = parent / f'chr{args.chrom}'
    parent.mkdir(parents=True, exist_ok=True)
    receipt = destination / 'summary.json'
    if receipt.is_file():
        summary = json.loads(receipt.read_text())
        for group, stats in summary['ancestries'].items():
            if stats.get('status') == 'complete':
                validate_zip(destination / group / 'ld' / f'chr_{args.chrom}.zip')
        print(f'Already complete: {destination}')
        return
    lock = parent / f'chr{args.chrom}.lock'
    lock.mkdir()  # Exclusive ownership. A SIGKILL can leave this lock behind.
    try:
        scratch = cfg['scratch_dir'] or os.environ.get('SLURM_TMPDIR')
        if not scratch:
            raise ValueError('Set scratch_dir or run on SLURM with SLURM_TMPDIR')
        with tempfile.TemporaryDirectory(prefix=f'cag-chr{args.chrom}-', dir=path(scratch)) as tmp:
            work = Path(tmp)
            summary = {'fingerprint': fingerprint, 'ancestries': {}}
            def command(arguments):
                arguments = list(map(str, arguments))
                print('RUN', ' '.join(arguments), flush=True)
                with (work / 'commands.jsonl').open('a') as handle:
                    handle.write(json.dumps(arguments) + '\n')
                subprocess.run(arguments, check=True)
            memory = max(512, int(os.environ.get('SLURM_MEM_PER_NODE', '131072')) * 3 // 4)
            base = [args.plink, '--threads', args.threads, '--memory', memory]
            available = sample_ids(psam)
            keep = load_keep(path(cfg['sample_keep'])) if cfg['sample_keep'] else None
            selected, unmatched = samples(path(cfg['metadata']), available, cfg['ancestries'], keep)
            summary['unmatched_metadata'] = unmatched
            summary['psam_samples'] = len(available)
            summary['requested_panel_variants'] = len(load_panel(path(cfg['panel']), args.chrom))
            all_pairs = [pair for pairs in selected.values() for pair in pairs]
            if not all_pairs:
                raise ValueError('No selected samples intersect PSAM; check FID/IID conventions')
            write_keep(work / 'all.keep', all_pairs)
            panel = load_panel(path(cfg['panel']), args.chrom)
            mapping, summary['panel_match'] = match_pvar(pvar, panel, args.chrom, cfg['info_field'], cfg['min_info'])
            (work / 'extract.txt').write_text(''.join(old + '\n' for old in mapping))
            shared = work / 'shared'
            command([*base, '--pgen', pgen, '--pvar', pvar, '--psam', psam,
                     '--chr', args.chrom, '--keep', work / 'all.keep', '--extract', work / 'extract.txt',
                     '--snps-only', 'just-acgt', '--min-alleles', 2, '--max-alleles', 2,
                     '--hard-call-threshold', cfg['hard_call_threshold'], '--make-bed', '--out', shared])
            summary['shared_variants'] = line_count(work / 'shared.bim')
            summary['shared_samples'] = line_count(work / 'shared.fam')
            for group, pairs in selected.items():
                stats = summary['ancestries'][group] = {'metadata_matched': len(pairs)}
                if len(pairs) < cfg['min_samples']:
                    stats.update(status='skipped', reason='Too few matched samples')
                    continue
                folder = work / group
                folder.mkdir()
                write_keep(folder / 'samples.keep', pairs)
                # Sample filtering on BED, then variant QC on the retained samples.
                command([*base, '--bfile', shared, '--keep', folder / 'samples.keep',
                         '--mind', cfg['mind'], '--make-bed', '--out', folder / 'samples'])
                stats['samples_after_qc'] = line_count(folder / 'samples.fam')
                if stats['samples_after_qc'] < cfg['min_samples']:
                    stats.update(status='skipped', reason='Too few samples after missingness QC')
                    continue
                command([*base, '--bfile', folder / 'samples', '--nonfounders',
                         '--maf', cfg['maf'], '--mac', cfg['mac'], '--geno', cfg['geno'],
                         '--hwe', cfg['hwe'], 'midp', 'keep-fewhet', '--make-bed', '--out', folder / 'qc'])
                stats['variants'] = line_count(folder / 'qc.bim')
                if stats['variants'] < 2:
                    stats.update(status='skipped', reason='Fewer than two QC variants')
                    continue
                # Record source alleles before replacing IDs in BIM.
                with (folder / 'qc.bim').open() as handle, (folder / 'variants.tsv').open('w') as output:
                    output.write('ORIGINAL_ID\tRSID\tCHROM\tPOS\tREF\tALT\n')
                    for line in handle:
                        old = line.split()[1]
                        (pos, ref, alt), rsid = mapping[old]
                        output.write(f'{old}\t{rsid}\t{args.chrom}\t{pos}\t{ref}\t{alt}\n')
                stats['map'] = annotate(folder / 'qc.bim', path(cfg['genetic_map'], args.chrom), mapping)
                compute_ld(folder / 'qc', folder / 'ld', folder / 'temp', cfg, args.chrom, group, args.threads)
                stats['status'] = 'complete'
                if stats['samples_after_qc'] < 1000:
                    stats['warning'] = 'Small ancestry sample; assess LD stability'
            # Keep logs, QC metadata and optionally QC BED. Shared BED stays in scratch.
            for f in list(work.rglob('*')):
                if not f.is_file():
                    continue
                if f.suffix == '.bed' and not (cfg['retain_bed'] and f.name == 'qc.bed'):
                    f.unlink()
                elif f.parent == work and f.suffix not in ('.log', '.jsonl'):
                    f.unlink()
            (work / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
            staging = parent / f'.chr{args.chrom}.partial'
            if staging.exists():
                shutil.rmtree(staging)
            shutil.copytree(work, staging)
            staging.replace(destination)
            print(f'Complete: {destination}', flush=True)
    finally:
        lock.rmdir()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', default='config/pipeline.json')
    parser.add_argument('--chrom', type=int, choices=range(1, 23), required=True)
    parser.add_argument('--plink', default='plink2')
    parser.add_argument('--threads', type=int, default=1)
    run(parser.parse_args())


if __name__ == '__main__':
    main()
