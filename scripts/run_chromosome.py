#!/usr/bin/env python3
"""One chromosome per SLURM task; ancestries processed sequentially."""
import argparse
from collections import Counter
import contextlib
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

import pandas as pd

from annotate_bim import annotate, read_bim, variant_table
from prepare_panel import chromosome, identity, sha256, RSID
from ld_blocks import parse_blocks

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


def read_config(filename):
    defaults = json.loads((ROOT / 'config/pipeline.json').read_text())
    config = json.loads(filename.read_text())
    if set(config) - set(defaults):
        raise ValueError(f'Unknown configuration keys: {set(config) - set(defaults)}')
    defaults.update(config)
    return defaults


def validate_config(cfg):
    if not isinstance(cfg['rsids'], list) or any(not isinstance(p, str) or not p for p in cfg['rsids']):
        raise ValueError('rsids must be a list of input file paths')
    groups = cfg['ancestries']
    if not groups or len(set(groups)) != len(groups) or not set(groups) <= GROUPS:
        raise ValueError('Invalid ancestries')
    if cfg['estimator'] not in {'block', 'windowed'}:
        raise ValueError('Invalid LD estimator')
    if cfg['storage_dtype'] not in {'int8', 'int16', 'float32', 'float64'}:
        raise ValueError('Invalid LD storage dtype')
    for key in ('maf', 'geno', 'mind', 'hard_call_threshold'):
        if not 0 <= cfg[key] < 0.5:
            raise ValueError(f'Invalid {key}')
    if (not 0 < cfg['hwe'] <= 1 or not math.isfinite(cfg['mac']) or cfg['mac'] < 1 or cfg['min_samples'] < 2
            or not math.isfinite(cfg['window_cm']) or cfg['window_cm'] <= 0):
        raise ValueError('Invalid QC/LD parameters')
    if not 0 <= cfg['min_info'] <= 1:
        raise ValueError('Invalid minimum imputation quality')
    if cfg['estimator'] == 'block':
        if not cfg['blocks'] or '{ancestry}' not in cfg['blocks']:
            raise ValueError('Block file template must contain {ancestry}')
        sources = cfg['block_ancestry_map']
        if any(sources.get(g) not in {'AFR', 'EAS', 'EUR', 'SAS'} for g in groups):
            raise ValueError('Missing or invalid block_ancestry_map entry')


def sample_ids(psam):
    table = pd.read_csv(psam, sep=r'\s+', dtype=str, keep_default_na=False)
    table.columns = table.columns.str.lstrip('#')
    if 'IID' not in table:
        raise ValueError('PSAM missing IID')
    if 'FID' not in table:
        table['FID'] = '0'
    if table.duplicated(['FID', 'IID']).any():
        raise ValueError('Duplicate PSAM sample')
    return set(table[['FID', 'IID']].itertuples(index=False, name=None))


def samples(metadata, available, groups, keep=None):
    table = pd.read_csv(metadata, sep='\t', dtype=str, keep_default_na=False)
    if not {'FID', 'IID', 'Ancestry'} <= set(table.columns):
        raise ValueError('Metadata requires FID, IID, Ancestry headers')
    if not table[['FID', 'IID']].apply(lambda col: col.str.fullmatch(r'\S+')).all().all():
        raise ValueError('Invalid sample ID')
    if table.IID.duplicated().any():
        raise ValueError('Duplicate metadata ID')
    if not table.Ancestry.isin(GROUPS).all():
        raise ValueError('Unknown ancestry in metadata')
    pairs = pd.MultiIndex.from_frame(table[['FID', 'IID']])
    matched = pairs.isin(available)
    allowed = matched if keep is None else matched & pairs.isin(keep)
    return {group: list(table.loc[allowed & table.Ancestry.eq(group), ['FID', 'IID']]
                        .itertuples(index=False, name=None)) for group in groups}, int((~matched).sum())


def load_keep(filename):
    table = pd.read_csv(filename, sep=r'\s+', dtype=str, keep_default_na=False)
    table.columns = table.columns.str.lstrip('#')
    if list(table.columns) != ['FID', 'IID']:
        raise ValueError('Sample keep file requires FID IID header')
    if not table.apply(lambda col: col.str.fullmatch(r'\S+')).all().all():
        raise ValueError('Invalid sample keep ID')
    return set(table.itertuples(index=False, name=None))


def load_panel(filename, chrom):
    # Read only this chromosome into memory, even for multi-million-variant panels.
    with pd.read_csv(filename, sep='\t', keep_default_na=False, chunksize=250_000,
                     dtype={'CHROM': 'int64', 'POS': 'int64', 'REF': str, 'ALT': str, 'RSID': str}) as chunks:
        table = pd.concat([chunk.loc[chunk.CHROM.eq(chrom)] for chunk in chunks], ignore_index=True)
    if table.empty:
        raise ValueError(f'No panel variants on chromosome {chrom}')
    valid = (table.POS.gt(0) & table.REF.isin(list('ACGT')) & table.ALT.isin(list('ACGT'))
             & table.REF.ne(table.ALT) & table.RSID.str.fullmatch(RSID.pattern))
    if not valid.all():
        raise ValueError('Invalid panel row')
    if table.duplicated(['POS', 'REF', 'ALT']).any() or table.RSID.duplicated().any():
        raise ValueError('Ambiguous panel; use prepare_panel.py')
    return dict(zip(table[['POS', 'REF', 'ALT']].itertuples(index=False, name=None), table.RSID))


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
    pd.DataFrame(pairs, columns=['#FID', 'IID']).to_csv(filename, sep='\t', index=False)


def line_count(filename):
    with filename.open() as handle:
        return sum(1 for _ in handle)


def validate_zip(filename):
    with zipfile.ZipFile(filename) as zipped:
        if '.zmetadata' not in zipped.namelist() or zipped.testzip() is not None:
            raise ValueError(f'Invalid LD ZIP {filename}')


def local_references(cfg, chrom):
    """Check local inputs before genotype processing; never download on workers."""
    references = [path(cfg['metadata']), path(cfg['panel']), path(cfg['genetic_map'], chrom)]
    if cfg['sample_keep']:
        references.append(path(cfg['sample_keep']))
    if cfg['estimator'] == 'block':
        references += [path(cfg['blocks'].replace('{ancestry}', group), chrom)
                       for group in cfg['ancestries']]
    missing = [str(f) for f in references if not f.is_file() or f.stat().st_size == 0]
    if missing:
        raise FileNotFoundError('Missing local references: ' + ', '.join(missing)
                                + '. Run bash scripts/prepare_references.sh on the login node first.')
    if cfg['estimator'] == 'block':
        for group in cfg['ancestries']:
            blockfile = path(cfg['blocks'].replace('{ancestry}', group), chrom)
            parse_blocks(blockfile.read_text(), chrom=chrom)
    return references


def compute_ld(bfile, output, temporary, cfg, chrom, group, threads):
    import magenpy as mgp
    import zarr
    if cfg['estimator'] == 'windowed':
        kwargs = {'cm_window_size': cfg['window_cm']}
    else:
        blockfile = path(cfg['blocks'].replace('{ancestry}', group), chrom)
        parse_blocks(blockfile.read_text(), chrom=chrom)
        kwargs = {'ld_blocks_file': str(blockfile)}
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
            ld_variants = int(matrix.shape[0])
            if cfg['estimator'] == 'block':
                matrix.set_store_attr('LD block file SHA256', sha256(blockfile))
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
    return ld_variants


def run(args):
    cfg = read_config(path(args.config))
    validate_config(cfg)
    if args.threads < 1:
        raise ValueError('Thread count must be positive')
    references = local_references(cfg, args.chrom)
    prefix = path(cfg['input_prefix'], args.chrom)
    pgen, psam = Path(str(prefix) + '.pgen'), Path(str(prefix) + '.psam')
    pvars = [Path(str(prefix) + ext) for ext in ('.pvar', '.pvar.gz', '.pvar.zst')]
    pvars = [p for p in pvars if p.is_file()]
    if len(pvars) != 1:
        raise ValueError(f'Expected exactly one PVAR: {pvars}')
    pvar = pvars[0]
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
            all_pairs = [pair for pairs in selected.values() for pair in pairs]
            if not all_pairs:
                raise ValueError('No selected samples intersect PSAM; check FID/IID conventions')
            write_keep(work / 'all.keep', all_pairs)
            panel = load_panel(path(cfg['panel']), args.chrom)
            summary['requested_panel_variants'] = len(panel)
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
                variants = variant_table(mapping).loc[read_bim(folder / 'qc.bim').SNP].reset_index()
                variants.insert(2, 'CHROM', args.chrom)
                variants.to_csv(folder / 'variants.tsv', sep='\t', index=False)
                stats['map'] = annotate(folder / 'qc.bim', path(cfg['genetic_map'], args.chrom), mapping)
                stats['ld_variants'] = compute_ld(folder / 'qc', folder / 'ld', folder / 'temp', cfg, args.chrom, group, args.threads)
                stats['variants_excluded_by_blocks'] = stats['variants'] - stats['ld_variants']
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
