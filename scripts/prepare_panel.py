#!/usr/bin/env python3
"""Stream GRCh38 dbSNP once, retaining unambiguous requested biallelic SNPs."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
import os
from pathlib import Path
import re

import pandas as pd

COLUMNS = ['CHROM', 'POS', 'REF', 'ALT', 'RSID']
RSID = re.compile(r'rs[0-9]+')


def chromosome(value):
    value = value.removeprefix('chr')
    if value.startswith('NC_'):
        value = value[3:].split('.')[0].lstrip('0')
    return int(value) if value.isdigit() and 1 <= int(value) <= 22 else None


def sha256(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def identity(path):
    path = Path(path).resolve()
    stat = path.stat()
    return {'path': str(path), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}


def resolve_panel(vcf, wanted):
    matches = {}
    ambiguous = set()
    opener = gzip.open if str(vcf).endswith(('.gz', '.bgz')) else open
    counts = Counter()
    with opener(vcf, 'rt') as handle:
        for line in handle:
            if line.startswith('#'):
                continue
            f = line.rstrip().split('\t', 5)
            if len(f) < 5:
                raise ValueError('Malformed VCF record')
            ids = wanted.intersection(f[2].split(';'))
            if not ids:
                continue
            counts['requested_records'] += 1
            chrom = chromosome(f[0])
            ref, alt = f[3].upper(), f[4].upper()
            # A multiallelic record is excluded, even if one ALT is a SNP.
            if chrom is None or ref not in 'ACGT' or alt not in 'ACGT' or len(ref) != 1 or len(alt) != 1 or ref == alt:
                continue
            key = (chrom, int(f[1]), ref, alt)
            if key[1] < 1:
                raise ValueError('Invalid VCF position')
            for rsid in ids:
                if rsid in matches and matches[rsid] != key:
                    ambiguous.add(rsid)
                matches[rsid] = key
    # Multiple requested rsIDs at the same allele-aware locus are ambiguous too.
    loci = Counter(key for rsid, key in matches.items() if rsid not in ambiguous)
    result = [(key, rsid) for rsid, key in matches.items()
              if rsid not in ambiguous and loci[key] == 1]
    counts.update(requested=len(wanted), resolved=len(result),
                  ambiguous_rsids=len(ambiguous),
                  unresolved_or_excluded=len(wanted) - len(result))
    return sorted(result), dict(counts)


def read_rsids(paths):
    wanted = set()
    for path in paths:
        for line in path.read_text().splitlines():
            value = line.strip()
            if not value or value.startswith('#'):
                continue
            if not RSID.fullmatch(value):
                raise ValueError(f'Expected one rsID per line in {path}: {value!r}')
            wanted.add(value)
    if not wanted:
        raise ValueError('Empty panel')
    return wanted


def verify_dbsnp(vcf):
    """Verify a local checksum on the compute node, never scan the VCF on login."""
    checksum_file = Path(str(vcf) + '.md5')
    if not checksum_file.is_file():
        return  # Independently supplied local VCFs may have been verified upstream.
    expected = re.search(r'\b[a-fA-F0-9]{32}\b', checksum_file.read_text())
    if expected is None:
        raise ValueError(f'Invalid checksum file: {checksum_file}')
    print(f'Verifying dbSNP MD5: {vcf}', flush=True)
    with vcf.open('rb') as handle:
        observed = hashlib.file_digest(handle, 'md5').hexdigest()
    if observed.lower() != expected.group().lower():
        raise ValueError(f'dbSNP checksum mismatch: {vcf}; replace/redownload the invalid file')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assembly', choices=['GRCh38'], required=True)
    parser.add_argument('--dbsnp-vcf', type=Path, required=True)
    parser.add_argument('--rsids', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path, default=Path('data/reference/panel.grch38.tsv'))
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Panel preparation requires a SLURM job. Use bash scripts/run_pipeline.sh on the login node.')
    before = [identity(p) for p in [args.dbsnp_vcf, *args.rsids]]
    verify_dbsnp(args.dbsnp_vcf)
    wanted = read_rsids(args.rsids)
    rows, summary = resolve_panel(args.dbsnp_vcf, wanted)
    if not rows:
        raise ValueError('No usable panel variants resolved')
    if before != [identity(p) for p in [args.dbsnp_vcf, *args.rsids]]:
        raise ValueError('Panel inputs changed during preparation; rerun with fixed inputs')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix('.tmp')
    pd.DataFrame.from_records([(*key, rsid) for key, rsid in rows], columns=COLUMNS).to_csv(
        temporary, sep='\t', index=False)
    temporary.replace(args.output)
    summary.update(assembly=args.assembly, dbsnp=identity(args.dbsnp_vcf),
                   lists=[dict(identity(p), sha256=sha256(p)) for p in args.rsids],
                   panel=identity(args.output), panel_sha256=sha256(args.output),
                   chromosomes=dict(Counter(key[0] for key, _ in rows)))
    receipt = args.output.with_suffix('.json.tmp')
    receipt.write_text(json.dumps(summary, indent=2) + '\n')
    receipt.replace(args.output.with_suffix('.json'))
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
