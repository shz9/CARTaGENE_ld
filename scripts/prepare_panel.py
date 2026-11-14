#!/usr/bin/env python3
"""Stream GRCh38 dbSNP once, retaining unambiguous requested biallelic SNPs."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
from pathlib import Path
import re

COLUMNS = ['CHROM', 'POS', 'REF', 'ALT', 'RSID']
RSID = re.compile(r'rs[0-9]+')


def chromosome(value):
    value = value.removeprefix('chr')
    if value.startswith('NC_'):
        value = value[3:].split('.')[0].lstrip('0')
    return int(value) if value.isdigit() and 1 <= int(value) <= 22 else None


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--assembly', choices=['GRCh38'], required=True)
    parser.add_argument('--dbsnp-vcf', type=Path, required=True)
    parser.add_argument('--rsids', type=Path, nargs='+', required=True)
    parser.add_argument('--output', type=Path, default=Path('data/reference/panel.grch38.tsv'))
    args = parser.parse_args()
    wanted = set()
    for path in args.rsids:
        for line in path.read_text().splitlines():
            value = line.strip()
            if not value or value.startswith('#'):
                continue
            if not RSID.fullmatch(value):
                raise ValueError(f'Expected one rsID per line in {path}: {value!r}')
            wanted.add(value)
    if not wanted:
        raise ValueError('Empty panel')
    rows, summary = resolve_panel(args.dbsnp_vcf, wanted)
    if not rows:
        raise ValueError('No usable panel variants resolved')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_suffix('.tmp')
    with temporary.open('w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(COLUMNS)
        for key, rsid in rows:
            writer.writerow([*key, rsid])
    temporary.replace(args.output)
    summary.update(assembly=args.assembly, dbsnp=identity(args.dbsnp_vcf),
                   lists=[dict(identity(p), sha256=sha256(p)) for p in args.rsids],
                   panel_sha256=sha256(args.output))
    args.output.with_suffix('.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
