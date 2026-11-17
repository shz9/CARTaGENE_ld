#!/usr/bin/env python3
"""Annotate filtered BIM tables with rsIDs and interpolated GRCh38 cM."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

BIM_COLUMNS = ['CHR', 'SNP', 'CM', 'POS', 'A1', 'A2']


def read_bim(filename):
    table = pd.read_csv(filename, sep=r'\s+', header=None, keep_default_na=False,
                        dtype={0: str, 1: str, 3: 'int64', 4: str, 5: str})
    table.columns = BIM_COLUMNS
    return table


def validate_genetic_map(table):
    if list(table.columns) != ['position', 'cm'] or len(table) < 2:
        raise RuntimeError('Genetic map requires position/cm columns and at least two markers')
    if not np.isfinite(table.to_numpy()).all() or table.position.le(0).any():
        raise RuntimeError('Invalid genetic-map values')
    if (table.position % 1 != 0).any() or table.position.diff().le(0).any():
        raise RuntimeError('Genetic-map positions must be integers and strictly increasing')
    if table.cm.diff().lt(0).any():
        raise RuntimeError('Genetic-map cM values are decreasing')
    return table


def load_genetic_map(filename):
    table = validate_genetic_map(pd.read_csv(filename, sep=r'\s+'))
    return table.position.to_numpy(), table.cm.to_numpy()


def variant_table(mapping):
    return pd.DataFrame.from_records(
        [(old, rsid, *key) for old, (key, rsid) in mapping.items()],
        columns=['ORIGINAL_ID', 'RSID', 'POS', 'REF', 'ALT']).set_index('ORIGINAL_ID')


def annotate(bim, genetic_map, mapping):
    positions, cms = load_genetic_map(genetic_map)
    table = read_bim(bim)
    joined = table.join(variant_table(mapping), on='SNP', rsuffix='_expected', validate='one_to_one')
    alleles_match = ((joined.A1.eq(joined.REF) & joined.A2.eq(joined.ALT))
                     | (joined.A1.eq(joined.ALT) & joined.A2.eq(joined.REF)))
    if not (joined.POS.eq(joined.POS_expected) & alleles_match).all():
        raise ValueError('BIM allele/position mismatch or missing annotation')
    if table.POS.diff().lt(0).any() or joined.RSID.duplicated().any():
        raise ValueError('Unsorted/duplicate BIM')
    table['SNP'] = joined.RSID
    table['CM'] = np.interp(table.POS, positions, cms)
    below, above = table.POS.lt(positions[0]), table.POS.gt(positions[-1])
    exact = table.POS.isin(positions)
    counts = {'below': int(below.sum()), 'above': int(above.sum()), 'exact': int(exact.sum()),
              'interpolated': int((~(below | above | exact)).sum())}
    temporary = bim.with_suffix('.bim.tmp')
    table.to_csv(temporary, sep='\t', header=False, index=False, float_format='%.10g')
    temporary.replace(bim)
    return {key: value for key, value in counts.items() if value}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bim', type=Path, required=True)
    parser.add_argument('--genetic-map', type=Path, required=True)
    parser.add_argument('--rsid-map', type=Path, required=True)
    args = parser.parse_args()
    names = pd.read_csv(args.rsid_map, sep=r'\s+', header=None, dtype=str, keep_default_na=False)
    names.columns = ['old', 'new']
    if names.old.duplicated().any():
        raise ValueError('Duplicate source variant IDs')
    table = read_bim(args.bim)
    if not names.old.isin(table.SNP).all():
        raise ValueError('rsID map contains variants absent from BIM')
    renamed = table.SNP.map(names.set_index('old').new).fillna(table.SNP)
    mapping = {row.SNP: ((row.POS, row.A1, row.A2), rsid)
               for row, rsid in zip(table.itertuples(index=False), renamed)}
    print(annotate(args.bim, args.genetic_map, mapping))


if __name__ == '__main__':
    main()
