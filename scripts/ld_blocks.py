"""Read, validate and write local GRCh38 LD-block tables."""
from io import StringIO

import numpy as np
import pandas as pd

COLUMNS = ['chr', 'start', 'stop']


def parse_blocks(content, chrom=None, require_all=False):
    try:
        table = pd.read_csv(StringIO(content), sep=r'\s+', header=None,
                            comment='#', dtype=str, keep_default_na=False)
    except pd.errors.EmptyDataError as exc:
        raise ValueError('Empty LD block file') from exc
    header = table.iloc[0].str.strip().str.lower().replace({'end': 'stop'})
    if {'chr', 'start'} <= set(header):
        if 'stop' not in set(header):
            raise ValueError('Missing end/stop column')
        table.columns = header
        table = table.iloc[1:][COLUMNS].copy()
    else:
        table = table.iloc[:, :3].copy()
        table.columns = COLUMNS
    if table.empty:
        raise ValueError('Empty LD block file')
    table['chr'] = table['chr'].str.lower().str.removeprefix('chr')
    table = table.apply(pd.to_numeric, errors='raise')
    if not np.isfinite(table.to_numpy()).all() or (table % 1 != 0).any().any():
        raise ValueError('LD block coordinates must be finite integers')
    table = table.astype('int64').sort_values(COLUMNS).reset_index(drop=True)
    if not (table['chr'].between(1, 22) & table.start.ge(0) & table.stop.gt(table.start)).all():
        raise ValueError('Invalid LD block interval')
    if (table['chr'].eq(table['chr'].shift()) & table.start.lt(table.stop.shift())).any():
        raise ValueError('Overlapping LD blocks')
    chromosomes = set(table['chr'])
    if require_all and chromosomes != set(range(1, 23)):
        raise ValueError('LD block reference must cover all 22 autosomes')
    if chrom is not None and chrom not in chromosomes:
        raise ValueError(f'No LD blocks for chromosome {chrom}')
    return table


def write_blocks(destination, blocks):
    temporary = destination.with_name(destination.name + '.tmp')
    blocks.to_csv(temporary, sep='\t', index=False)
    temporary.replace(destination)
