#!/usr/bin/env python3
"""Download and normalize ancestry-matched GRCh38 blocks on the LOGIN node."""
import argparse
import json
import os
from pathlib import Path
import urllib.request

from ld_blocks import parse_blocks, write_blocks
from prepare_panel import sha256

ROOT = Path(__file__).resolve().parent.parent
RAW = 'https://raw.githubusercontent.com/jmacdon/LDblocks_GRCh38/master/data/'


def prepare(output, populations, sources, method='pyrho', force=False):
    if not isinstance(sources, dict):
        raise ValueError('block_ancestry_map must be an ancestry-to-source object')
    for group in populations:
        if group not in sources:
            raise ValueError(f'Missing block_ancestry_map entry for {group}')
        if sources[group] not in {'AFR', 'EAS', 'EUR', 'SAS'}:
            raise ValueError(f'Unsupported source block ancestry for {group}: {sources[group]!r}')
    output.mkdir(parents=True, exist_ok=True)
    manifest = output / 'summary.json'
    existing = json.loads(manifest.read_text()) if manifest.is_file() else {}
    result = {'assembly': 'GRCh38', 'method': method, 'ancestries': {}}
    cache = {}
    for group in populations:
        source = sources[group]
        url = RAW + f'{method}_{source}_LD_blocks.bed'
        destination = output / f'{group}.bed'
        previous = existing.get('ancestries', {}).get(group, {})
        if (not force and destination.is_file() and previous.get('url') == url
                and previous.get('sha256') == sha256(destination)):
            blocks = parse_blocks(destination.read_text(), require_all=True)
        else:
            if source not in cache:
                request = urllib.request.Request(url, headers={'User-Agent': 'cartagene-ld-preparation'})
                with urllib.request.urlopen(request, timeout=60) as response:
                    cache[source] = parse_blocks(response.read().decode('utf-8-sig'), require_all=True)
            blocks = cache[source]
            write_blocks(destination, blocks)
        result['ancestries'][group] = {
            'source_ancestry': source, 'proxy': group != source, 'url': url,
            'sha256': sha256(destination), 'blocks': len(blocks),
            'chromosomes': blocks['chr'].value_counts().sort_index().to_dict()}
        print(f'{group} <- {source}: {len(blocks)} blocks, {destination}', flush=True)
    temporary = manifest.with_suffix('.json.tmp')
    temporary.write_text(json.dumps(result, indent=2) + '\n')
    temporary.replace(manifest)
    return result


def main():
    if os.environ.get('SLURM_JOB_ID'):
        raise SystemExit('Run reference downloads directly on the login node, outside a SLURM job.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT / 'config/pipeline.json')
    parser.add_argument('--output-dir', type=Path, default=ROOT / 'data/blocks')
    parser.add_argument('--populations', help='Comma-separated subset; defaults to config ancestries')
    parser.add_argument('--method', choices=['pyrho', 'decode'], default='pyrho')
    parser.add_argument('--force', action='store_true')
    args = parser.parse_args()
    config = json.loads(args.config.read_text())
    populations = args.populations.split(',') if args.populations else config['ancestries']
    if not populations or not set(populations) <= {'AFR', 'AMR', 'CSA', 'EAS', 'EUR', 'MID'}:
        parser.error('Unsupported ancestry; use AFR,AMR,CSA,EAS,EUR,MID')
    prepare(args.output_dir, populations, config['block_ancestry_map'], args.method, args.force)


if __name__ == '__main__':
    main()
