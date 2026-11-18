#!/usr/bin/env python3
"""Move legacy fingerprint folders aside and expose ancestry/chr_N.zip (no LD recomputation)."""
import argparse
import json
import re
from pathlib import Path


def migrate(output, work):
    output, work = output.resolve(), work.resolve()
    if not output.is_dir():
        return
    legacy = [p for p in output.iterdir() if p.is_dir() and re.fullmatch(r'[0-9a-f]{20}', p.name)]
    if not legacy:
        print('No legacy fingerprint folders found.')
        return
    planned, targets = [], set()
    for folder in legacy:
        if any(folder.glob('*.lock')):
            raise ValueError(f'Legacy jobs may still be running in {folder}; wait for completion before migrating')
        destination = work / folder.name
        if destination.exists():
            raise ValueError(f'Migration destination already exists: {destination}')
        for archive in folder.glob('chr*/*/ld/chr_*.zip'):
            target = output / archive.parent.parent.name / archive.name
            if target.exists() or target in targets:
                raise ValueError(f'Conflicting results for {target}; select the desired legacy run before migrating')
            targets.add(target)
            planned.append((destination / archive.relative_to(folder), target))
    work.mkdir(parents=True, exist_ok=True)
    if output.stat().st_dev != work.stat().st_dev:
        raise ValueError('Use a work directory on the same filesystem for lightweight rename-only migration')
    for folder in legacy:
        folder.rename(work / folder.name)
    published = {}
    for source, target in planned:
        target.parent.mkdir(parents=True, exist_ok=True)
        source.rename(target)
        stat = target.stat()
        published.setdefault(source.parents[2], {})[target.parent.name] = {
            'path': str(target), 'size': stat.st_size, 'mtime_ns': stat.st_mtime_ns}
        print(target)
    for chromosome, archives in published.items():
        (chromosome / 'published.json').write_text(json.dumps(archives, indent=2) + '\n')
    print(f'Moved private legacy diagnostics to {work}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path, default=Path('results'))
    parser.add_argument('--work-dir', type=Path, default=Path('data/work/legacy'))
    args = parser.parse_args()
    migrate(args.output_dir, args.work_dir)


if __name__ == '__main__':
    main()
