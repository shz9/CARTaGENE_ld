"""Verify reference preparation and local-only compute preflight."""
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from ld_blocks import parse_blocks
from prepare_blocks import prepare
from run_chromosome import local_references

ROOT = Path(__file__).resolve().parents[1]
BLOCK_SOURCES = json.loads((ROOT / 'config/pipeline.json').read_text())['block_ancestry_map']


class BlockTests(unittest.TestCase):
    def test_prior_pipeline_mappings(self):
        self.assertEqual(BLOCK_SOURCES, {'AFR': 'AFR', 'AMR': 'AFR', 'CSA': 'SAS',
                                        'EAS': 'EAS', 'EUR': 'EUR', 'MID': 'EUR'})
        cfg = json.loads((ROOT / 'config/pipeline.json').read_text())
        self.assertEqual((cfg['estimator'], cfg['storage_dtype'], cfg['blocks']),
                         ('block', 'int8', 'data/blocks/{ancestry}.bed'))

    def test_invalid_intervals_and_missing_chromosomes(self):
        with self.assertRaisesRegex(ValueError, 'Overlapping'):
            parse_blocks('chr\tstart\tstop\n1\t10\t30\n1\t20\t40\n')
        with self.assertRaisesRegex(ValueError, '22 autosomes'):
            parse_blocks('1\t10\t30\n', require_all=True)
        with self.assertRaisesRegex(ValueError, 'chromosome 22'):
            parse_blocks('1\t10\t30\n', chrom=22)

    def test_block_table_formats_and_fractional_coordinates(self):
        for content in ('# comment\nchr\tstart\tend\nchr2\t30\t40\nchr1\t10\t20\n',
                        'chr2 30 40\nchr1 10 20\n'):
            with self.subTest(content=content):
                table = parse_blocks(content)
                self.assertEqual(list(table.columns), ['chr', 'start', 'stop'])
                self.assertEqual(list(table.itertuples(index=False, name=None)),
                                 [(1, 10, 20), (2, 30, 40)])
        for content in ('1 1.5 20\n', '1 10 inf\n', 'chr\tstart\tfoo\n1\t10\t20\n'):
            with self.subTest(content=content), self.assertRaises(ValueError):
                parse_blocks(content)

    def test_prepare_uses_matching_sources_and_reuses_verified_cache(self):
        content = ('chr\tstart\tend\n' + ''.join(f'{c}\t10\t30\n' for c in range(1, 23))).encode()
        with tempfile.TemporaryDirectory() as tmp:
            output = Path(tmp)
            with patch('prepare_blocks.urllib.request.urlopen', side_effect=lambda *a, **k: io.BytesIO(content)) as download:
                summary = prepare(output, list(BLOCK_SOURCES), BLOCK_SOURCES)
                urls = [call.args[0].full_url for call in download.call_args_list]
                self.assertEqual(len(urls), 4)  # AFR/EUR proxies reuse downloaded sources.
                self.assertTrue(any('pyrho_SAS_' in url for url in urls))
                self.assertFalse(any('pyrho_AMR_' in url for url in urls))
            with patch('prepare_blocks.urllib.request.urlopen', side_effect=AssertionError('Unexpected download')):
                prepare(output, list(BLOCK_SOURCES), BLOCK_SOURCES)
            self.assertTrue(summary['ancestries']['MID']['proxy'])
            self.assertEqual((output / 'AMR.bed').read_bytes(), (output / 'AFR.bed').read_bytes())
            (output / 'CSA.bed').write_text('corrupt')
            with patch('prepare_blocks.urllib.request.urlopen', side_effect=lambda *a, **k: io.BytesIO(content)) as download:
                prepare(output, list(BLOCK_SOURCES), BLOCK_SOURCES)
                self.assertEqual(download.call_count, 1)
            modified = dict(BLOCK_SOURCES, MID='SAS')
            with patch('prepare_blocks.urllib.request.urlopen', side_effect=lambda *a, **k: io.BytesIO(content)) as download:
                summary = prepare(output, list(modified), modified)
                self.assertEqual(download.call_count, 1)
                self.assertIn('pyrho_SAS_', download.call_args.args[0].full_url)
                self.assertEqual(summary['ancestries']['MID']['source_ancestry'], 'SAS')


    def test_preflight_requires_local_blocks_for_each_ancestry(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            for filename in ('metadata', 'panel', 'map'):
                (root / filename).write_text('data')
            cfg = {'metadata': str(root / 'metadata'), 'panel': str(root / 'panel'),
                   'genetic_map': str(root / 'map'), 'sample_keep': None,
                   'estimator': 'block', 'blocks': str(root / '{ancestry}.bed'),
                   'ancestries': ['EUR', 'CSA']}
            with self.assertRaisesRegex(FileNotFoundError, 'login node'):
                local_references(cfg, 22)
            for group in cfg['ancestries']:
                (root / f'{group}.bed').write_text('chr\tstart\tstop\n22\t10\t30\n')
            self.assertEqual(len(local_references(cfg, 22)), 5)
            (root / 'CSA.bed').write_text('chr\tstart\tstop\n1\t10\t30\n')
            with self.assertRaisesRegex(ValueError, 'chromosome 22'):
                local_references(cfg, 22)

    def test_login_scripts_reject_slurm_invocation(self):
        env = dict(os.environ, SLURM_JOB_ID='123')
        commands = [['bash', str(ROOT / 'scripts/run_pipeline.sh')],
                    [sys.executable, str(ROOT / 'scripts/run_pipeline.py')],
                    ['bash', str(ROOT / 'scripts/prepare_references.sh')],
                    ['bash', str(ROOT / 'scripts/setup_env.sh')],
                    [sys.executable, str(ROOT / 'scripts/prepare_blocks.py')],
                    [sys.executable, str(ROOT / 'scripts/prepare_genetic_map.py')]]
        for command in commands:
            with self.subTest(command=command):
                result = subprocess.run(command, env=env, capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('login node', result.stderr)


if __name__ == '__main__':
    unittest.main()
