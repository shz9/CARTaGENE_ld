"""Synthetic login-node preparation and submission-gate integration tests."""
import argparse
from contextlib import ExitStack
import gzip
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
import run_pipeline as pipeline
import run_chromosome
from prepare_blocks import RAW
from prepare_panel import sha256


class SubmissionTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        cfg = json.loads((pipeline.ROOT / 'config/pipeline.json').read_text())
        cfg.update(input_prefix=str(self.root / 'chr{chrom}'),
                   metadata=str(self.root / 'metadata.tsv'), panel=str(self.root / 'panel.tsv'),
                   genetic_map=str(self.root / 'map{chrom}.tsv'), blocks=str(self.root / 'data/blocks/{ancestry}.bed'),
                   output_dir=str(self.root / 'results'), ancestries=['EUR'], min_samples=2)
        self.cfg = cfg
        self.config = self.root / 'config.json'
        self.config.write_text(json.dumps(cfg))
        self.write('metadata.tsv', 'FID\tIID\tAncestry\n001\t001\tEUR\n002\t002\tEUR\n')
        self.write('chr22.psam', '#FID IID SEX\n001 001 1\n002 002 2\n')
        (self.root / 'chr22.pgen').write_bytes(b'\x6c\x1b\x10test')
        self.write('chr22.pvar', '#CHROM POS ID REF ALT\n22 100 chr22:100:A:G A G\n')
        self.write('panel.tsv', 'CHROM\tPOS\tREF\tALT\tRSID\n22\t100\tA\tG\trs1\n22\t200\tC\tT\trs2\n')
        self.write('map22.tsv', 'position\tcm\n10\t0\n300\t1\n')
        blocks = self.write('data/blocks/EUR.bed', 'chr\tstart\tstop\n' + ''.join(f'{c}\t1\t1000\n' for c in range(1, 23)))
        self.write('data/blocks/summary.json', json.dumps({'ancestries': {'EUR': {
            'source_ancestry': 'EUR', 'url': RAW + 'pyrho_EUR_LD_blocks.bed', 'sha256': sha256(blocks)}}}))

    def write(self, filename, content):
        filename = self.root / filename
        filename.parent.mkdir(parents=True, exist_ok=True)
        filename.write_text(content)
        return filename

    def invoke(self, *args):
        with ExitStack() as stack:
            stack.enter_context(patch.object(pipeline, 'ROOT', self.root))
            stack.enter_context(patch.object(run_chromosome, 'ROOT', self.root))
            # read_config normally uses the repo's defaults; this fixture is a complete config.
            stack.enter_context(patch.object(pipeline, 'read_config', return_value=self.cfg))
            stack.enter_context(patch.object(sys, 'argv', ['run_pipeline.py', '--config', str(self.config), '--chromosomes', '22', *args]))
            stack.enter_context(patch.dict(os.environ, {}, clear=False))
            if 'SLURM_JOB_ID' in os.environ:
                stack.enter_context(patch.dict(os.environ, {'SLURM_JOB_ID': ''}))
            stack.enter_context(patch('run_pipeline.shutil.which', side_effect=lambda name: f'/bin/{name}'))
            stack.enter_context(patch('run_pipeline.subprocess.check_output', return_value='PLINK v2.00'))
            stack.enter_context(patch('urllib.request.urlopen', side_effect=AssertionError('Unexpected network access')))
            calls = stack.enter_context(patch.object(pipeline, 'command', return_value=subprocess.CompletedProcess([], 0, stdout='12345\n')))
            self.commands = calls
            previous = Path.cwd()
            try:
                pipeline.main()
            finally:
                os.chdir(previous)
            return calls

    def test_cached_inputs_submit_once_with_config_snapshot(self):
        calls = self.invoke()
        calls.assert_called_once()
        args = calls.call_args.args
        self.assertEqual(args[:3], ('sbatch', '--parsable', '--array=22%4'))
        self.assertEqual(json.loads(Path(args[-1]).read_text()), self.cfg)
        report = json.loads((self.root / 'logs/preflight.json').read_text())
        self.assertEqual((report['status'], report['job_id']), ('submitted', '12345'))

    def test_check_and_prepare_only_do_not_submit(self):
        for flag in ('--check-only', '--prepare-only'):
            with self.subTest(flag=flag):
                self.invoke(flag).assert_not_called()

    def test_one_missing_chromosome_prevents_entire_submission(self):
        with self.assertRaisesRegex(ValueError, 'chr21'):
            self.invoke('--chromosomes', '21-22')
        self.commands.assert_not_called()
        self.assertFalse((self.root / 'logs/preflight.json').exists())

    def test_invalid_pgen_prevents_submission(self):
        (self.root / 'chr22.pgen').write_bytes(b'bad')
        with self.assertRaisesRegex(ValueError, 'PGEN magic'):
            self.invoke()

    def test_stale_block_mapping_fails_check_only(self):
        self.cfg['block_ancestry_map']['EUR'] = 'AFR'
        with self.assertRaisesRegex(ValueError, 'Stale/mismatched blocks'):
            self.invoke('--check-only')

    def test_missing_panel_requests_lists_without_downloading(self):
        (self.root / 'panel.tsv').unlink()
        with self.assertRaisesRegex(ValueError, 'Supply --rsids'):
            self.invoke()

    def test_gzipped_pvar_header(self):
        pvar = self.root / 'chr22.pvar'
        with gzip.open(str(pvar) + '.gz', 'wt') as handle:
            handle.write(pvar.read_text())
        pvar.unlink()
        self.invoke('--check-only').assert_not_called()

    def test_chromosome_lists(self):
        self.assertEqual(pipeline.chromosomes('1-3,22,2'), [1, 2, 3, 22])
        for value in ('22-1', '0', '23', '1%4', 'all'):
            with self.subTest(value=value), self.assertRaises(argparse.ArgumentTypeError):
                pipeline.chromosomes(value)

    def test_download_checksum_failure_does_not_publish(self):
        destination = self.root / 'dbsnp.gz'
        partial = self.root / 'dbsnp.gz.part'
        partial.write_bytes(b'not the expected VCF')
        with patch('urllib.request.urlopen', return_value=io.BytesIO(b'0' * 32)), patch.object(pipeline, 'command'):
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                pipeline.download_dbsnp('https://example.org/dbsnp.gz', destination)
        self.assertFalse(destination.exists())


if __name__ == '__main__':
    unittest.main()
