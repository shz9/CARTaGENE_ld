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
import zipfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import run_pipeline as pipeline
import run_chromosome
from prepare_blocks import RAW
from prepare_panel import identity, sha256, verify_dbsnp, PANEL_POLICY


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
        self.write('scripts/prepare_panel.py', '# synthetic code fingerprint\n')
        self.write('metadata.tsv', 'FID\tIID\tAncestry\n001\t001\tEUR\n002\t002\tEUR\n')
        self.write('chr22.psam', '#FID IID SEX\n001 001 1\n002 002 2\n')
        (self.root / 'chr22.pgen').write_bytes(b'\x6c\x1b\x10test')
        self.write('chr22.pvar', '#CHROM POS ID REF ALT\n22 100 chr22:100:A:G A G\n')
        self.write('panel.tsv', 'CHROM\tPOS\tREF\tALT\tRSID\n22\t100\tA\tG\trs1\n22\t200\tC\tT\trs2\n')
        self.write('map22.tsv', 'position\tcm\n10\t0\n300\t1\n')
        ids = self.write('metadata/rsids.txt', 'rs1\nrs2\n')
        self.write('panel.json', json.dumps({
            'assembly': 'GRCh38', 'variant_policy': PANEL_POLICY, 'lists': [dict(identity(ids), sha256=sha256(ids))],
            'panel': identity(self.root / 'panel.tsv'), 'panel_sha256': sha256(self.root / 'panel.tsv')}))
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

    def legacy_result(self, chrom=22):
        folder = self.root / ('results/' + 'a' * 20) / f'chr{chrom}'
        archive = folder / f'EUR/ld/chr_{chrom}.zip'
        archive.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(archive, 'w') as zipped:
            zipped.writestr('.zmetadata', '{}')
        references = [self.root / p for p in ('metadata.tsv', 'panel.tsv', f'map{chrom}.tsv', 'data/blocks/EUR.bed')]
        old_config = {k: v for k, v in self.cfg.items() if k != 'work_dir'}
        summary = {'fingerprint': {'chrom': chrom, 'config': old_config,
                   'inputs': [identity(self.root / f'chr{chrom}.{ext}') for ext in ('pgen', 'pvar', 'psam')],
                   'references': {str(p): sha256(p) for p in references}},
                   'ancestries': {'EUR': {'status': 'complete'}}}
        (folder / 'summary.json').write_text(json.dumps(summary))
        return archive

    def test_launcher_migrates_and_does_not_recompute_completed_outputs(self):
        self.legacy_result()
        self.invoke().assert_not_called()
        self.assertTrue((self.root / 'results/EUR/chr_22.zip').is_file())
        self.invoke().assert_not_called()  # Recognized on subsequent runs as well.

    def test_changed_analysis_settings_do_not_reuse_migrated_results(self):
        self.legacy_result()
        self.cfg['mac'] = 10
        self.invoke().assert_called_once()

    def test_changed_inputs_do_not_reuse_migrated_results(self):
        self.legacy_result()
        with (self.root / 'chr22.pgen').open('ab') as handle:
            handle.write(b'changed')
        self.invoke().assert_called_once()

    def test_only_missing_chromosomes_are_submitted(self):
        self.legacy_result()
        for ext in ('psam', 'pgen', 'pvar'):
            source = self.root / f'chr22.{ext}'
            target = self.root / f'chr21.{ext}'
            target.write_bytes(source.read_bytes().replace(b'22', b'21') if ext == 'pvar' else source.read_bytes())
        self.write('map21.tsv', 'position\tcm\n10\t0\n300\t1\n')
        calls = self.invoke('--chromosomes', '21-22')
        self.assertIn('--array=21%4', calls.call_args.args)

    def test_organize_only_migrates_without_preparation_or_submission(self):
        self.legacy_result()
        (self.root / 'panel.tsv').unlink()
        self.invoke('--organize-only').assert_not_called()
        self.assertTrue((self.root / 'results/EUR/chr_22.zip').is_file())

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
        self.cfg['rsids'] = []
        with self.assertRaisesRegex(ValueError, 'Supply --rsids'):
            self.invoke()

    def test_configured_list_changes_invalidate_panel(self):
        (self.root / 'metadata/rsids.txt').write_text('rs1\nrs3\n')
        with self.assertRaisesRegex(ValueError, 'missing or stale'):
            self.invoke('--check-only')
        self.commands.assert_not_called()

    def test_explicit_rsids_override_configured_list(self):
        override = self.write('alternate.txt', 'rs1\nrs2\n')
        self.write('panel.json', json.dumps({
            'assembly': 'GRCh38', 'variant_policy': PANEL_POLICY, 'lists': [dict(identity(override), sha256=sha256(override))],
            'panel': identity(self.root / 'panel.tsv'), 'panel_sha256': sha256(self.root / 'panel.tsv')}))
        (self.root / 'metadata/rsids.txt').unlink()
        self.invoke('--check-only', '--rsids', str(override)).assert_not_called()

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

    def test_download_defers_checksum_scan_to_compute(self):
        destination = self.root / 'dbsnp.gz'
        (self.root / 'dbsnp.gz.part').write_bytes(b'not the expected VCF')
        with patch('urllib.request.urlopen', return_value=io.BytesIO(b'0' * 32)), patch.object(pipeline, 'command'):
            pipeline.download_dbsnp('https://example.org/dbsnp.gz', destination)
        self.assertTrue(destination.is_file())
        self.assertTrue(Path(str(destination) + '.md5').is_file())
        with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
            verify_dbsnp(destination)

    def test_missing_panel_submits_only_preparation_then_exits(self):
        (self.root / 'panel.tsv').unlink()
        vcf = self.write('dbsnp.vcf', 'local reference')
        with patch('prepare_panel.read_rsids', side_effect=AssertionError('Heavy rsID loading on login')):
            calls = self.invoke('--dbsnp-vcf', str(vcf))
        calls.assert_called_once()
        args = calls.call_args.args
        self.assertEqual(args[:2], ('sbatch', '--parsable'))
        self.assertEqual(Path(args[2]).name, 'prepare_panel.sbatch')
        self.assertEqual(args[3:5], (vcf, self.root / 'panel.tsv'))
        report = json.loads((self.root / 'logs/preflight.json').read_text())
        self.assertEqual((report['status'], report['job_id']), ('awaiting_panel', '12345'))

    def test_prepare_only_can_submit_panel_job(self):
        (self.root / 'panel.tsv').unlink()
        vcf = self.write('dbsnp.vcf', 'local reference')
        calls = self.invoke('--prepare-only', '--dbsnp-vcf', str(vcf))
        calls.assert_called_once()
        self.assertEqual(Path(calls.call_args.args[2]).name, 'prepare_panel.sbatch')

    def test_active_panel_job_is_not_submitted_again(self):
        (self.root / 'panel.tsv').unlink()
        vcf = self.write('dbsnp.vcf', 'local reference')
        self.invoke('--dbsnp-vcf', str(vcf))
        with patch.object(pipeline, 'panel_job_active', return_value=True):
            self.invoke('--dbsnp-vcf', str(vcf)).assert_not_called()

    def test_completed_or_failed_panel_job_can_be_retried(self):
        (self.root / 'panel.tsv').unlink()
        vcf = self.write('dbsnp.vcf', 'local reference')
        self.invoke('--dbsnp-vcf', str(vcf))
        with patch.object(pipeline, 'panel_job_active', return_value=False):
            self.invoke('--dbsnp-vcf', str(vcf)).assert_called_once()

    def test_panel_job_status_handles_completed_ids(self):
        with patch('run_pipeline.subprocess.run', return_value=subprocess.CompletedProcess([], 0, stdout='54321\n')):
            self.assertFalse(pipeline.panel_job_active('12345'))
            self.assertTrue(pipeline.panel_job_active('54321'))

    def test_old_snp_only_panel_schedules_rebuild(self):
        receipt = self.root / 'panel.json'
        report = json.loads(receipt.read_text())
        report.pop('variant_policy')
        receipt.write_text(json.dumps(report))
        vcf = self.write('dbsnp.vcf', 'local reference')
        calls = self.invoke('--dbsnp-vcf', str(vcf))
        calls.assert_called_once()
        self.assertEqual(Path(calls.call_args.args[2]).name, 'prepare_panel.sbatch')

    def test_check_only_never_submits_missing_panel(self):
        (self.root / 'panel.tsv').unlink()
        with self.assertRaisesRegex(ValueError, 'missing or stale'):
            self.invoke('--check-only')
        self.commands.assert_not_called()


if __name__ == '__main__':
    unittest.main()
