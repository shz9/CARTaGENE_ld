"""Final-output publication, legacy migration, and Zenodo package layout."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tarfile
import tempfile
import unittest
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_chromosome import publish_results
from migrate_outputs import migrate

ROOT = Path(__file__).resolve().parents[1]


class OutputTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)

    def archive(self, filename):
        filename.parent.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(filename, 'w') as zipped:
            zipped.writestr('.zmetadata', '{}')
            zipped.writestr('matrix/data/0', 'synthetic chunk')
        return filename

    def test_only_ld_archives_are_published(self):
        work, output = self.root / 'scratch', self.root / 'results'
        self.archive(work / 'EUR/ld/chr_1.zip')
        (work / 'EUR/qc.fam').write_text('private sample ID')
        metadata = self.root / 'work/chr_1/fingerprint'
        metadata.parent.mkdir(parents=True)
        summary = {'ancestries': {'EUR': {'status': 'complete'}}}
        publish_results(work, output, metadata, 1, 'fingerprint', summary, {})
        self.assertEqual([p.relative_to(output).as_posix() for p in output.rglob('*') if p.is_file()], ['EUR/chr_1.zip'])
        self.assertTrue((metadata / 'EUR/qc.fam').exists())
        self.assertFalse((metadata / 'EUR/ld').exists())
        previous = json.loads((metadata.parent / 'current.json').read_text())
        summary['ancestries']['EUR']['status'] = 'skipped'
        publish_results(work, output, metadata, 1, 'fingerprint', summary, previous)
        self.assertFalse((output / 'EUR/chr_1.zip').exists())

    def test_legacy_migration_preserves_diagnostics(self):
        output = self.root / 'results'
        old = output / ('a' * 20)
        self.archive(old / 'chr1/EUR/ld/chr_1.zip')
        (old / 'chr1/summary.json').write_text('{}')
        work = self.root / 'work'
        migrate(output, work)
        self.assertTrue((output / 'EUR/chr_1.zip').exists())
        self.assertFalse(old.exists())
        self.assertTrue((work / ('a' * 20) / 'chr1/summary.json').exists())

    def test_migration_conflicts_do_not_move_anything(self):
        output = self.root / 'results'
        for token in ('a' * 20, 'b' * 20):
            self.archive(output / token / 'chr1/EUR/ld/chr_1.zip')
        with self.assertRaisesRegex(ValueError, 'Conflicting'):
            migrate(output, self.root / 'work')
        self.assertTrue((output / ('a' * 20)).exists())

    @unittest.skipUnless(all(shutil.which(c) for c in ('tar', 'unzip', 'sha256sum')), 'Packaging tools unavailable')
    def test_zenodo_package_contains_zarr_not_private_files(self):
        source, output, scratch = self.root / 'ld', self.root / 'zenodo', self.root / 'scratch'
        scratch.mkdir()
        self.archive(source / 'EUR/chr_1.zip')
        self.archive(source / 'EUR/chr_2.zip')
        (source / 'EUR/qc.fam').write_text('private sample ID')
        env = dict(os.environ, TMPDIR=str(scratch))
        env.pop('SLURM_JOB_ID', None)
        env.pop('SLURM_TMPDIR', None)
        env.pop('INPUT_DIR', None)
        env.pop('OUTPUT_DIR', None)
        subprocess.run(['bash', str(ROOT / 'scripts/prepare_for_zenodo.sh'), str(source), str(output)],
                       env=env, check=True, capture_output=True)
        with tarfile.open(output / 'EUR.tar.gz') as archive:
            names = archive.getnames()
        self.assertIn('EUR/chr_1/.zmetadata', names)
        self.assertIn('EUR/chr_2/matrix/data/0', names)
        self.assertFalse(any('qc.fam' in name or name.endswith('.zip') for name in names))
        self.assertTrue((output / 'EUR.tar.gz.sha256').exists())
        self.assertEqual(list(scratch.iterdir()), [])


if __name__ == '__main__':
    unittest.main()
