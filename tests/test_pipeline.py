"""Synthetic correctness tests; no private data needed."""
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from prepare_panel import chromosome, resolve_panel
from run_chromosome import annotate, load_panel, match_pvar, samples, sample_ids, validate_zip, load_keep, write_keep
from prepare_genetic_map import normalize_map, map_statistics


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.addCleanup(self.tmp.cleanup)

    def write(self, name, text):
        p = self.root / name
        p.write_text(text)
        return p

    def test_chromosomes(self):
        self.assertEqual(chromosome('NC_000022.11'), 22)
        self.assertEqual(chromosome('chr10'), 10)
        self.assertIsNone(chromosome('NC_000023.11'))
        self.assertIsNone(chromosome('chr1_random'))

    def test_rsid_resolution_excludes_indels_multiallelic_and_ambiguity(self):
        vcf = self.write('db.vcf', '##reference=GRCh38\n#CHROM\tPOS\tID\tREF\tALT\n'
            'NC_000001.11\t10\trs1\tA\tG\n'
            'NC_000001.11\t20\trs2\tC\tT\n'
            'NC_000001.11\t30\trs2\tC\tG\n'
            'NC_000001.11\t40\trs3;rs4\tA\tC\n'
            'NC_000001.11\t50\trs5\tA\tAC\n'
            'NC_000001.11\t60\trs6\tA\tC,G\n')
        result, stats = resolve_panel(vcf, {f'rs{i}' for i in range(1, 8)})
        self.assertEqual(result, [((1, 10, 'A', 'G'), 'rs1')])
        self.assertEqual(stats['unresolved_or_excluded'], 6)

    def test_pvar_exact_alleles_and_quality(self):
        pvar = self.write('chr1.pvar', '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tINFO\n'
            'chr1\t10\tchr1:10:A:G\tA\tG\tR2=0.9\n'
            'chr1\t20\tchr1:20:G:A\tG\tA\tR2=0.99\n'
            'chr1\t30\tchr1:30:A:C\tA\tC\tR2=NaN\n'
            'chr1\t40\tchr1:40:A:C\tA\tC\tR2=0.4\n')
        panel = {(10, 'A', 'G'): 'rs1', (20, 'A', 'G'): 'rs2',
                 (30, 'A', 'C'): 'rs3', (40, 'A', 'C'): 'rs4'}
        found, stats = match_pvar(pvar, panel, 1, 'R2', 0.8)
        self.assertEqual(list(found), ['chr1:10:A:G'])
        self.assertEqual(stats['info_failed'], 2)

    def test_duplicate_source_records_excluded(self):
        pvar = self.write('chr1.pvar', '#CHROM POS ID REF ALT\n1 10 id1 A G\n1 10 id2 A G\n1 20 id3 C T\n')
        found, stats = match_pvar(pvar, {(10, 'A', 'G'): 'rs1', (20, 'C', 'T'): 'rs2'}, 1)
        self.assertEqual(list(found), ['id3'])
        self.assertEqual(stats['selected'], 1)

    def test_missing_info_fails(self):
        pvar = self.write('chr1.pvar', '#CHROM POS ID REF ALT\n1 10 id1 A G\n')
        with self.assertRaisesRegex(ValueError, 'no INFO'):
            match_pvar(pvar, {(10, 'A', 'G'): 'rs1'}, 1, 'R2')

    def test_samples_preserve_strings_and_intersect_pairs(self):
        psam = self.write('a.psam', '#FID IID SEX\n001 001 2\n002 002 1\n004 004 2\n')
        meta = self.write('labels.tsv', 'FID\tIID\tAncestry\n001\t001\tEUR\n002\t002\tAFR\n003\t003\tEUR\n004\t004\tOTH\n005\t005\tOTH\n')
        groups, missing = samples(meta, sample_ids(psam), ['EUR', 'AFR'], {('001', '001')})
        self.assertEqual(groups, {'EUR': [('001', '001')], 'AFR': []})
        self.assertEqual(missing, 1)

    def test_psam_without_fid(self):
        psam = self.write('a.psam', '#IID SEX\n001 2\n')
        self.assertEqual(sample_ids(psam), {('0', '001')})

    def test_duplicate_metadata_rejected(self):
        meta = self.write('labels.tsv', 'FID\tIID\tAncestry\n1\t1\tEUR\n2\t1\tAFR\n')
        with self.assertRaisesRegex(ValueError, 'Duplicate'):
            samples(meta, {('1', '1')}, ['EUR'])

    def test_cm_interpolation_and_endpoint_clamping_preserve_allele_order(self):
        genetic_map = self.write('map.tsv', 'position\tcm\n10\t1\n30\t3\n')
        bim = self.write('test.bim', '1 id1 0 5 G A\n1 id2 0 20 C T\n1 id3 0 40 A C\n')
        stats = annotate(bim, genetic_map,
                         {'id1': ((5, 'A', 'G'), 'rs1'), 'id2': ((20, 'T', 'C'), 'rs2'), 'id3': ((40, 'C', 'A'), 'rs3')})
        rows = [line.split() for line in bim.read_text().splitlines()]
        self.assertEqual([r[2] for r in rows], ['1', '2', '3'])
        self.assertEqual(rows[0][4:], ['G', 'A'])
        self.assertEqual(stats, {'below': 1, 'interpolated': 1, 'above': 1})

    def test_decreasing_map_rejected(self):
        genetic_map = self.write('map.tsv', 'position\tcm\n10\t3\n30\t1\n')
        bim = self.write('test.bim', '1 id1 0 20 A G\n')
        with self.assertRaisesRegex(RuntimeError, 'decreasing'):
            annotate(bim, genetic_map, {})

    def test_panel_duplicate_rsid_rejected(self):
        panel = self.write('panel.tsv', 'CHROM\tPOS\tREF\tALT\tRSID\n1\t10\tA\tG\trs1\n1\t20\tC\tG\trs1\n')
        with self.assertRaisesRegex(ValueError, 'Ambiguous'):
            load_panel(panel, 1)

    def test_reference_cli_and_load(self):
        vcf = self.write('db.vcf', '#CHROM\tPOS\tID\tREF\tALT\n1\t10\trs1\tA\tG\n')
        ids = self.write('ids', 'rs1\nrs2\n')
        output = self.root / 'panel.tsv'
        script = Path(__file__).resolve().parents[1] / 'scripts/prepare_panel.py'
        subprocess.run([sys.executable, str(script), '--assembly', 'GRCh38', '--dbsnp-vcf', str(vcf),
                        '--rsids', str(ids), '--output', str(output)], check=True, capture_output=True)
        self.assertEqual(load_panel(output, 1), {(10, 'A', 'G'): 'rs1'})
        self.assertEqual(json.loads(output.with_suffix('.json').read_text())['resolved'], 1)

    def test_map_normalization_roundtrip_and_validation(self):
        output = self.root / 'map.tsv'
        summary = normalize_map(io.BytesIO(b'1 marker1 0.5 10\n1 marker2 1.5 30\n'), output, 1)
        self.assertEqual(summary, map_statistics(output))
        self.assertEqual(summary['first'], {'position': 10, 'cm': 0.5})
        json.dumps(summary)
        with self.assertRaisesRegex(RuntimeError, 'strictly increasing'):
            normalize_map(io.BytesIO(b'1 marker1 0.5 10\n1 marker2 1.5 10\n'), output, 1)
        self.assertEqual(summary, map_statistics(output))  # Invalid input must not overwrite output.

    def test_keep_roundtrip_preserves_literal_na_and_leading_zeros(self):
        pairs = [('001', '001'), ('NA', 'NA')]
        filename = self.root / 'samples.keep'
        write_keep(filename, pairs)
        self.assertEqual(load_keep(filename), set(pairs))
        meta = self.write('labels.tsv', 'FID\tIID\tAncestry\n001\t001\tEUR\nNA\tNA\tEUR\n')
        selected, missing = samples(meta, set(pairs), ['EUR', 'AFR'], set())
        self.assertEqual(selected, {'EUR': [], 'AFR': []})
        self.assertEqual(missing, 0)

    def test_annotation_rejects_allele_mismatch_without_overwriting(self):
        genetic_map = self.write('map.tsv', 'position\tcm\n10\t1\n30\t3\n')
        bim = self.write('test.bim', '1 id1 0 20 A G\n')
        original = bim.read_bytes()
        with self.assertRaisesRegex(ValueError, 'allele/position mismatch'):
            annotate(bim, genetic_map, {'id1': ((20, 'A', 'C'), 'rs1')})
        self.assertEqual(bim.read_bytes(), original)

    def test_zip_validation(self):
        archive = self.root / 'ld.zip'
        with zipfile.ZipFile(archive, 'w') as handle:
            handle.writestr('data', b'123')
        with self.assertRaisesRegex(ValueError, 'Invalid LD ZIP'):
            validate_zip(archive)
        with zipfile.ZipFile(archive, 'a') as handle:
            handle.writestr('.zmetadata', '{}')
        validate_zip(archive)


if __name__ == '__main__':
    unittest.main()
