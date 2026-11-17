"""Optional numerical integration check using the real magenpy native backend."""
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from run_chromosome import compute_ld


@unittest.skipUnless(importlib.util.find_spec('magenpy'), 'magenpy is not installed')
class NativeLDTest(unittest.TestCase):
    def test_signed_ld_matches_numpy_and_zip_reopens(self):
        import magenpy as mgp
        import numpy as np
        rng = np.random.default_rng(17)
        dosage = rng.binomial(2, 0.35, size=(40, 6))
        dosage[:, 1] = 2 - dosage[:, 0]  # Check negative as well as positive LD.
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            prefix = root / 'test'
            prefix.with_suffix('.fam').write_text(''.join(f'{i} {i} 0 0 0 -9\n' for i in range(40)))
            prefix.with_suffix('.bim').write_text(''.join(f'22 rs{j+1} {j/100} {100+j} A G\n' for j in range(6)))
            # Variant-major PLINK BED: 00=two A1, 10=one A1, 11=zero A1.
            encoding = {0: 3, 1: 2, 2: 0}
            with prefix.with_suffix('.bed').open('wb') as handle:
                handle.write(bytes([0x6c, 0x1b, 0x01]))
                for j in range(6):
                    for i in range(0, 40, 4):
                        handle.write(bytes([sum(encoding[int(dosage[i+k, j])] << (2*k) for k in range(4))]))
            cases = [('windowed', 'int16', 'EUR', None),
                     ('block', 'int8', 'EUR', 103), ('block', 'int8', 'AFR', 102)]
            for estimator, dtype, group, split in cases:
                with self.subTest(estimator=estimator, dtype=dtype, ancestry=group):
                    if split is not None:
                        (root / f'{group}.bed').write_text(
                            f'chr\tstart\tstop\n22\t90\t{split}\n22\t{split}\t110\n')
                    cfg = {'estimator': estimator, 'window_cm': 3, 'storage_dtype': dtype,
                           'blocks': str(root / '{ancestry}.bed')}
                    output = root / f'{estimator}-{group}'
                    compute_ld(prefix, output, root / 'tmp', cfg, 22, group, 2)
                    archive = output / 'chr_22.zip'
                    ld = mgp.LDMatrix.from_path(str(archive))
                    observed = ld.to_csr(dtype='float64', return_symmetric=True).toarray()
                    expected = np.corrcoef(dosage, rowvar=False)
                    if split is not None:
                        nleft = split - 100
                        expected[:nleft, nleft:] = 0
                        expected[nleft:, :nleft] = 0
                    np.testing.assert_allclose(observed, expected, atol=0.008 if dtype == 'int8' else 0.0001)
                    self.assertEqual(str(ld.stored_dtype), dtype)
                    self.assertEqual(ld.get_store_attr('Ancestry'), group)
                    self.assertEqual(list(output.iterdir()), [archive])



if __name__ == '__main__':
    unittest.main()
