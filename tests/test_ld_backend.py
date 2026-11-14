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
            cfg = {'estimator': 'windowed', 'window_cm': 3, 'storage_dtype': 'int16'}
            compute_ld(prefix, root / 'ld', root / 'tmp', cfg, 22, 'EUR', 2)
            archive = root / 'ld/chr_22.zip'
            ld = mgp.LDMatrix.from_path(str(archive))
            observed = ld.to_csr(dtype='float64', return_symmetric=True).toarray()
            np.testing.assert_allclose(observed, np.corrcoef(dosage, rowvar=False), atol=0.0001)
            self.assertEqual(ld.get_store_attr('Ancestry'), 'EUR')
            self.assertEqual(list((root / 'ld').iterdir()), [archive])


if __name__ == '__main__':
    unittest.main()
