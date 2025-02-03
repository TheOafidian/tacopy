import tempfile
import unittest
import os
import tacopy as tp
import pandas as pd

class TestIOMethods(unittest.TestCase):

    def test_read_taco(self):
        taco = tp.read_tacos("tests/data/urt")
        self.assertIsInstance(taco, tp.Tacopy)
        self.assertEqual(len(taco.samples.index), 217)

    def test_write_taco(self):
        samples = pd.DataFrame(dict(sample_id=["s1","s2"],sample=["test1","test2"]))
        taxa = pd.DataFrame(dict(taxon_id=["t1","t2"], taxon=["tax1","tax2"]))
        counts = pd.DataFrame(dict(taxon_id=["t1","t1","t2"],sample_id=["s1","s2","s1"], count=[1,2,3]))
        tt = tp.Tacopy(samples,taxa,counts)
        
        with tempfile.TemporaryDirectory() as tmpdir:
            tt.write(tmpdir)
            assert os.path.exists(os.path.join(tmpdir, "samples.csv"))                    
            assert os.path.exists(os.path.join(tmpdir, "taxa.csv"))
            assert os.path.exists(os.path.join(tmpdir, "counts.csv"))

if __name__ == '__main__':
    unittest.main()