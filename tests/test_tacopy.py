import unittest
import os
import tacopy as tp
import pandas as pd

class TestTacoMethods(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._urt = tp.read_tacos("tests/data/urt")
        return super().setUpClass()

    def test_add_rel_abundance(self):
        urt_rel = self._urt.add_rel_abundance()
        rel_all = urt_rel.counts.groupby("sample_id")["rel_abundance"].sum().sum()
        self.assertEqual(rel_all, urt_rel.samples.sample_id.count())

    def test_drop_empty_samples(self):
        
        nsamples_after = self._urt\
            .drop_empty_samples()\
            .describe().loc["samples", "count"]
        self.assertEqual(self._urt.describe().loc["samples", "count"]-3, nsamples_after)

if __name__ == '__main__':
    unittest.main()