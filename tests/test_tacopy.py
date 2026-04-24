import unittest
import os
import tacopy as tp
import pandas as pd

class TestTacoMethods(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls._urt = tp.read_tacos("tests/data/urt")
        return super().setUpClass()

    def test_summary(self):
        self.assertEqual(self._urt.summary(), "Pyditaco(samples=217, taxa=1957, counts=7693)")

    def test_describe(self):
        outcome = pd.DataFrame({"count": [217, 1957, 7693]})
        outcome.index = ["samples","taxa","counts"]
        pd.testing.assert_frame_equal(self._urt.describe(), outcome)
   
    def test_matrix(self):
        matrix = self._urt.matrix() 
        self.assertEqual(matrix.shape, (1957, 214))       
  
    def test_matrix_taxa(self):
        matrix_kingdom = self._urt.matrix(taxon="kingdom")
        self.assertEqual(matrix_kingdom.shape, (2,214))
    
    def test_matrix_sample(self):
        matrix_sample = self._urt.matrix(sample="location")
        self.assertEqual(matrix_sample.shape, (1957, 2))
        
    def test_matrix_taxa_sample(self):
        matrix_sample_kingdom = self._urt.matrix(taxon="kingdom", sample="location")
        self.assertEqual(matrix_sample_kingdom.shape, (2,2))
        
    
    def test_matrix_taxa_sample_relative(self):    
        matrix_sample_kingdom_rel = self._urt.matrix(taxon="kingdom", sample="location", values="relab")
        self.assertEqual(matrix_sample_kingdom_rel.shape, (2,2))

    
    def test_add_rel_abundance(self):
        urt_rel = self._urt.add_rel_abundance()
        rel_all = urt_rel.counts.groupby("sample_id")["rel_abundance"].sum().sum()
        self.assertEqual(rel_all, urt_rel.samples.sample_id.count())


    def test_add_mean_rel_abundance(self):
        self.assertTrue((self._urt.add_mean_rel_abundance()\
            .taxa.mean_rel_abundance < 1).all())
    
    def test_add_prevalence(self):
        pass
    
    def test_add_taxon_name(self):
        urtg = self._urt.add_taxon_name(rank="genus")
    
    def test_cluster_samples(self):
        clusts=3
        tpc = self._urt.cluster_samples(n=clusts)
        # +1 because urt contains empty samples and these are not assigned a cluster
        self.assertEqual(len(tpc.samples.cluster.unique()), clusts+1)
    # todo: do for diff ranks
    def test_aggregate_taxa_genus(self):
        urt_agg = tp.aggregate_taxa(self._urt, rank="genus")
        self.assertNotIn("species", urt_agg.taxa.columns)
        self.assertNotIn("genus", urt_agg.counts.columns)
        self.assertEqual(len(urt_agg.taxa.taxon_id), 576)
        
    def test_aggregate_taxa_family(self):
        urt_agg = tp.aggregate_taxa(self._urt, rank="family")
        self.assertNotIn("genus", urt_agg.taxa.columns)
        self.assertNotIn("family", urt_agg.counts.columns)
        self.assertEqual(len(urt_agg.taxa.taxon_id), 235)

    def test_aggregate_taxa_kingdom(self):
        urt_agg = tp.aggregate_taxa(self._urt, rank="kingdom")
        self.assertEqual(len(urt_agg.taxa.taxon_id), 2) # urt has bacteria and unknown

    def test_drop_empty_samples(self):
        
        nsamples_after = self._urt\
            .drop_empty_samples()\
            .describe().loc["samples", "count"]
        self.assertEqual(self._urt.describe().loc["samples", "count"]-3, nsamples_after)


if __name__ == '__main__':
    unittest.main()