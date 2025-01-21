import pandas as pd
import numpy as np
from pyditacos.util._decorators import set_module

@set_module("pyditaco")
class Pyditaco():

    def __init__(self, samples, taxa, counts):
        if any([samples.empty, taxa.empty, counts.empty]):
            raise ValueError("One of the tables is empty")
        self._ensure_sample_id(samples)
        self._ensure_taxon_id(taxa)
        self.counts = counts
        self.rank_names = self._extract_rank_names()

    def _ensure_sample_id(self, df:pd.DataFrame):
        if not "sample_id" in df.columns:
            df = df.reset_index(drop=True)
            df = df.reset_index(names=['sample_id'])
            df['sample_id'] = 's' + df['sample_id'].astype('str')
        self.samples = df 


    def _ensure_taxon_id(self, df:pd.DataFrame):
        if not "taxon_id" in df.columns:
            df = df.reset_index(drop=True)
            df = df.reset_index(names=['taxon_id'])
            df['taxon_id'] = 's' + df['taxon_id'].astype('str')
        self.taxa = df 


    def _extract_rank_names(self):
        taxa_cols = self.taxa.columns
        rank_names = [c for c in taxa_cols if c not in ['taxon','taxon_id', 'sequence']]
        return rank_names

    def matrix(self, taxon="taxon_id", sample="sample_id", dtype=int) -> np.matrix:
        return self.counts.pivot(
                   index=taxon,columns=sample).fillna(0).astype(dtype)

    def write(self, path:str):
        from pyditacos.core.io import write_pyditacos
        write_pyditacos(self, path)