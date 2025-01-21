import pandas as pd
import numpy as np
from pyditacos.util._decorators import set_module
from pyditacos.core.add_counts import add_relative_abundances
from pyditacos.core.add_taxa import calculate_prevalence

@set_module("pyditaco")
class Pyditaco:
    def __init__(self, samples, taxa, counts):
        if any([samples.empty, taxa.empty, counts.empty]):
            raise ValueError("One of the tables is empty")
        self._ensure_sample_id(samples)
        self._ensure_taxon_id(taxa)
        self.__dict__["counts"] = counts
        self.__dict__["rank_names"] = self._extract_rank_names()

    def _ensure_sample_id(self, df: pd.DataFrame):
        if not "sample_id" in df.columns:
            df = df.reset_index(drop=True)
            df = df.reset_index(names=["sample_id"])
            df["sample_id"] = "s" + df["sample_id"].astype("str")
        self.__dict__["samples"] = df

    def _ensure_taxon_id(self, df: pd.DataFrame):
        if not "taxon_id" in df.columns:
            df = df.reset_index(drop=True)
            df = df.reset_index(names=["taxon_id"])
            df["taxon_id"] = "t" + df["taxon_id"].astype("str")
        self.__dict__["taxa"] = df

    def _extract_rank_names(self):
        taxa_cols = self.taxa.columns
        rank_names = [
            c for c in taxa_cols if c not in ["taxon", "taxon_id", "sequence"]
        ]
        return rank_names
    
    def __getattribute__(self, name):
        return super().__getattribute__(name)
    
    
    def __setattr__ (self, name, value):
        
        if name == "samples":
            if not "sample_id" in value.columns:
                return self._ensure_sample_id(value)
        
        if name == "taxa":
            if not "taxon_id" in value.columns:
                return self._ensure_taxon_id(value)
        
        if name == "counts":
            if not "taxon_id" in value.columns:
                raise ValueError("The taxon identifier 'taxon_id' is missing from the counts table")
            if not "sample_id" in value.columns:
                raise ValueError("The sample identifier 'sample_id' is missing from the counts table")
            if not "count" in value.columns:
                raise ValueError("The count column is missing")
        
        return super().__setattr__(name, value)
                

    def __repr__(self):
        """Return a string representation of the object"""
        samples = repr(self.samples)
        taxa = repr(self.taxa)
        counts = repr(self.counts)
        return "\n".join([samples, taxa, counts])

    def _repr_html_(self):
        """Return a HTML representation of the object"""
        samples = self.samples._repr_html_()
        taxa = self.taxa._repr_html_()
        counts = self.counts._repr_html_()
        return "\n".join([samples, taxa, counts])
    
    def summary(self):
        """Return a summary of the object"""
        return f"Pyditaco(samples={self.samples.shape[0]}, taxa={self.taxa.shape[0]}, counts={self.counts.shape[0]})"

    def matrix(self, taxon="taxon_id", sample="sample_id", dtype=int) -> np.matrix:
        return self.counts.pivot(index=taxon, columns=sample).fillna(0).astype(dtype)

    def add_rel_abundance(self, sample="sample_id"):
        if not "count" in self.counts.columns:
            raise ValueError("The count column is missing")
        if not sample in self.counts.columns:
            raise ValueError(
                f"The sample identifier '{sample}' is missing from the counts table"
            )

        self.counts = add_relative_abundances(self.counts, sample)
        return self

    def add_prevalence(self, rank=None):

        if rank:
            # TODO: aggregate taxa
            pass

        prevalence = calculate_prevalence(self)
        self.taxa = self.taxa.join(prevalence, on="taxon_id")
        return self

    def add_taxon_name_color(
        self,
        method="mean_rel_abundance",
        n=12,
        samples="sample_id",
        taxon="taxon_id",
        rank=None,
    ):
        if not method in ["mean_rel_abundance", "prevalence"]:
            raise ValueError(f"Unknown method '{method}'")

        if not rank:
            rank = self.rank_names[-1]

        if method == "mean_rel_abundance":
            pass
        elif method == "prevalence":
            pass
        
    def write(self, path: str):
        from pyditacos.core.io import write_pyditacos

        write_pyditacos(self, path)
