import copy
import pandas as pd
import numpy as np
import logging
import tacopy as pt
from tacopy.util._decorators import set_module
from tacopy.core.add_counts import add_relative_abundances
from tacopy.core.add_taxa import (
    calculate_prevalence,
    name_n_most_prevalent,
    name_n_most_abundant,
)
from tacopy.core.cluster import cluster_samples


@set_module("tacopy")
class Tacopy:
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
            df["sample_id"] = df.index + 1
            df["sample_id"] = "s" + df["sample_id"].astype("str")
        self.__dict__["samples"] = df

    def _ensure_taxon_id(self, df: pd.DataFrame):
        if not "taxon_id" in df.columns:
            df = df.reset_index(drop=True)
            df["taxon_id"] = df.index + 1
            df["taxon_id"] = "t" + df["taxon_id"].astype("str")
        self.__dict__["taxa"] = df

    def _extract_rank_names(self):
        taxa_cols = self.taxa.columns
        rank_names = [
            c
            for c in taxa_cols
            if c not in ["taxon", "taxon_id", "sequence", "prevalence"]
        ]
        return rank_names

    def __getattribute__(self, name):
        return super().__getattribute__(name)

    def __getitem__(self, key):
        return getattr(self, key)

    def __setattr__(self, name, value):
        if name == "samples":
            if not "sample_id" in value.columns:
                return self._ensure_sample_id(value)

            # Adapt counts table
            self.__dict__["counts"] = self.counts[
                self.counts.sample_id.isin(value.sample_id)
            ]
            # Adapt taxa table
            self.__dict__["taxa"] = self.taxa[
                self.taxa.taxon_id.isin(self.counts.taxon_id)
            ]

        if name == "taxa":
            if not "taxon_id" in value.columns:
                return self._ensure_taxon_id(value)
            # Adapt counts table
            self.__dict__["counts"] = self.counts[
                self.counts.taxon_id.isin(value.taxon_id)
            ]
            # Adapt samples table
            self.__dict__["samples"] = self.samples[
                self.samples.sample_id.isin(self.counts.sample_id)
            ]

        if name == "counts":
            if not "taxon_id" in value.columns:
                raise ValueError(
                    "The taxon identifier 'taxon_id' is missing from the counts table"
                )
            if not "sample_id" in value.columns:
                raise ValueError(
                    "The sample identifier 'sample_id' is missing from the counts table"
                )
            if not "count" in value.columns:
                raise ValueError("The count column is missing")

            # Adapt samples & taxa tables
            self.__dict__["samples"] = self.samples[
                self.samples.sample_id.isin(self.counts.sample_id)
            ]
            self.__dict__["taxa"] = self.taxa[
                self.taxa.taxon_id.isin(self.counts.taxon_id)
            ]

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

    def describe(self):
        """Return a summary of the object in dataframe format"""
        desc = pd.DataFrame(
            {"count": [self.samples.shape[0], self.taxa.shape[0], self.counts.shape[0]]}
        )
        desc.index = pd.Index(["samples", "taxa", "counts"])
        return desc

    def summary(self):
        """Return a summary of the object"""
        return f"Pyditaco(samples={self.samples.shape[0]}, taxa={self.taxa.shape[0]}, counts={self.counts.shape[0]})"

    def matrix(self, taxon="taxon_id", sample="sample_id", values="count") -> np.matrix:
        kwargs = dict()
        tt = copy.copy(self)

        if taxon not in tt.counts.columns:
            tt = pt.aggregate_taxa(tt, rank=taxon)
            tt.counts = tt.counts.merge(
                tt.taxa[["taxon_id", taxon]], how="left", on="taxon_id"
            )

        if sample not in tt.counts.columns:
            tt = pt.aggregate_samples(tt, sample)
            sample = "sample"  # gets renamed internally
            tt.counts = tt.counts.merge(
                tt.samples[["sample_id", "sample"]], how="left", on="sample_id"
            )

        if values.lower() in ["rel_abundance", "relab", "rel_ab"]:
            values = "rel_abundance"
            if not values in self.taxa.columns:
                tt = tt.add_rel_abundance(overwrite=True)

            matrix = (
                tt.counts.pivot_table(
                    index="taxon_id", columns="sample_id", values="count", **kwargs
                )
                .fillna(0)
                .astype(float)
            )

        elif values.lower() in ["count", "counts"]:
            matrix = (
                tt.counts.pivot_table(
                    index="taxon_id", columns="sample_id", values="count", **kwargs
                )
                .fillna(0)
                .astype(int)
            )

        else:
            raise ValueError(
                f"{values} not supported; please select from `counts` or `rel_abundance`"
            )

        if taxon != "taxon_id":
            matrix = matrix.rename(
                index={k: v for k, v in zip(tt.taxa.taxon_id, tt.taxa[taxon])}
            )
        if sample != "sample_id":
            matrix = matrix.rename(
                columns={k: v for k, v in zip(tt.samples.sample_id, tt.samples[sample])}
            )
        return matrix

    def add_rel_abundance(self, sample="sample_id", overwrite=False):
        if not "count" in self.counts.columns:
            raise ValueError("The count column is missing")
        if not sample in self.counts.columns:
            raise ValueError(
                f"The sample identifier '{sample}' is missing from the counts table"
            )

        taco = copy.copy(self)
        taco.counts = add_relative_abundances(taco.counts, sample, overwrite)
        return taco

    def add_mean_rel_abundance(self, sample="sample_id", taxon="taxon_id"):
        keep_rel_abund = True
        taco = copy.copy(self)
        if not "rel_abundance" in self.counts.columns:
            keep_rel_abund = False
            taco = taco.add_rel_abundance(sample=sample)
        mean_rel_abund = taco.counts.groupby(taxon).mean("rel_abundance")[
            "rel_abundance"
        ]
        mean_rel_abund.name = "mean_rel_abundance"
        taco.taxa = pd.merge(taco.taxa, mean_rel_abund, left_on=taxon, right_on=taxon)
        if not keep_rel_abund:
            taco.counts = taco.counts.drop(columns=["rel_abundance"])
        return taco

    def add_prevalence(self):
        taco = copy.copy(self)
        prevalence = calculate_prevalence(taco)
        taco.taxa = pd.merge(taco.taxa, prevalence, on="taxon_id")
        return taco

    def drop_empty_samples(self):
        taco = copy.copy(self)
        taco.samples = taco.samples[
            taco.samples["sample_id"].isin(taco.counts["sample_id"])
        ]
        logging.info(f"{self.samples.count() - taco.samples.count()} samples dropped.")
        return taco

    def add_taxon_name(
        self,
        method="mean_rel_abundance",
        n=12,
        rank="taxon",
    ):
        if "taxon_name" in self.taxa.columns:
            self.taxa = self.taxa.drop(columns=["taxon_name"])

        if not method in ["mean_rel_abundance", "prevalence"]:
            raise ValueError(f"Unknown method '{method}'")

        taco = copy.copy(self)
        if not rank == "taxon":
            taco = pt.aggregate_taxa(taco, rank)
            taco.taxa = taco.taxa.drop(columns="taxon_name")

        if method == "mean_rel_abundance":
            taco = name_n_most_abundant(taco, n)
        elif method == "prevalence":
            taco = name_n_most_prevalent(taco, n)

        return taco

    def cluster_samples(self, n=None, method="kmeans", *args, **kwargs):
        taco = copy.copy(self)
        if "cluster" in taco.samples.columns:
            taco.samples = taco.samples.drop(columns=["cluster"])

        if not n:
            n = int(np.round(np.log(taco.samples["sample_id"].count())))

        if method == "agglomerative":
            model = cluster_samples(
                taco, method=method, n_clusters=n, linkage="average", *args, **kwargs
            )
        elif method == "kmeans":
            model = cluster_samples(
                taco, method=method, n_clusters=n, random_state=42, *args, **kwargs
            )
        clusters = taco.drop_empty_samples().samples[["sample_id"]]

        clusters["cluster"] = np.char.add("cluster_", model.labels_.astype("str"))
        taco.samples = taco.samples.set_index("sample_id").join(
            clusters.set_index("sample_id")
        )
        return taco
    
    def everything(self):
        return self.counts\
            .merge(self.samples, on="sample_id", how="left")\
                .merge(self.taxa, on="taxon_id", how="left")

    def write(self, path: str):
        from tacopy.core.io import write_tacos

        write_tacos(self, path)
