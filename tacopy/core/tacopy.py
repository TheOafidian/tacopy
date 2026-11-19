import copy
import pandas as pd
import numpy as np
import logging
import re
from skbio.diversity import (
    alpha_diversity,
    get_alpha_diversity_metrics,
    beta_diversity,
    get_beta_diversity_metrics,
)
from skbio.stats.composition import clr
import tacopy as pt
from tacopy.util._decorators import set_module
from tacopy.core.add_counts import add_relative_abundances
from tacopy.core.add_taxa import (
    calculate_prevalence,
    name_n_most_prevalent,
    name_n_most_abundant,
)
from tacopy.core.cluster import cluster_samples, add_hier_clust_order

DEFAULT_ALPHA_METRICS = [
    "inv_simpson",
    "shannon",
    "simpson",
    "observed_features",
    "chao1",
    "ace",
    "pielou_e",
]


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

    def describe(self) -> pd.DataFrame:
        """Return a summary of the object in dataframe format.

        Returns
        -------
        pd.DataFrame
            A DataFrame with a 'count' column indexed by 'samples', 'taxa',
            and 'counts', reporting the number of rows in each respective table.
        """
        desc = pd.DataFrame(
            {"count": [self.samples.shape[0], self.taxa.shape[0], self.counts.shape[0]]}
        )
        desc.index = pd.Index(["samples", "taxa", "counts"])
        return desc

    def summary(self) -> str:
        """Return a short human-readable summary string of the object

        Returns
        -------
        str
            A string of the form ``Pyditaco(samples=N, taxa=M, counts=K)``
            reporting the number of rows in each respective table.
        """
        return f"Pyditaco(samples={self.samples.shape[0]}, taxa={self.taxa.shape[0]}, counts={self.counts.shape[0]})"

    def matrix(self, taxon="taxon_id", sample="sample_id", values="count") -> np.matrix:
        """
        Pivot the counts table into a sample-by-taxon abundance matrix.

        If the requested ``taxon`` or ``sample`` column is not already present
        in the counts table, taxa or samples are aggregated to that level before
        pivoting. The returned matrix is indexed by the taxon identifier and has
        samples as columns.

        Parameters
        ----------
        taxon : str, optional
            Column name to use as the row index of the matrix. If not
            ``'taxon_id'``, taxa are first aggregated to the requested rank
            (e.g. ``'genus'``, ``'family'``). Defaults to ``'taxon_id'``.
        sample : str, optional
            Column name to use as the column index of the matrix. If not
            ``'sample_id'``, samples are first aggregated to the requested
            level. Defaults to ``'sample_id'``.
        values : str, optional
            The values to fill the matrix with. Accepts ``'count'`` /
            ``'counts'`` (integer raw counts) or ``'rel_abundance'`` /
            ``'relab'`` / ``'rel_ab'`` (float relative abundances).
            Defaults to ``'count'``.

        Returns
        -------
        np.matrix
            A matrix of shape ``(n_taxa, n_samples)`` filled with either raw
            counts (``int``) or relative abundances (``float``). Missing
            combinations are filled with ``0``.

        Raises
        ------
        ValueError
            If ``values`` is not one of the supported options.
        """
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
                tt.add_rel_abundance(overwrite=True)

            matrix = (
                tt.counts.pivot_table(
                    index="taxon_id", columns="sample_id", values="rel_abundance", **kwargs
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

    def clr_matrix(
        self, taxon="taxon_id", sample="sample_id", values="count", pseudocount=1
    ) -> np.matrix:
        """
        Return a centred log-ratio (CLR) transformed sample-by-taxon matrix.

        A pseudocount is added before transformation to avoid log(0). The
        resulting DataFrame shares the index and column labels of the untransformed
        matrix produced by :meth:`matrix`.

        Parameters
        ----------
        taxon : str, optional
            Row index for the matrix. See :meth:`matrix` for details.
            Defaults to ``'taxon_id'``.
        sample : str, optional
            Column index for the matrix. See :meth:`matrix` for details.
            Defaults to ``'sample_id'``.
        values : str, optional
            Values to use before CLR transformation. See :meth:`matrix` for
            supported options. Defaults to ``'count'``.
        pseudocount : int or float, optional
            Value added to every cell before transformation to handle zeros.
            Defaults to ``1``.

        Returns
        -------
        pd.DataFrame
            CLR-transformed matrix with the same index and columns as
            :meth:`matrix`.
        """
        mat = self.matrix(taxon=taxon, sample=sample, values=values) + pseudocount
        clr_mat = pd.DataFrame(clr(mat))
        clr_mat.columns = mat.columns
        clr_mat.index = mat.index
        return clr_mat

    def add_sample_clustered(self):
        return add_hier_clust_order(self)

    def add_rel_abundance(self, sample:str="sample_id", overwrite:bool=False):
        """
        Add a relative abundance column to the counts table.

        Parameters
        ----------
        sample : str, optional
            Column in the counts table identifying samples, used to compute
            per-sample totals. Defaults to ``'sample_id'``.
        overwrite : bool, optional
            Passed to :func:`add_relative_abundances`. If ``True``, the counts
            table is modified in-place; if ``False`` (default) a copy is made.

        Returns
        -------
        Tacopy
            A new Tacopy object whose counts table has a ``'rel_abundance'``
            column added.

        Raises
        ------
        ValueError
            If the ``'count'`` column or the specified ``sample`` column is
            absent from the counts table.
        """
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
        """
        Add mean relative abundance per taxon to the taxa table.

        Relative abundance is computed if not already present and is removed
        again afterward unless it was present on entry. The mean is taken
        across all samples for each taxon and stored in a new
        ``'mean_rel_abundance'`` column in ``taco.taxa``.

        Parameters
        ----------
        sample : str, optional
            Column identifying samples in the counts table. Defaults to
            ``'sample_id'``.
        taxon : str, optional
            Column identifying taxa to group by. Defaults to ``'taxon_id'``.

        Returns
        -------
        Tacopy
            A new Tacopy object with ``'mean_rel_abundance'`` added to the
            taxa table.
        """
        keep_rel_abund = True
        taco = copy.copy(self)
        if not "rel_abundance" in self.counts.columns:
            keep_rel_abund = False
            taco = taco.add_rel_abundance(sample=sample)
        mean_rel_abund = taco.counts.groupby(taxon)["rel_abundance"].mean()
    
        mean_rel_abund.name = "mean_rel_abundance"
        taco.taxa = pd.merge(taco.taxa, mean_rel_abund, left_on=taxon, right_on=taxon)
        if not keep_rel_abund:
            taco.counts = taco.counts.drop(columns=["rel_abundance"])
        return taco

    def add_prevalence(self):
        """
        Add prevalence (proportion of samples in which a taxon is detected)
        to the taxa table.

        Returns
        -------
        Tacopy
            A new Tacopy object with a ``'prevalence'`` column added to the
            taxa table.
        """
        taco = copy.copy(self)
        prevalence = calculate_prevalence(taco)
        taco.taxa = pd.merge(taco.taxa, prevalence, on="taxon_id")
        return taco

    def drop_empty_samples(self):
        """
        Remove samples from the samples table that have no counts.

        Samples whose ``'sample_id'`` does not appear in the counts table are
        dropped. The number of dropped samples is logged at INFO level.

        Returns
        -------
        Tacopy
            A new Tacopy object with empty samples removed.
        """
        taco = copy.copy(self)
        taco.samples = taco.samples[
            taco.samples["sample_id"].isin(taco.counts["sample_id"])
        ]
        logging.info(f"{self.samples.count() - taco.samples.count()} samples dropped.")
        return taco

    def add_taxon_name(
        self,
        method="mean_rel_abundance",
        n=None,
        rank="taxon",
    ):
        """
        Assign human-readable taxon names based on abundance or prevalence.

        The top ``n`` taxa (by the chosen metric) receive their rank-level name;
        the remainder are labelled as ``'other'``. Any existing ``'taxon_name'``
        column is replaced. If ``rank`` is not ``'taxon'``, taxa are first
        aggregated to that rank before naming.

        Parameters
        ----------
        method : str, optional
            Ranking metric used to select the top ``n`` taxa. Either
            ``'mean_rel_abundance'`` (default) or ``'prevalence'``.
        n : int, optional
            Number of top taxa to name individually. Defaults to ``12``.
        rank : str, optional
            Taxonomic rank to aggregate to before naming. Defaults to
            ``'taxon'`` (no aggregation).

        Returns
        -------
        Tacopy
            A new Tacopy object with ``'taxon_name'`` set in the taxa table.

        Raises
        ------
        ValueError
            If ``method`` is not ``'mean_rel_abundance'`` or ``'prevalence'``.
        """
        if "taxon_name" in self.taxa.columns:
            self.taxa = self.taxa.drop(columns=["taxon_name"])

        if not method in ["mean_rel_abundance", "prevalence"]:
            raise ValueError(f"Unknown method '{method}'")

        taco = copy.copy(self)
        if not rank == "taxon":
            taco = pt.aggregate_taxa(taco, rank)
            taco.taxa = taco.taxa.drop(columns="taxon_name")
 
        if not n:
            n = len(taco.taxa.taxon_id)

        if method == "mean_rel_abundance":
            taco = name_n_most_abundant(taco, n)
        elif method == "prevalence":
            taco = name_n_most_prevalent(taco, n)

        return taco

    def cluster_samples(self, n=None, method="kmeans", *args, **kwargs):
        """
        Cluster samples and add cluster labels to the samples table.

        If ``n`` is not provided, the number of clusters is estimated as
        ``round(log(n_samples))``. Any existing ``'cluster'`` column is
        replaced. Cluster labels follow the pattern ``'cluster_0'``,
        ``'cluster_1'``, etc.

        Parameters
        ----------
        n : int or None, optional
            Number of clusters. If ``None`` (default), estimated automatically
            from sample count.
        method : str, optional
            Clustering algorithm. Either ``'kmeans'`` (default) or
            ``'agglomerative'`` (average linkage).
        *args, **kwargs
            Additional arguments forwarded to the underlying clustering
            function.

        Returns
        -------
        Tacopy
            A new Tacopy object with a ``'cluster'`` column added to the
            samples table.
        """
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

    def add_alpha(self, metric="shannon", *args, **kwargs):
        """
        Compute an alpha diversity metric and add it to the samples table.

        Parameters
        ----------
        metric : str, optional
            Alpha diversity metric to compute. Defaults to ``'shannon'``.
            Call :func:`get_alpha_diversity_metrics` for supported options.
        *args, **kwargs
            Additional arguments forwarded to :func:`get_alpha_diversity_metrics`.

        Returns
        -------
        Tacopy
            A new Tacopy object with the computed metric added as a column
            to the samples table.

        Raises
        ------
        ValueError
            If ``metric`` is not in the list of supported metrics.
        """
        taco = copy.copy(self)
        supported_metrics = get_alpha_diversity_metrics(*args, **kwargs)

        if metric not in supported_metrics:
            raise ValueError(
                f"Currently only the following metrics are supported:\n{', '.join(supported_metrics)}"
            )

        alpha = alpha_diversity(metric, taco.matrix().T)
        alpha.name = metric

        taco.samples = taco.samples.set_index("sample_id").join(alpha).reset_index()
        return taco

    def add_alphas(self, list_metrics=DEFAULT_ALPHA_METRICS):
        """
        Compute multiple alpha diversity metrics and add them to the samples table.

        Calls :meth:`add_alpha` iteratively for each metric in ``list_metrics``.

        Parameters
        ----------
        list_metrics : list of str, optional
            Alpha diversity metrics to compute. Defaults to
            ``DEFAULT_ALPHA_METRICS``.

        Returns
        -------
        Tacopy
            A new Tacopy object with one column per metric added to the
            samples table.
        """
        taco = copy.copy(self)

        for m in list_metrics:
            taco = taco.add_alpha(m)

        return taco

    def betas(self, metric="aitchison", pairs_col=None):
        """
        Compute pairwise beta diversity between all samples.

        For the ``'aitchison'`` metric, CLR-transformed counts are used with
        Euclidean distance. All other metrics are passed directly to
        :func:`beta_diversity`. Empty samples are dropped before computation.
        Sample metadata is joined onto both sides of each pair.

        Parameters
        ----------
        metric : str, optional
            Beta diversity metric. Defaults to ``'aitchison'``. Call
            :func:`get_beta_diversity_metrics` for other supported options.
        pairs_col : str or None, optional
            If provided, retain only pairs where both samples share the same
            value in this column (e.g. to compare within subjects). Shared
            metadata columns are deduplicated and the ``'_1'`` suffix is
            stripped. Defaults to ``None`` (all pairs returned).

        Returns
        -------
        pd.DataFrame
            Long-format DataFrame with one row per sample pair, containing a
            ``'beta'`` column, a ``'beta_metric'`` column, and all sample
            metadata columns for both members of each pair.

        Raises
        ------
        ValueError
            If ``metric`` is not in the list of supported metrics.
        """
        taco = copy.copy(self).drop_empty_samples()

        supported_metrics = get_beta_diversity_metrics()
        supported_metrics.append("aitchison")
        validate = True
        if metric not in supported_metrics:
            raise ValueError(
                f"Currently only the following metrics are supported:\n{', '.join(supported_metrics)}"
            )

        if metric == "aitchison":
            counts = taco.clr_matrix().T
            metric = "euclidean"
            validate = False

        else:
            counts = taco.matrix().T

        betas = beta_diversity(metric, counts, validate=validate)

        betas = betas.to_series().reset_index()
        betas.columns = ["sample_id_1", "sample_id_2", "beta"]
        if validate:
            betas["beta_metric"] = metric
        else:
            betas["beta_metric"] = "aitchison"
            
        samps1 = taco.samples.copy()
        samps1.columns = [s + "_1" for s in samps1.columns]

        samps2 = taco.samples.copy()
        samps2.columns = [s + "_2" for s in samps2.columns]

        betas = betas.merge(samps1).merge(samps2)
        if pairs_col:
            betas = (
                betas[betas[f"{pairs_col}_1"] == betas[f"{pairs_col}_2"]]
                .T.drop_duplicates()
                .T
            )
            varies_per_sample = betas.columns[
                betas.columns.str.endswith("_2")
            ].str.replace("_2$", "_1", regex=True)
            to_rename = {
                c: re.sub("_1$", "", c)
                for c in betas.columns
                if not c in varies_per_sample and c.endswith("_1")
            }
            betas = betas.rename(columns=to_rename)

        return betas.drop(columns=["sample_id_1", "sample_id_2"])

    def everything(self):
        """
        Return a fully joined DataFrame of counts, samples, and taxa.

        Performs a left join of the counts table onto both the samples and
        taxa tables.

        Returns
        -------
        pd.DataFrame
            A long-format DataFrame with one row per count, enriched with all
            sample and taxon metadata columns.
        """
        return self.counts\
            .merge(self.samples, on="sample_id", how="left")\
                .merge(self.taxa, on="taxon_id", how="left")

    def read(self, path:str):
        from tacopy.core.io import read_tacos
        return read_tacos(path)


    def write(self, path: str):
        """
        Write the Tacopy object to a tidytacos directory.

        Delegates to :func:`~tacopy.core.io.write_tacos`.

        Parameters
        ----------
        path : str
            Path to the output directory. Created if it does not exist.
        """
        from tacopy.core.io import write_tacos

        write_tacos(self, path)
