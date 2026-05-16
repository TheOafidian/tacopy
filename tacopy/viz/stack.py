import numpy as np
import pandas as pd
import tacopy as pt
from tacopy import Tacopy
from tacopy.core.decorators import ensure_taco_columns
import matplotlib.pyplot as plt
import seaborn as sns
import seaborn.objects as so

def _auto_figsize(n_samples: int, n_taxa: int) -> tuple[float, float]:
    """Estimate a good figsize based on data dimensions."""
    width = max(8, n_samples * 0.05)
    height = max(5, n_taxa * 0.3 + 3)
    return (width, height)

@ensure_taco_columns(
        ("counts", "rel_abundance", "add_rel_abundance"),
        ("samples", "cluster", "cluster_samples"),
        ("taxa", "taxon_name", "add_taxon_name"),
        ("taxa", "prevalence", "add_prevalence"),
)
def tacoplot_stack(taco: Tacopy, x="sample", **kwargs):
    if x not in taco.samples.columns:
        raise ValueError(
            f"The sample identifier '{x}' is missing from the samples table"
        )
    
    if taco.taxa.taxon_name.nunique() > 20:
        taco = taco.add_taxon_name()

    counts = taco.counts.merge(taco.samples, on="sample_id")

    counts = counts.sort_values(["cluster", x])
    ordered_samples = counts[x].drop_duplicates().tolist()
    counts = counts.set_index("taxon_id").join(taco.taxa.set_index("taxon_id"))

    agg_counts = (
        counts.groupby([x, "taxon_name"])["rel_abundance"].sum().reset_index()
    )

    agg_counts["rel_abundance"] = 100 * agg_counts.rel_abundance

    ##
    # Order taxa: most abundant at bottom, Other at top
    mean_abundance = (
        agg_counts.groupby("taxon_name")["rel_abundance"].mean().sort_values(ascending=False)
    )
    taxon_order = mean_abundance.index.tolist()
    if "Other" in taxon_order:
        taxon_order.remove("Other")
        taxon_order.append("Other")  # moves to top of stack

    agg_counts["taxon_name"] = pd.Categorical(
        agg_counts["taxon_name"], categories=taxon_order, ordered=True
    )
    ##

    agg_counts[x] = pd.Categorical(
        agg_counts[x], categories=ordered_samples, ordered=True
    )
    #agg_counts["taxon_name"] = pd.Categorical(agg_counts["taxon_name"])
    agg_counts = agg_counts.sort_values(x)

    n_samples = agg_counts[x].nunique()
    n_taxa = agg_counts["taxon_name"].nunique()

    figsize = kwargs.get("figsize", _auto_figsize(n_samples, n_taxa))
    dpi = kwargs.get("dpi", 400)

    figure = (agg_counts
     .pivot(index="sample", columns="taxon_name", values="rel_abundance")
     .plot.bar(stacked=True)
    )

    return figure



@ensure_taco_columns(
        ("counts", "rel_abundance", "add_rel_abundance"),
        ("samples", "cluster", "cluster_samples"),
        ("taxa", "taxon_name", "add_taxon_name"),
        ("taxa", "prevalence", "add_prevalence"),
)
def tacoplot_stack2(taco: Tacopy, x="sample", how="bar", *args, **kwargs):
    
    if x not in taco.samples.columns:
        raise ValueError(
            f"The sample identifier '{x}' is missing from the samples table"
        )
    
    if taco.taxa.taxon_name.nunique() > 20:
        taco = taco.add_taxon_name()

    counts = taco.counts.merge(taco.samples, on="sample_id")

    counts = counts.sort_values(["cluster", "sample_id"])
    ordered_samples = counts["sample_id"].drop_duplicates().tolist()
    counts = counts.set_index("taxon_id").join(taco.taxa.set_index("taxon_id"))

    agg_counts = (
        counts.groupby(["sample_id", "taxon_name"])["rel_abundance"].sum().reset_index()
    )

    agg_counts["rel_abundance"] = 100 * agg_counts.rel_abundance

    ##
    # Order taxa: most abundant at bottom, Other at top
    mean_abundance = (
        agg_counts.groupby("taxon_name")["rel_abundance"]
        .mean().sort_values(ascending=False)
    )
    taxon_order = mean_abundance.index.tolist()
    if "Other" in taxon_order:
        taxon_order.remove("Other")
        taxon_order.append("Other")  # moves to top of stack

    agg_counts["taxon_name"] = pd.Categorical(
        agg_counts["taxon_name"], categories=taxon_order, ordered=True
    )
    ##

    agg_counts["sample_id"] = pd.Categorical(
        agg_counts["sample_id"], categories=ordered_samples, ordered=True
    )
    #agg_counts["taxon_name"] = pd.Categorical(agg_counts["taxon_name"])
    agg_counts = agg_counts.sort_values("sample_id")

    n_samples = agg_counts["sample_id"].nunique()
    n_taxa = agg_counts["taxon_name"].nunique()

    figsize = kwargs.get("figsize", _auto_figsize(n_samples, n_taxa))
    dpi = kwargs.get("dpi", 400)

    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    plot_options = []
    if how == "bar":
        plot_options.append(so.Bar(edgewidth=0))
    else:
        plot_options.append(so.Area(alpha=0.7))
    plot_options.append(so.Stack())

    figure = (
        so.Plot(agg_counts, x="sample_id", y="rel_abundance", color="taxon_name")
        .add(*plot_options)
        .label(x="Sample", y="Relative abundance (%)")
        .scale(color="Paired")
        .on(ax)
    )
    ax.tick_params(axis="x", rotation=90)
    figure = figure.plot()
    plt.tight_layout(pad=1.2)

    return figure
