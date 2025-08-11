import numpy as np
from tacopy import Tacopy
import tacopy as pt
import matplotlib.pyplot as plt
import seaborn as sns
import seaborn.objects as so


def tacoplot_stack(taco: Tacopy, x="sample", how="bar", *args, **kwargs):
    if x not in taco.samples.columns:
        raise ValueError(
            f"The sample identifier '{x}' is missing from the samples table"
        )

    if "rank" in kwargs:
        taco = pt.aggregate_taxa(rank=kwargs["rank"])

    if "rel_abundance" not in taco.counts.columns:
        taco = taco.add_rel_abundance()

    if "prevalence" not in taco.taxa.columns:
        taco = taco.add_prevalence()

    if "figsize" not in kwargs:
        kwargs["figsize"] = (10, 10)
    if "dpi" not in kwargs:
        kwargs["dpi"] = 300

    fig, ax = plt.subplots(figsize=kwargs["figsize"], dpi=kwargs["dpi"])
    
    # sort by clustering
    if not "cluster" in taco.samples.columns:
        n = int(np.round(np.log(taco.samples["sample_id"].count())))
        if "n_clusters" in kwargs:
            n = kwargs.get("n_clusters")
    
        taco = taco.cluster_samples(n=n)
    counts = taco.counts.merge(taco.samples, left_on="sample_id", right_on="sample_id")
    counts["sample_id"] = counts["cluster"] + counts["sample_id"]
    counts = counts.sort_values("sample_id")
    counts = counts.set_index("taxon_id").join(taco.taxa.set_index("taxon_id"))
    
    plot_options = []
    if how == "bar":
        plot_options.append(so.Bar())
    else:
        plot_options.append(so.Area(alpha=.7))
    plot_options.append(so.Stack())
    figure = so.Plot(
        counts, 
        x="sample_id", 
        y="rel_abundance", 
        color="taxon_name"
    ).add(*plot_options)
    return figure