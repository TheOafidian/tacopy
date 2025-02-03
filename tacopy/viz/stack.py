from pyditacos import Pyditaco
import pyditacos as pt
import matplotlib.pyplot as plt
import seaborn as sns


def tacoplot_stack(taco: Pyditaco, x="sample", *args, **kwargs):
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
    ax.stackplot(
        pt.samples[x],
    )
