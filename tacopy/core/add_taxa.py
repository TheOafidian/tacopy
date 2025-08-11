import pandas as pd


def calculate_prevalence(taco):
    """Add the prevalence of a taxa to the taxa table in a Pyditaco object"""

    prevalence = taco.counts["taxon_id"].value_counts() / taco.samples.shape[0]
    prevalence.name = "prevalence"
    return prevalence.reset_index()

def _rename_taxa_by_most_x(taco, x, n, rank):
    most_x_taxa = taco.taxa.sort_values(x, ascending=False).head(
        n - 1
    )[["taxon_id", rank]]

    new_names = pd.merge(
        taco.taxa,
        most_x_taxa,
        on="taxon_id",
        how="left",
        suffixes=("_prev", None),
    )[["taxon_id", rank]].\
        set_index("taxon_id").\
            fillna("Other").\
                rename(columns={rank: "taxon_name"})

    taco.taxa = pd.merge(taco.taxa, new_names, left_on="taxon_id", right_on="taxon_id", how="left")
    return taco

def name_n_most_abundant(taco, n):
    """Sorts taxa by mean relative abundance and assigns a name to the n most abundant taxa under the taxon_name column"""
    keep_mean_rel_abund = True
    if not "mean_rel_abundance" in taco.taxa.columns:
        taco = taco.add_mean_rel_abundance()
        keep_mean_rel_abund = False
    
    lca_rank = taco.rank_names[-1]
    
    taco = _rename_taxa_by_most_x(taco, "mean_rel_abundance", n, lca_rank)
    if not keep_mean_rel_abund:
        taco.taxa = taco.taxa.drop(columns=["mean_rel_abundance"])

    return taco

def name_n_most_prevalent(taco, n):
    """Sorts taxa by prevalence and assigns a name to the n most prevalent taxa under the taxon_name column"""
    keep_prevalence = True
    if not "prevalence" in taco.taxa.columns:
        taco = taco.add_prevalence()
        keep_prevalence = False

    lca_rank = taco.rank_names[-1]
    taco = _rename_taxa_by_most_x(taco, "prevalence", n, lca_rank)

    if not keep_prevalence:
        taco.taxa = taco.taxa.drop(columns=["prevalence"])
    
    return taco
    