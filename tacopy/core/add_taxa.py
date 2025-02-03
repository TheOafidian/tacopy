import pandas as pd


def calculate_prevalence(taco):
    """Add the prevalence of a taxa to the taxa table in a Pyditaco object"""

    prevalence = taco.counts["taxon_id"].value_counts() / taco.samples.shape[0]
    prevalence.name = "prevalence"
    return prevalence.reset_index()

def name_n_most_abundant(taco, n):
    pass

def name_n_most_prevalent(taco, n):
    """Sorts taxa by prevalence and assigns a name to the n most prevalent taxa under the taxon_name column"""
    keep_prevalence = True
    if not "prevalence" in taco.taxa.columns:
        taco = taco.add_prevalence()
        keep_prevalence = False

    lca_rank = taco.rank_names[-1]

    most_abundant_taxa = taco.taxa.sort_values("prevalence", ascending=False).head(
        n - 1
    )[["taxon_id", lca_rank]]

    new_names = pd.merge(
        taco.taxa,
        most_abundant_taxa,
        on="taxon_id",
        how="left",
        suffixes=("_prev", None),
    )[["taxon_id", lca_rank]].\
        set_index("taxon_id").\
            fillna("Other").\
                rename(columns={lca_rank: "taxon_name"})

    taco.taxa = taco.taxa.join(new_names, on="taxon_id")
    
    if not keep_prevalence:
        taco.taxa = taco.taxa.drop(columns=["prevalence"])
    
    return taco
    
def name_n_most_abundant(taco, n):
    pass