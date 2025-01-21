def calculate_prevalence(taco):
    """Add the prevalence of a taxa to the taxa table in a Pyditaco object"""

    prevalence = taco.counts["taxon_id"].value_counts() / taco.samples.shape[0]
    prevalence.name = "prevalence"
    return prevalence.reset_index()


def name_n_most_abundant(taco, n, taxon="genus"):
    if not "prevalence" in taco.taxa.columns:
        taco = taco.add_prevalence(taxon=taxon)

    most_abundant_taxa = taco.taxa.sort_values("prevalence", ascending=False).head(n)
    
