import copy


def aggregate_taxa(taco, rank="genus"):
    """Aggregate taxa at a given rank in a Pyditaco object"""

    if not rank in taco.taxa.columns:
        raise ValueError(f"The rank '{rank}' is missing from the taxa table")

    taco_agg = copy.copy(taco)
    tax_columns = ["taxon_id"]
    tax_columns.extend(taco_agg.rank_names[: taco_agg.rank_names.index(rank) + 1])
    taxa = taco_agg.taxa[tax_columns]

    # Fill up NAN values from left to right
    taxa = taxa.ffill(axis=1)

    # Aggregate counts
    counts = taco_agg.counts.merge(taxa[["taxon_id", rank]], on="taxon_id")
    counts = counts.groupby(["sample_id", rank])["count"].sum().reset_index()

    # Aggregate the taxa table
    taco_agg.taxa = taxa.drop(columns=["taxon_id"]).drop_duplicates()

    # Assign new taxon_ids to counts
    taco_agg.counts = counts.set_index(rank).join(taco_agg.taxa.set_index(rank))[
        ["sample_id", "count", "taxon_id"]
    ]
    # Update rank names
    taco_agg.rank_names = taco_agg._extract_rank_names()

    return taco_agg
