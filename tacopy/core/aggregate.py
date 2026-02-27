import copy


def aggregate_taxa(taco, rank="genus"):
    """Aggregate taxa at a given rank in a Pyditaco object"""

    if not rank in taco.taxa.columns:
        raise ValueError(f"The rank '{rank}' is missing from the taxa table")

    taco_agg = copy.copy(taco)
    if rank == "taxon_name":
        taxa = taco_agg.taxa[["taxon_id","taxon_name"]]
    else:
        tax_columns = ["taxon_id"]
        tax_columns.extend(taco_agg.rank_names[: taco_agg.rank_names.index(rank) + 1])
        taxa = taco_agg.taxa[tax_columns]

    # Fill up NAN values from left to right
    # and add "unknown" if highest level is not classified
    top_level = taco_agg._extract_rank_names()[0]
    taxa.loc[taxa[top_level].isna(), top_level] = "unknown"
    taxa = taxa.ffill(axis=1)

    # Aggregate counts
    counts = taco_agg.counts.merge(taxa[["taxon_id", rank]], on="taxon_id")
    counts = counts.groupby(["sample_id", rank])["count"].sum().reset_index()

    # Aggregate the taxa table
    taco_agg.taxa = taxa.drop(columns=["taxon_id"]).drop_duplicates()

    # Assign new taxon_ids to counts
    taco_agg.counts = counts.set_index(rank)\
        .join(taco_agg.taxa.set_index(rank))\
        .reset_index()\
            [
        ["sample_id", "count", "taxon_id"]
    ]
    
    if "taxon_name" in taco_agg.taxa.columns and not rank == "taxon_name":
            taco_agg.taxa = taco_agg.taxa.drop(columns=["taxon_name"])
    # Update rank names
    taco_agg.rank_names = taco_agg._extract_rank_names()

    return taco_agg
