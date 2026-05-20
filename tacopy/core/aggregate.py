import copy
from tacopy import Tacopy

def aggregate_taxa(taco, rank="genus"):
    """Aggregate taxa at a given rank in a Pyditaco object"""

    if not rank in taco.taxa.columns:
        raise ValueError(f"The rank '{rank}' is missing from the taxa table")

    taco_agg = copy.copy(taco)
    if rank == "taxon_name":
        taxa = taco_agg.taxa[["taxon_id", "taxon_name"]]
        ranks = ["taxon_id", "taxon_name"]
        samps = ["sample_id","taxon_name"]
    else:
        tax_columns = ["taxon_id"]
        tax_columns.extend(taco_agg.rank_names[: taco_agg.rank_names.index(rank) + 1])
        taxa = taco_agg.taxa[tax_columns]
        ranks = ["taxon_id"]
        samps = ["sample_id"]
        ranks.extend(taco_agg.rank_names[: taco_agg.rank_names.index(rank) + 1])
        samps.extend(taco_agg.rank_names[: taco_agg.rank_names.index(rank) + 1])

    # Fill up NAN values from left to right
    # and add "unknown" if highest level is not classified
    top_level = taco_agg._extract_rank_names()[0]
    taxa.loc[taxa[top_level].isna(), top_level] = "unknown"
    taxa = taxa.ffill(axis=1)

    taxids_new = taxa\
      .groupby(ranks[1:])["taxon_id"]\
      .apply(list)\
      .reset_index(drop=True)\
      .reset_index()\
      .rename(
        columns = {
          "taxon_id": "taxon_id_old", 
          "index" : "taxon_id"})
    
    taxids_new["taxon_id"] = "t" + taxids_new["taxon_id"].astype(str)
    taxids_new = taxids_new.explode("taxon_id_old")
    taxa = taxa.rename(columns={"taxon_id":"taxon_id_old"})\
      .merge(taxids_new, on="taxon_id_old", validate="one_to_many")

    # Aggregate counts
    counts = taco_agg.counts.rename(columns={"taxon_id":"taxon_id_old"})

    counts = counts.merge(taxa, on="taxon_id_old", how="left", validate="many_to_one")
    counts = counts[["sample_id","taxon_id","count"]].groupby(["sample_id","taxon_id"])["count"].sum().reset_index()
 
    taco_agg = Tacopy(
        taco_agg.samples, 
        # Aggregate the taxa table
        taxa.drop(columns=["taxon_id_old"]).drop_duplicates(),
        counts)

    if "taxon_name" in taco_agg.taxa.columns and not rank == "taxon_name":
        taco_agg.taxa = taco_agg.taxa.drop(columns=["taxon_name"])
    # Update rank names
    taco_agg.rank_names = taco_agg._extract_rank_names()

    taco_agg.taxa["taxon_name"] = taco_agg.taxa[rank]

    return taco_agg



def aggregate_samples(taco, group="sample"):
    """Aggregate samples to a given group in a Pyditaco object"""

    if not group in taco.samples.columns:
        raise ValueError(f"The column '{group}' is missing from the sample table")

    taco_agg = copy.copy(taco)
    samples = taco_agg.samples[["sample_id", group]]

    # Aggregate counts
    counts = taco_agg.counts.merge(samples, on="sample_id")
    counts = counts.groupby(["taxon_id", group])["count"].sum().reset_index()

    taco_agg.samples = samples[[group]].drop_duplicates()
    taco_agg.counts = counts.merge(taco_agg.samples, on=group).drop(columns=[group])

    taco_agg.samples = taco_agg.samples.rename(columns={group: "sample"})
    return taco_agg
