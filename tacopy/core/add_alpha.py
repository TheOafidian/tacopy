import copy
import numpy as np

def add_obs(taco):

    taco_alpha = copy.copy(taco)
    observed = taco_alpha.counts.groupby("sample_id")["taxon_id"].nunique()
    observed.name = "obs"
    taco_alpha.samples = taco_alpha.samples.merge(observed, on="sample_id")
    return taco_alpha

def add_shannon(taco):

    taco_alpha = copy.copy(taco)
    rel = taco.add_rel_abundance().counts.copy()
    rel["pilnpi"] = - rel["rel_abundance"] * np.log(rel["rel_abundance"])
    shannon = rel.groupby("sample_id")["pilnpi"].sum()

    shannon.name = "shannon"
    taco_alpha.samples = taco_alpha.samples.merge(shannon, on="sample_id", how="left")

    return taco_alpha
 
    