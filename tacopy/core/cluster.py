import copy
import numpy as np
import pandas as pd
from sklearn.cluster import AgglomerativeClustering, KMeans
from skbio.diversity import (
    beta_diversity,
 )

def add_hier_clust_order(taco, taxon_name='taxon_id', sample_name='sample_id'):
    taco = copy.copy(taco)
    matrix = taco.matrix(values="rel_abundance").T
    matrix_beta = beta_diversity("braycurtis",matrix)
    model = AgglomerativeClustering(distance_threshold=0, 
                                    n_clusters=None, linkage="average",
                                    metric="precomputed")
    model.fit(matrix_beta.data)
    taco.samples["sample_clustered"] = pd.Categorical(
        taco.samples.sample_id, categories=np.array(matrix_beta.ids)[model.labels_], ordered=True
    )
    return taco

def cluster_samples(taco, method="hierarchical", *args, **kwargs):
    models = {"kmeans": KMeans, "agglomerative": AgglomerativeClustering}
    model = models[method](*args, **kwargs)
    counts = taco.counts.pivot_table(
        columns="taxon_id", values="count", index="sample_id"
    ).fillna(0)
    model = model.fit(counts)

    return model
