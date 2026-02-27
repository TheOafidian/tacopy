from sklearn.cluster import AgglomerativeClustering

def cluster_samples(taco, *args, **kwargs):
    
    model = AgglomerativeClustering(*args, **kwargs)
    counts = taco.counts.\
        pivot_table(columns="taxon_id", values="count", index="sample_id").\
            fillna(0)
    model = model.fit(counts)
    
    return model
    