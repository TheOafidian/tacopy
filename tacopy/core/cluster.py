from sklearn.cluster import AgglomerativeClustering, KMeans


def cluster_samples(taco, method="kmeans", *args, **kwargs):
    models = {"kmeans": KMeans, "agglomerative": AgglomerativeClustering}
    model = models[method](*args, **kwargs)
    counts = taco.counts.pivot_table(
        columns="taxon_id", values="count", index="sample_id"
    ).fillna(0)
    model = model.fit(counts)

    return model
