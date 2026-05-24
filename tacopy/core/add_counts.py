import pandas as pd


def add_relative_abundances(
    count_table: pd.DataFrame, sample:str="sample_id", overwrite:bool=False
):
    """
        Add per-sample relative abundances to a counts DataFrame.

    Relative abundance is computed as each count divided by the total count
    for its sample, and stored in a new ``'rel_abundance'`` column.

    Parameters
    ----------
    count_table : pd.DataFrame
        A long-format counts table with at least a ``'count'`` column and a
        column matching ``sample``.
    sample : str, optional
        The column name identifying samples to group by when computing
        totals. Defaults to ``'sample_id'``.
    overwrite : bool, optional
        If ``False`` (default), operate on a copy of the DataFrame so the
        original is left unchanged. If ``True``, add the ``'rel_abundance'``
        column in-place.

    Returns
    -------
    pd.DataFrame
        A DataFrame with a ``'rel_abundance'`` column added. When
        ``overwrite=False`` this is a new copy; when ``overwrite=True`` it
        is the same object as ``count_table``.
    """
    if overwrite:
        count_table["rel_abundance"] = count_table["count"] / count_table.groupby(sample)[
            "count"
        ].transform("sum")
    else:
        count_table = count_table.copy()
        count_table["rel_abundance"] = count_table["count"] / count_table.groupby(
            sample
        )["count"].transform("sum")
    return count_table
