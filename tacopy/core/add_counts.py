import pandas as pd


def add_relative_abundances(
    count_table: pd.DataFrame, sample="sample_id", overwrite=False
):
    """Return the relative abundance to the count table"""
    if overwrite:
        count_table["count"] = count_table["count"] / count_table.groupby(sample)[
            "count"
        ].transform("sum")
    else:
        count_table["rel_abundance"] = count_table["count"] / count_table.groupby(
            sample
        )["count"].transform("sum")
    return count_table
