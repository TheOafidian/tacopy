import pandas as pd


def add_relative_abundances(count_table: pd.DataFrame, sample="sample_id"):
    """Return the relative abundance to the count table"""
    count_table["rel_abundance"] = rel_ct = count_table["count"] / count_table.groupby(
        sample
    )["count"].transform("sum")
    return count_table
