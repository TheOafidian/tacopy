import copy
from pathlib import Path
import pandas as pd
from tacopy.core.tacopy import Tacopy

def read_tacos(tacodir):
    """
    Read a tidytacos or tidyamplicons directory into a Tacopy object.

    Expects the directory to contain ``samples.csv``, ``taxa.csv``, and
    either ``counts.csv`` or the legacy ``abundances.csv``.

    Parameters
    ----------
    tacodir : str or Path
        Path to the tidytacos or tidyamplicons directory to read from.

    Returns
    -------
    Tacopy
        A Tacopy object populated with the samples, taxa, and counts tables.

    Raises
    ------
    FileNotFoundError
        If neither ``counts.csv`` nor ``abundances.csv`` is found in the
        directory.
    """
    samples = pd.read_csv(Path(tacodir) / "samples.csv")
    taxa = pd.read_csv(Path(tacodir) / "taxa.csv")
    try:
        counts = pd.read_csv(Path(tacodir) / "counts.csv")
    except FileNotFoundError:
        counts = pd.read_csv(Path(tacodir) / "abundances.csv")
        counts.rename({"count": "abundance"})
    return Tacopy(samples, taxa, counts)


def write_tacos(tacopy, tacodir):
    """
    Write a Tacopy object to a tidytacos directory.

    Creates the target directory if it does not already exist, then writes
    the samples, taxa, and counts tables as CSV files.

    Parameters
    ----------
    tacopy : Tacopy
        The Tacopy object to serialise.
    tacodir : str or Path
        Path to the output directory. Created automatically if absent.

    Returns
    -------
    str or Path
        The path to the directory the files were written to (echoes
        ``tacodir``).
    """
    Path(tacodir).mkdir(exist_ok=True)
    tacopy.samples.to_csv(Path(tacodir) / "samples.csv", index=False)
    tacopy.taxa.to_csv(Path(tacodir) / "taxa.csv", index=False)
    tacopy.counts.to_csv(Path(tacodir) / "counts.csv", index=False)
    return tacodir


def write_formula(fmls:list, outf:str):

    with open(Path(outf), "w") as f:
        f.write("\n".join(fmls))


def write_mda_inputs(ta:Tacopy, count_matrix_outf:str ="counts.tsv", metadata_outf:str = "meta.tsv"):

    tt = copy.copy(ta)
    if not "taxon_name" in tt.taxa.columns:
        tt = tt.add_taxon_name(rank=tt.rank_names[-1], n=None)

    M = tt.matrix(taxon="taxon_name").T
    M.index.name = "SampleID"
    M.to_csv(Path(count_matrix_outf), sep="\t")
    ta.samples.set_index("sample_id").loc[M.index]\
        .to_csv(Path(metadata_outf), sep="\t")
    return count_matrix_outf, metadata_outf
