from pathlib import Path
import pandas as pd
from tacopy.core.tacopy import Tacopy


def read_tacos(tacodir):
    """Read a tidytacos folder as a tacopy object"""
    samples = pd.read_csv(Path(tacodir) / "samples.csv")
    taxa = pd.read_csv(Path(tacodir) / "taxa.csv")
    try:
        counts = pd.read_csv(Path(tacodir) / "counts.csv")
    except FileNotFoundError:
        counts = pd.read_csv(Path(tacodir) / "abundances.csv")
        counts.rename({"count": "abundance"})
    return Tacopy(samples, taxa, counts)


def write_tacos(tacopy, tacodir):
    """Write a tacopy object to a tidytacos folder"""
    Path(tacodir).mkdir(exist_ok=True)
    tacopy.samples.to_csv(Path(tacodir) / "samples.csv", index=False)
    tacopy.taxa.to_csv(Path(tacodir) / "taxa.csv", index=False)
    tacopy.counts.to_csv(Path(tacodir) / "counts.csv", index=False)
    return tacodir
