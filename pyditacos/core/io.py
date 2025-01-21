from pathlib import Path
import pandas as pd
from pyditacos.core.pyditaco import Pyditaco

def read_pyditacos(tacodir):
    """Read a tidytacos folder as a pyditacos object"""
    samples = pd.read_csv(Path(tacodir) / 'samples.csv')
    taxa = pd.read_csv(Path(tacodir) / 'taxa.csv')
    try:
        counts = pd.read_csv(Path(tacodir) / 'counts.csv')
    except FileNotFoundError:
        counts = pd.read_csv(Path(tacodir) / 'abundances.csv')
        counts.rename({"count":"abundance"})
    return Pyditaco(samples, taxa, counts)

def write_pyditacos(pyditaco, tacodir):
    """Write a pyditacos object to a tidytacos folder"""
    Path(tacodir).mkdir(exist_ok=True)
    pyditaco.samples.to_csv(Path(tacodir) / 'samples.csv', index=False)
    pyditaco.taxa.to_csv(Path(tacodir) / 'taxa.csv', index=False)
    pyditaco.counts.to_csv(Path(tacodir) / 'counts.csv', index=False)
    return tacodir
