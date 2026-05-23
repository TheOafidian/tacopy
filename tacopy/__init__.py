from tacopy.core.tacopy import Tacopy
from tacopy.core.io import read_tacos, write_tacos
from tacopy.core.aggregate import aggregate_taxa, aggregate_samples
from tacopy.viz.stack import tacoplot_stack

_dependencies_ = ["pandas", "skbio", "sklearn"]
_missing_dependencies_ = []

for _dep in _dependencies_:
    try:
        __import__(_dep)
    except ImportError:
        _missing_dependencies_.append(_dep)

if _missing_dependencies_:
    raise ImportError(f"Missing dependencies: {_missing_dependencies_}")
del _dep, _dependencies_, _missing_dependencies_

__all__ = ["Tacopy", "read_pyditacos", "write_pyditacos"]
