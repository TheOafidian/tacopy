from pyditacos.core.pyditaco import Pyditaco
from pyditacos.core.io import read_pyditacos, write_pyditacos
from pyditacos.core.aggregate import aggregate_taxa


_dependencies_ = ["pandas"]
_missing_dependencies_ = []

for _dep in _dependencies_:
    try:
        __import__(_dep)
    except ImportError:
        _missing_dependencies_.append(_dep)

if _missing_dependencies_:
    raise ImportError(f"Missing dependencies: {_missing_dependencies_}")
del _dep, _dependencies_, _missing_dependencies_

__all__ = ["Pyditaco", "read_pyditacos", "write_pyditacos"]
