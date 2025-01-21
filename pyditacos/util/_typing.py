from typing import TypeVar, Any
from collections.abc import Callable

FuncType = Callable[..., Any]
F = TypeVar("F", bound=FuncType)