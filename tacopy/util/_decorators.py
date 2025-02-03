from collections.abc import Callable
from tacopy.util._typing import F


def set_module(module) -> Callable[[F], F]:
    """Private decorator for overriding __module__ on a function or class.

    Example usage::

        @set_module("pandas")
        def example():
            pass


        assert example.__module__ == "pandas"
    """

    def decorator(func: F) -> F:
        if module is not None:
            func.__module__ = module
        return func

    return decorator
