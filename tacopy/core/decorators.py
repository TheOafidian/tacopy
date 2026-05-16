import functools
import inspect

def ensure_taco_columns(*checks):
    """
    Decorator that ensures a taco object has required columns/attributes
    before running the wrapped function.

    Each check is a tuple of (column_name, method_to_call_if_missing), e.g.:
        ("taxa","prevalence", "add_prevalence")
        ("counts", "rel_abundance", "add_rel_abundance")

    The decorator auto-detects which parameter is the taco object by
    looking for a parameter named 'taco', or falls back to the first argument.

    Usage:
        @ensure_taco_columns(
            ("taxa","prevalence", "add_prevalence"),
            ("counts", "rel_abundance", "add_rel_abundance"),
        )
        def my_analysis(taco, ...):
            ...
    """
    def decorator(func):
        sig = inspect.signature(func)
        param_names = list(sig.parameters.keys())

        taco_idx = next(
            (i for i, name in enumerate(param_names) if name == "taco"), 0
        )
        taco_param_name = param_names[taco_idx]

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            if taco_param_name in kwargs:
                taco = kwargs[taco_param_name]
                def set_taco(t): kwargs[taco_param_name] = t
            elif taco_idx  < len(args):
                taco = args[taco_idx]
                args = list(args)
                def set_taco(t): args[taco_idx] = t
            else:
                raise ValueError(
                    f"Could not find taco argument '{taco_param_name}' "
                    f"in call to {func.__name__}()"
                )
            
            for table, column, method in checks:
                if column not in taco[table].columns:
                    taco = getattr(taco, method)()
                    set_taco(taco)
            
            return func(*args, **kwargs)
        return wrapper
    return decorator