"""Docking methods - simplified function-based approach."""

# Lazy imports to avoid loading methods that are not needed
# Methods are imported only when accessed

# Dictionary of methods with their functions
# All methods are imported lazily to avoid errors if modules are incomplete
_methods_cache = None


def _get_methods():
    """Get METHODS dictionary with lazy loading."""
    global _methods_cache
    if _methods_cache is not None:
        return _methods_cache

    methods = {}

    # Gnina - lazy import
    try:
        from . import gnina

        methods["gnina"] = {
            "dock": gnina.dock_gnina,
            "extract_metrics": gnina.extract_metrics_gnina,
        }
    except (ImportError, AttributeError):
        pass  # Gnina not available

    # QVina - lazy import
    try:
        from . import qvina

        methods["qvina"] = {
            "dock": qvina.dock_qvina,
            "extract_metrics": qvina.extract_metrics_qvina,
        }
    except (ImportError, AttributeError):
        pass  # QVina not available

    # PLAPT - lazy import, only if available and complete
    try:
        from . import plapt

        if hasattr(plapt, "dock_plapt") and hasattr(plapt, "extract_metrics_plapt"):
            methods["plapt"] = {
                "dock": plapt.dock_plapt,
                "extract_metrics": plapt.extract_metrics_plapt,
            }
    except (ImportError, AttributeError):
        pass  # PLAPT not available or incomplete

    # DynamicBind - lazy import, only if available and complete
    try:
        from . import dynamicbind

        if (
            hasattr(dynamicbind, "dock_dynamicbind")
            and hasattr(dynamicbind, "extract_metrics_dynamicbind")
            and hasattr(dynamicbind, "preprocess_dynamicbind")
        ):
            methods["dynamicbind"] = {
                "dock": dynamicbind.dock_dynamicbind,
                "extract_metrics": dynamicbind.extract_metrics_dynamicbind,
                "preprocess": dynamicbind.preprocess_dynamicbind,
            }
    except (ImportError, AttributeError):
        pass  # DynamicBind not available or incomplete

    # Interformer (bench variant) - lazy import, only if available and complete
    try:
        from . import interformer

        if (
            hasattr(interformer, "dock_interformer")
            and hasattr(interformer, "extract_metrics_interformer")
            and hasattr(interformer, "preprocess_interformer")
        ):
            methods["interformer"] = {
                "dock": interformer.dock_interformer,
                "extract_metrics": interformer.extract_metrics_interformer,
                "preprocess": interformer.preprocess_interformer,
            }
    except (ImportError, AttributeError):
        pass  # Interformer not available or incomplete

    # NOTE: Method integration removed by request.

    _methods_cache = methods
    return methods


# Create METHODS as a property-like accessor
class _MethodsDict:
    """Lazy-loaded methods dictionary."""

    def __getitem__(self, key):
        return _get_methods()[key]

    def __contains__(self, key):
        return key in _get_methods()

    def keys(self):
        return _get_methods().keys()

    def values(self):
        return _get_methods().values()

    def items(self):
        return _get_methods().items()

    def get(self, key, default=None):
        return _get_methods().get(key, default)


METHODS = _MethodsDict()

__all__ = ["METHODS"]

