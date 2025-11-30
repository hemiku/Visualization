"""Visualization backend selection and factory.

Default backend: PyVista (simple pip install)
Optional backend: Mayavi (requires conda, enable via VIZ_BACKEND=mayavi)

Usage:
    from visualization.backends import get_backend

    backend = get_backend()  # Uses default (pyvista) or VIZ_BACKEND env var
    backend = get_backend('mayavi')  # Explicitly use mayavi
"""

import os
from typing import Optional

from .base import VisualizationBackend

__all__ = ["get_backend", "set_backend", "VisualizationBackend"]

_backend_instance: Optional[VisualizationBackend] = None


def get_backend(name: Optional[str] = None) -> VisualizationBackend:
    """Get a visualization backend instance.

    Args:
        name: Backend name ('pyvista' or 'mayavi').
              If None, uses VIZ_BACKEND env var, defaulting to 'pyvista'.

    Returns:
        VisualizationBackend instance.

    Raises:
        ImportError: If the requested backend is not installed.
        ValueError: If an unknown backend name is provided.
    """
    global _backend_instance

    if name is None:
        name = os.environ.get("VIZ_BACKEND", "pyvista").lower()

    # Return cached instance if it matches requested backend
    if _backend_instance is not None:
        if getattr(_backend_instance, "name", None) == name:
            return _backend_instance

    # Create new backend instance
    if name == "mayavi":
        from .mayavi_backend import MayaviBackend

        _backend_instance = MayaviBackend()
    elif name == "pyvista":
        from .pyvista_backend import PyVistaBackend

        _backend_instance = PyVistaBackend()
    else:
        raise ValueError(
            f"Unknown backend: '{name}'. Available: 'pyvista', 'mayavi'"
        )

    return _backend_instance


def set_backend(name: str) -> VisualizationBackend:
    """Set and return a new backend, clearing any cached instance.

    Args:
        name: Backend name ('pyvista' or 'mayavi').

    Returns:
        New VisualizationBackend instance.
    """
    global _backend_instance
    _backend_instance = None
    return get_backend(name)
