"""Abstract base class for visualization backends."""

from abc import ABC, abstractmethod
from typing import Tuple, List, Union
import numpy as np


class VisualizationBackend(ABC):
    """Abstract base for visualization backends (PyVista, Mayavi, etc.)."""

    @property
    @abstractmethod
    def name(self) -> str:
        """Return backend name."""
        ...

    @abstractmethod
    def create_figure(
        self,
        title: str = "Figure",
        bgcolor: Tuple[float, float, float] = (1.0, 1.0, 1.0),
        size: Tuple[int, int] = (600, 600),
    ) -> None:
        """Create a new figure/plotter window."""
        ...

    @abstractmethod
    def clear(self) -> None:
        """Clear the current figure."""
        ...

    @abstractmethod
    def add_contour3d(
        self,
        X: np.ndarray,
        Y: np.ndarray,
        Z: np.ndarray,
        scalars: np.ndarray,
        contours: Union[int, List[float]],
        opacity: float = 0.5,
    ) -> None:
        """Add isosurface contours to the scene."""
        ...

    @abstractmethod
    def add_sphere(
        self,
        x: float,
        y: float,
        z: float,
        radius: float,
        color: Tuple[float, float, float],
        resolution: int = 20,
    ) -> None:
        """Add a sphere (for atoms) at the given position."""
        ...

    @abstractmethod
    def add_tube(
        self,
        points: np.ndarray,
        radius: float,
        color: Tuple[float, float, float],
    ) -> None:
        """Add a tube/cylinder (for bonds) connecting points."""
        ...

    @abstractmethod
    def add_text3d(
        self,
        x: float,
        y: float,
        z: float,
        text: str,
        color: Tuple[float, float, float],
        scale: float,
    ) -> None:
        """Add 3D text label at the given position."""
        ...

    @abstractmethod
    def add_outline(self) -> None:
        """Add bounding box outline to the scene."""
        ...

    @abstractmethod
    def add_scalar_bar(self, title: str = "") -> None:
        """Add color bar / scalar bar to the scene."""
        ...

    @abstractmethod
    def show(self) -> None:
        """Display the visualization window."""
        ...

    @abstractmethod
    def save(self, filename: str) -> None:
        """Save the visualization to an image file."""
        ...

    @abstractmethod
    def close(self) -> None:
        """Close the visualization window."""
        ...
