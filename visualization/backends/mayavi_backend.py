"""Mayavi visualization backend (optional).

Mayavi is an optional backend that requires:
- conda install mayavi -c conda-forge
- Enable via: VIZ_BACKEND=mayavi environment variable

Use when:
- You already have Mayavi installed
- You need specific Mayavi features
- You want visual consistency with existing Mayavi-based workflows
"""

from typing import List, Tuple, Union

import numpy as np

from .base import VisualizationBackend


class MayaviBackend(VisualizationBackend):
    """Mayavi-based visualization backend."""

    def __init__(self):
        # Import mayavi lazily to avoid import errors when not installed
        from mayavi import mlab

        self.mlab = mlab
        self._figure = None

    @property
    def name(self) -> str:
        return "mayavi"

    def create_figure(
        self,
        title: str = "Figure",
        bgcolor: Tuple[float, float, float] = (1.0, 1.0, 1.0),
        size: Tuple[int, int] = (600, 600),
    ) -> None:
        """Create a new Mayavi figure."""
        self._figure = self.mlab.figure(title, bgcolor=bgcolor, size=size)
        self._size = size
        self.mlab.clf()

    def clear(self) -> None:
        """Clear the current figure."""
        self.mlab.clf()

    def add_contour3d(
        self,
        X: np.ndarray,
        Y: np.ndarray,
        Z: np.ndarray,
        scalars: np.ndarray,
        contours: Union[int, List[float]],
        opacity: float = 0.5,
    ) -> None:
        """Add isosurface contours using Mayavi's contour3d."""
        self.mlab.contour3d(X, Y, Z, scalars, contours=contours, opacity=opacity)

    def add_sphere(
        self,
        x: float,
        y: float,
        z: float,
        radius: float,
        color: Tuple[float, float, float],
        resolution: int = 20,
    ) -> None:
        """Add a sphere at the given position."""
        # Mayavi's points3d uses scale_factor as diameter, not radius
        self.mlab.points3d(
            x,
            y,
            z,
            scale_factor=radius * 2,
            resolution=resolution,
            color=color,
            scale_mode="none",
        )

    def add_tube(
        self,
        points: np.ndarray,
        radius: float,
        color: Tuple[float, float, float],
    ) -> None:
        """Add a tube connecting the given points."""
        if len(points) < 2:
            return

        x = points[:, 0]
        y = points[:, 1]
        z = points[:, 2]
        self.mlab.plot3d(x, y, z, tube_radius=radius, color=color)

    def add_text3d(
        self,
        x: float,
        y: float,
        z: float,
        text: str,
        color: Tuple[float, float, float],
        scale: float,
    ) -> None:
        """Add a 3D text label at the given position."""
        self.mlab.text3d(x, y, z, text, color=color, scale=(scale, scale, scale))

    def add_outline(self) -> None:
        """Add a bounding box outline to the scene."""
        self.mlab.outline()

    def add_scalar_bar(self, title: str = "") -> None:
        """Add a scalar bar (color legend) to the scene."""
        self.mlab.scalarbar(title=title)

    def show(self) -> None:
        """Display the visualization window."""
        self.mlab.show()

    def save(self, filename: str) -> None:
        """Save the visualization to an image file."""
        # Use stored size to ensure correct output dimensions in offscreen mode
        if hasattr(self, "_size") and self._size:
            self.mlab.savefig(filename, size=self._size)
        else:
            self.mlab.savefig(filename)

    def close(self) -> None:
        """Close the visualization window."""
        self.mlab.close()
