"""PyVista visualization backend (default).

PyVista is the default backend because:
- Simple installation: pip install pyvista
- Excellent documentation
- Active community
- Pythonic API
"""

from typing import List, Tuple, Union

import numpy as np
import pyvista as pv

from .base import VisualizationBackend


class PyVistaBackend(VisualizationBackend):
    """PyVista-based visualization backend."""

    def __init__(self):
        self._plotter: pv.Plotter = None

    @property
    def name(self) -> str:
        return "pyvista"

    def create_figure(
        self,
        title: str = "Figure",
        bgcolor: Tuple[float, float, float] = (1.0, 1.0, 1.0),
        size: Tuple[int, int] = (600, 600),
    ) -> None:
        """Create a new PyVista plotter."""
        self._plotter = pv.Plotter(window_size=size)
        self._plotter.set_background(bgcolor)
        # Enable depth peeling for correct transparency rendering
        self._plotter.enable_depth_peeling()

    def clear(self) -> None:
        """Clear the current plotter."""
        if self._plotter is not None:
            self._plotter.clear()

    def add_contour3d(
        self,
        X: np.ndarray,
        Y: np.ndarray,
        Z: np.ndarray,
        scalars: np.ndarray,
        contours: Union[int, List[float]],
        opacity: float = 0.5,
    ) -> None:
        """Add isosurface contours using PyVista's contour filter."""
        if self._plotter is None:
            raise RuntimeError("No figure created. Call create_figure() first.")

        # Create structured grid from coordinate arrays
        grid = pv.StructuredGrid(X, Y, Z)
        # Flatten scalars in Fortran order for correct VTK interpretation
        grid["scalars"] = scalars.flatten(order="F")

        # Handle contour specification
        if isinstance(contours, int):
            # Generate evenly spaced contours
            vmin, vmax = scalars.min(), scalars.max()
            contour_values = np.linspace(vmin, vmax, contours + 2)[1:-1].tolist()
        else:
            contour_values = contours

        # Extract isosurfaces
        try:
            mesh = grid.contour(isosurfaces=contour_values, scalars="scalars")
            if mesh.n_points > 0:
                # Use coolwarm colormap for orbital +/- lobes
                self._plotter.add_mesh(mesh, opacity=opacity, cmap="coolwarm")
        except Exception:
            # Contour may fail for degenerate data
            pass

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
        if self._plotter is None:
            raise RuntimeError("No figure created. Call create_figure() first.")

        sphere = pv.Sphere(
            radius=radius,
            center=(x, y, z),
            theta_resolution=resolution,
            phi_resolution=resolution,
        )
        self._plotter.add_mesh(sphere, color=color, smooth_shading=True)

    def add_tube(
        self,
        points: np.ndarray,
        radius: float,
        color: Tuple[float, float, float],
        n_sides: int = 20,
    ) -> None:
        """Add a tube connecting the given points."""
        if self._plotter is None:
            raise RuntimeError("No figure created. Call create_figure() first.")

        if len(points) < 2:
            return

        # Create line and tube with higher resolution for smooth appearance
        line = pv.Line(points[0], points[-1])
        tube = line.tube(radius=radius, n_sides=n_sides)
        self._plotter.add_mesh(tube, color=color, smooth_shading=True)

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
        if self._plotter is None:
            raise RuntimeError("No figure created. Call create_figure() first.")

        # PyVista uses point labels (billboard-style, always face camera)
        self._plotter.add_point_labels(
            [(x, y, z)],
            [text],
            font_size=int(scale * 20),
            text_color=color,
            shape=None,
            always_visible=True,
        )

    def add_outline(self) -> None:
        """Add a bounding box to the scene."""
        if self._plotter is None:
            raise RuntimeError("No figure created. Call create_figure() first.")

        self._plotter.add_bounding_box()

    def add_scalar_bar(self, title: str = "") -> None:
        """Add a scalar bar (color legend) to the scene."""
        if self._plotter is None:
            raise RuntimeError("No figure created. Call create_figure() first.")

        self._plotter.add_scalar_bar(title=title)

    def show(self) -> None:
        """Display the visualization window."""
        if self._plotter is None:
            raise RuntimeError("No figure created. Call create_figure() first.")

        self._plotter.show()

    def save(self, filename: str) -> None:
        """Save the visualization to an image file."""
        if self._plotter is None:
            raise RuntimeError("No figure created. Call create_figure() first.")

        # For PyVista, we need to render before screenshot if not shown
        # show(screenshot=filename) handles this properly
        self._plotter.show(screenshot=filename, auto_close=False)

    def close(self) -> None:
        """Close the visualization window."""
        if self._plotter is not None:
            self._plotter.close()
            self._plotter = None
