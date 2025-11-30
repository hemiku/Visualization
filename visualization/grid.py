"""3D grid generation for molecular orbital calculations."""

from typing import List, Optional

import numpy as np

from visualization.constants import (
    DEFAULT_R_MAX_MULTIPLIER,
    DEFAULT_GRID_RESOLUTION,
    BOHR_EXTENT_FACTOR,
)


class Grid:
    """3D spatial grid for orbital calculations.

    The grid extends beyond the molecular boundaries based on the
    spatial extent of the basis functions.
    """

    def __init__(
        self,
        R_max_multip: float = DEFAULT_R_MAX_MULTIPLIER,
        x_min: Optional[float] = None,
        x_max: Optional[float] = None,
        x_n: int = DEFAULT_GRID_RESOLUTION,
        y_min: Optional[float] = None,
        y_max: Optional[float] = None,
        y_n: int = DEFAULT_GRID_RESOLUTION,
        z_min: Optional[float] = None,
        z_max: Optional[float] = None,
        z_n: int = DEFAULT_GRID_RESOLUTION,
    ):
        self.R_max_multip = R_max_multip
        self.x_min = x_min
        self.x_max = x_max
        self.x_n = x_n
        self.y_min = y_min
        self.y_max = y_max
        self.y_n = y_n
        self.z_min = z_min
        self.z_max = z_max
        self.z_n = z_n

    def return_grid_arrays(self) -> np.ndarray:
        """Return 3D meshgrid arrays for X, Y, Z coordinates."""
        return np.mgrid[
            self.x_min : self.x_max : self.x_n * 1j,
            self.y_min : self.y_max : self.y_n * 1j,
            self.z_min : self.z_max : self.z_n * 1j,
        ]

    def generate_grid_boundaries(self, basis: List, atoms_R: np.ndarray) -> None:
        """Calculate grid boundaries from basis set extent and atomic positions.

        Args:
            basis: List of basis set data per atom
            atoms_R: Atomic positions array (n_atoms x 3)
        """
        # Calculate maximum spatial extent of each basis function
        basis_extents = []
        for atom_basis in basis:
            for orbital_basis in atom_basis:
                max_exponent = orbital_basis.max()
                extent = BOHR_EXTENT_FACTOR / np.sqrt(max_exponent)
                basis_extents.append(extent)

        R_max = self.R_max_multip * np.max(basis_extents)

        self.x_max = np.max(atoms_R[:, 0]) + R_max
        self.x_min = np.min(atoms_R[:, 0]) - R_max

        self.y_max = np.max(atoms_R[:, 1]) + R_max
        self.y_min = np.min(atoms_R[:, 1]) - R_max

        self.z_max = np.max(atoms_R[:, 2]) + R_max
        self.z_min = np.min(atoms_R[:, 2]) - R_max
