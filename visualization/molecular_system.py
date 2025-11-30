"""Molecular system data container."""

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np


@dataclass
class MolecularSystem:
    """Container for molecular structure, basis set, and orbital data.

    This dataclass holds all data needed for molecular visualization:
    - Atomic structure (positions, charges, names, bonds)
    - Basis set information (basis functions, normalization)
    - Orbital coefficients and computed orbitals (AOs, MOs, geminals)

    Attributes can be accessed directly (e.g., `system.nAtoms`) or set
    directly (e.g., `system.nAtoms = 10`).
    """

    # Atomic structure
    nAtoms: int = 0
    atoms_R: Optional[np.ndarray] = None  # Atomic positions (nAtoms x 3)
    atoms_Charge: Optional[List] = None
    atoms_Name: Optional[List[str]] = None
    bonds: Optional[List] = None

    # Basis set
    spherical: bool = False
    nb: int = 0  # Number of basis functions
    basis: Optional[List] = None
    basis_norm: Optional[List] = None
    basis_norm2: Optional[List] = None

    # Orbital coefficients
    coeff: Optional[np.ndarray] = None
    Occ: Optional[np.ndarray] = None

    # Geminal data
    inactive: int = 0
    electrons: int = 0
    n_geminals: int = 0
    G_coeff: Optional[np.ndarray] = None
    Orb2Gem: Optional[np.ndarray] = None

    # Computed orbitals (filled by OrbitalsGenerator/GeminalGenerator)
    AOs: Optional[np.ndarray] = None
    MOs: Optional[np.ndarray] = None
    geminals: Optional[np.ndarray] = None

    # Backwards compatibility: keep Coeff as alias for coeff
    @property
    def Coeff(self) -> Optional[np.ndarray]:
        """Alias for coeff (backwards compatibility)."""
        return self.coeff

    @Coeff.setter
    def Coeff(self, value: np.ndarray) -> None:
        self.coeff = value

    def set_basis_and_norms(self, input: list) -> None:
        """Set basis, basis_norm, and basis_norm2 from a list.

        Args:
            input: List of [basis, basis_norm, basis_norm2]
        """
        self.basis = input[0]
        self.basis_norm = input[1]
        self.basis_norm2 = input[2]

    # Backwards compatibility methods (deprecated, use direct attribute access)
    def set_spherical(self, spherical: bool) -> None:
        self.spherical = spherical

    def set_nb(self, nb: int) -> None:
        self.nb = nb

    def set_nAtoms(self, nAtoms: int) -> None:
        self.nAtoms = nAtoms

    def set_inactive(self, inactive: int) -> None:
        self.inactive = inactive

    def set_electrons(self, electrons: int) -> None:
        self.electrons = electrons

    def set_occ(self, occ: np.ndarray) -> None:
        self.Occ = occ

    def set_coeff(self, coeff: np.ndarray) -> None:
        self.coeff = coeff

    def set_atoms_R(self, atoms_R: np.ndarray) -> None:
        self.atoms_R = atoms_R

    def set_atoms_Charge(self, atoms_Charge: List) -> None:
        self.atoms_Charge = atoms_Charge

    def set_atoms_Name(self, atoms_Name: List[str]) -> None:
        self.atoms_Name = atoms_Name

    def set_atoms(self, atoms_R: np.ndarray, atoms_Charge: List, atoms_Name: List[str]) -> None:
        self.atoms_R = atoms_R
        self.atoms_Charge = atoms_Charge
        self.atoms_Name = atoms_Name

    def set_bonds(self, bonds: List) -> None:
        self.bonds = bonds

    def get_spherical(self) -> bool:
        return self.spherical

    def get_nb(self) -> int:
        return self.nb

    def get_nAtoms(self) -> int:
        return self.nAtoms

    def get_inactive(self) -> int:
        return self.inactive

    def get_electrons(self) -> int:
        return self.electrons

    def get_occ(self) -> Optional[np.ndarray]:
        return self.Occ

    def get_coeff(self) -> Optional[np.ndarray]:
        return self.coeff

    def get_atoms_R(self) -> Optional[np.ndarray]:
        return self.atoms_R

    def get_atoms_Charge(self) -> Optional[List]:
        return self.atoms_Charge

    def get_atoms_Name(self) -> Optional[List[str]]:
        return self.atoms_Name

    def get_bonds(self) -> Optional[List]:
        return self.bonds

    def get_basis(self) -> Optional[List]:
        return self.basis

    def get_basis_norm(self) -> Optional[List]:
        return self.basis_norm
