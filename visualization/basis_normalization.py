"""Basis set normalization utilities for Gaussian-type orbitals.

This module provides normalization functions for contracted Gaussian
basis functions used in quantum chemistry calculations.

The normalization power follows: pow_val = (3 + 2*l) / 2
where l is the angular momentum quantum number:
- S-type: l=0, pow_val=3/2
- P-type: l=1, pow_val=5/2
- D-type: l=2, pow_val=7/2
- F-type: l=3, pow_val=9/2
- G-type: l=4, pow_val=11/2
- H-type: l=5, pow_val=13/2
"""

import numpy as np


def normalization_summation(data: np.ndarray, pow_val: float) -> np.ndarray:
    """Calculate normalization sum for contracted Gaussians.

    Args:
        data: Array where data[:,0] are exponents, data[:,1:] are coefficients
        pow_val: Power exponent (3/2 for S, 5/2 for P, etc.)

    Returns:
        Normalization factor(s)
    """
    if len(np.shape(data)) == 2:
        n_orb = np.shape(data)[-1] - 1
        n_expans = np.shape(data)[0]
        norm = np.zeros(n_orb, dtype=np.float64)

        for i in range(n_expans):
            for j in range(n_expans):
                norm += data[i, 1:] * data[j, 1:] / (data[i, 0] + data[j, 0]) ** pow_val
    else:
        norm = data[1] * data[1] / (data[0] + data[0]) ** pow_val

    return norm


def norm_s(data: np.ndarray) -> np.ndarray:
    """Normalize S-type (l=0) orbital."""
    pow_val = 3.0 / 2.0
    norm = np.pi ** (3.0 / 2.0) * normalization_summation(data, pow_val)
    return 1 / np.sqrt(norm)


def norm_p(data: np.ndarray) -> np.ndarray:
    """Normalize P-type (l=1) orbital."""
    pow_val = 5.0 / 2.0
    fact = 1.0 / 2.0
    norm = fact * np.pi ** (3.0 / 2.0) * normalization_summation(data, pow_val)
    return 1 / np.sqrt(norm)


def norm_d(data: np.ndarray) -> np.ndarray:
    """Normalize D-type (l=2) orbital."""
    pow_val = 7.0 / 2.0
    fact = 1.0 / 4.0
    norm = fact * np.pi ** (3.0 / 2.0) * normalization_summation(data, pow_val)
    return 1 / np.sqrt(norm)


def norm_f(data: np.ndarray) -> np.ndarray:
    """Normalize F-type (l=3) orbital."""
    pow_val = 9.0 / 2.0
    fact = 15.0 / 8.0
    norm = normalization_summation(data, pow_val) * np.pi ** (3.0 / 2.0) * fact
    return 1 / np.sqrt(norm)


def norm_g(data: np.ndarray) -> np.ndarray:
    """Normalize G-type (l=4) orbital."""
    pow_val = 11.0 / 2.0
    fact = 105.0 / 16.0
    norm = normalization_summation(data, pow_val) * np.pi ** (3.0 / 2.0) * fact
    return 1 / np.sqrt(norm)


def norm_h(data: np.ndarray) -> np.ndarray:
    """Normalize H-type (l=5) orbital."""
    pow_val = 13.0 / 2.0
    fact = 945.0 / 32.0
    norm = normalization_summation(data, pow_val) * np.pi ** (3.0 / 2.0) * fact
    return 1 / np.sqrt(norm)


def calc_norm_from_basis(basis: list) -> list:
    """Calculate normalization factors for all basis functions.

    Args:
        basis: List of atom basis sets, where each atom has orbital types
               (S, P, D, F, G) as sublists

    Returns:
        List of normalization factors with same structure as input
    """
    if not basis:
        return []

    norm_funcs = [norm_s, norm_p, norm_d, norm_f, norm_g, norm_h]
    basis_norm = []

    for atom_basis in basis:
        atom_basis_norm = []
        basis_norm.append(atom_basis_norm)

        for j, orbital_type_basis in enumerate(atom_basis):
            if j < len(norm_funcs):
                atom_basis_norm.append(norm_funcs[j](orbital_type_basis))

    return basis_norm
