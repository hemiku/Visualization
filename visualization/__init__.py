"""Quantum chemistry visualization package.

This package provides tools for visualizing molecular orbitals, atomic orbitals,
geminals, and dispersion interactions from quantum chemistry calculations.

Supported input formats:
- Dalton: Output files (.out) and tar archives (.tar.gz)
- Molpro: Output files for SAPT calculations
- Molden: Standard .molden format files
"""

from visualization.molecular_system import MolecularSystem
from visualization.grid import Grid
from visualization.inputs import Input, DaltonInput, MoldenInput
from visualization.input_molpro import MolproInput, MolproSaptInput
from visualization.orbitals import OrbitalsGenerator
from visualization.geminals import GeminalGenerator

# Import Visualization (uses PyVista by default, Mayavi optional)
from visualization.visualization import Visualization

__all__ = [
    'Visualization',
    'MolecularSystem',
    'Grid',
    'Input',
    'DaltonInput',
    'MoldenInput',
    'MolproInput',
    'MolproSaptInput',
    'OrbitalsGenerator',
    'GeminalGenerator',
]

__version__ = '0.0.1'
