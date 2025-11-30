# GVT - General Visualization Tool

Quantum chemistry visualization tool for molecular orbitals, dispersion interactions, and SAPT calculations.

## Features

- **Molecular Visualization**: Display molecular geometries with atoms and bonds
- **Orbital Visualization**: Visualize atomic orbitals (AOs) and molecular orbitals (MOs)
- **Dispersion Analysis**: Calculate and visualize dispersion interactions from EERPA-GVB calculations
- **Geminal Visualization**: Display geminal wavefunctions
- **SAPT Visualization**: Visualize SAPT(CAS) dispersion components
- **GPU Acceleration**: Optional CUDA support via CuPy for faster calculations

## Installation

### Option A: Fresh Environment (Recommended for Development)

```bash
git clone https://github.com/hemiku/Visualization.git
cd Visualization
uv sync                    # Creates .venv with all dependencies
uv run pytest tests/ -v    # Verify installation
```

### Option B: Into Existing Conda/Venv Environment

```bash
conda activate my-env      # Or: source .venv/bin/activate
git clone https://github.com/hemiku/Visualization.git
cd Visualization
uv pip install --system -e .
pytest tests/ -v           # Verify installation
```

### With Jupyter Notebook Support (Interactive 3D)

```bash
# Option A:
uv sync --extra notebook
uv run jupyter lab

# Option B:
uv pip install --system -e ".[notebook]"
jupyter lab
```

### For Development

```bash
uv sync --extra dev
```

## Quick Example

```python
import visualization.visualization as V

# Load water dimer
vis = V.Visualization(
    input_type='Dalton',
    input_sub_type='tar',
    input_name='Examples/H2O_H2O'
)

# Visualize geometry
vis.get_geometry()
vis.plot_geometry(
    plot_atoms=True,
    atom_names=True,
    plot_bonds=True
)
```

### In Jupyter Notebooks

```python
import pyvista as pv
pv.set_jupyter_backend('trame')  # Enable interactive 3D

import visualization.visualization as V
vis = V.Visualization(...)
vis.get_geometry()
vis.plot_geometry()  # Interactive 3D in notebook!
```

## Testing

```bash
# All tests
uv run pytest tests/ -v

# Quick tests only
uv run pytest tests/ -v -m "not slow"

# Specific categories
uv run pytest -m unit -v          # Unit tests
uv run pytest -m integration -v   # Integration tests
uv run pytest -m validation -v    # Scientific validation
```

## Supported Input Formats

- **Dalton**: Output files (.out) and tar archives (.tar.gz)
- **Molpro**: Output files for SAPT calculations
- **GAMMCOR**: EERPA-GVB results (.txt files, SAPTVIS binary)
- **SAPTVIS**: Binary files for SAPT(CAS) calculations

## Requirements

- Python >= 3.9
- NumPy >= 1.20
- SciPy >= 1.6
- PyVista >= 0.38 (default visualization backend)
- Optional: CuPy >= 10.0 (GPU acceleration)
- Optional: Mayavi >= 4.7 (alternative backend, install via conda)

## Package Structure

```
visualization/
├── visualization.py      # Main visualization class
├── dispersion_plot.py    # Dispersion calculations
├── inputs.py             # Input file parsers
├── input_molpro.py       # Molpro-specific parsers
├── orbitals.py           # Orbital generation
├── grid.py               # 3D grid handling
├── geminals.py           # Geminal calculations
├── molecular_system.py   # Molecular data structures
└── backends/             # Visualization backends (pyvista, mayavi)
```

## License

MIT License - see LICENSE file.

## Contributing

Contributions are welcome! Please ensure:
- All tests pass: `uv run pytest tests/`
- Code follows existing style
- New features include tests
