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

### Quick Start (Recommended)

Use the automated setup script for Conda + UV hybrid environment:

```bash
git clone https://github.com/hemiku/Visualization.git
cd Visualization
bash setup_environment.sh
```

This creates a conda environment with Mayavi and installs the package with UV.

### Manual Installation

See [UV_SETUP.md](UV_SETUP.md) for detailed installation options including:
- Pure UV installation (if Mayavi wheels available)
- Conda + UV hybrid (recommended)
- System site-packages approach

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
vis.plot_Geometry(
    plot_atoms=True,
    atom_names=True,
    plot_bonds=True
)
```

Or use the verification script:

```bash
conda activate viz-env
python show_water.py  # Opens 3D interactive window
```

## Testing

Run the test suite:

```bash
# All tests
pytest tests/ -v

# Quick tests only (skip slow tests)
pytest tests/ -v -k "not slow"

# Specific test category
pytest -m unit -v          # Unit tests
pytest -m integration -v   # Integration tests
pytest -m validation -v    # Scientific validation
pytest -m visual -v        # Visual output tests
```

See [TESTING.md](TESTING.md) for comprehensive testing documentation.

## Supported Input Formats

- **Dalton**: Output files (.out) and tar archives (.tar.gz)
- **Molpro**: Output files for SAPT calculations
- **GAMMCOR**: EERPA-GVB results (.txt files, SAPTVIS binary)
- **SAPTVIS**: Binary files for SAPT(CAS) calculations

## Documentation

- [UV_SETUP.md](UV_SETUP.md) - Installation and environment setup
- [TESTING.md](TESTING.md) - Testing guide and test categories
- `example_pl.ipynb` - Jupyter notebook with detailed examples

## Requirements

- Python ≥ 3.8
- NumPy ≥ 1.20
- SciPy ≥ 1.6
- Mayavi (for 3D visualization)
- Optional: CuPy ≥ 10.0 (for GPU acceleration)

## Development

### Running Tests

```bash
conda activate viz-env
pytest tests/ -v
```

### Package Structure

```
visualization/
├── visualization.py      # Main visualization class
├── dispersion_plot.py   # Dispersion calculations
├── inputs.py            # Input file parsers
├── input_molpro.py      # Molpro-specific parsers
├── orbitals.py          # Orbital generation
├── grid.py              # 3D grid handling
├── geminals.py          # Geminal calculations
├── molecular_system.py  # Molecular data structures
└── utils.py             # Utilities
```

### Test Categories

- **Unit tests** (`tests/unit/`): Individual component testing
- **Integration tests** (`tests/integration/`): Workflow testing
- **Validation tests** (`tests/validation/`): Scientific correctness
- **Visual tests** (`tests/visual/`): Rendering validation

## License

Specified in LICENSE file.

## Citation

If you use this tool in your research, please cite:

[Citation information to be added]

## Contributing

Contributions are welcome! Please ensure:
- All tests pass: `pytest tests/`
- Code follows existing style
- New features include tests
- Update documentation as needed

## Acknowledgments

Developed for visualizing quantum chemistry calculations from Dalton, Molpro, and GAMMCOR.
