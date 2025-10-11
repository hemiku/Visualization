# UV Setup Guide for Visualization Project

## The Mayavi Challenge

Mayavi has complex dependencies (VTK, Qt) that can be difficult to build from source. Since you're using Anaconda, here's the recommended hybrid approach:

## Recommended Setup: Conda + UV

### Option 1: Install Mayavi via Conda, Everything Else via UV

```bash
# 1. Create conda environment with mayavi
conda create -n viz-env python=3.12 mayavi -c conda-forge
conda activate viz-env

# 2. Install UV in the conda environment
pip install uv

# 3. Install the project with UV (excluding mayavi which is already installed)
uv pip install -e . --no-deps  # Install without dependencies
uv pip install numpy scipy      # Install the other core deps

# 4. Install dev dependencies
uv pip install pytest pytest-cov pytest-xdist black flake8 mypy
```

### Option 2: Pure UV with System Mayavi

If you already have mayavi installed in your conda base or another environment:

```bash
# Use UV with system site packages
uv venv --system-site-packages
source .venv/bin/activate
uv pip install -e .
```

### Option 3: Try UV-only (May Fail on Some Systems)

```bash
# This might work if pre-built wheels are available for your platform
uv sync --extra dev

# If it fails building mayavi, use Option 1 or 2 instead
```

## Current Status

The `pyproject.toml` has been configured with:
- **Core deps**: numpy, scipy only (mayavi excluded to avoid build issues)
- **Optional extras**:
  - `[gpu]` - CuPy for GPU acceleration
  - `[dev]` - Testing and code quality tools

## Installing Mayavi Manually

If UV sync fails to build mayavi, install it separately:

```bash
# Via conda (recommended)
conda install mayavi -c conda-forge

# Or via pip with pre-built wheels
pip install mayavi

# Then sync the rest with UV
uv sync --extra dev
```

## Testing Your Setup

```bash
# Test that mayavi is available
python -c "from mayavi import mlab; print('Mayavi OK')"

# Test the visualization package
uv run pytest tests/unit/test_smoke.py -v
```

## Why This Approach?

1. **Mayavi build complexity**: Requires VTK compilation which can take 30+ minutes and frequently fails
2. **Conda excels at mayavi**: Pre-built binaries with all dependencies
3. **UV excels at Python packages**: Fast, reliable for pure Python and simple C extensions
4. **Best of both worlds**: Use each tool for what it does best

## Updating Dependencies

```bash
# Update UV-managed packages
uv sync --upgrade

# Update conda packages
conda update mayavi vtk -c conda-forge
```
