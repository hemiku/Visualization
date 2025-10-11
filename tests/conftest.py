"""
Pytest configuration and shared fixtures for molecular visualization tests.

This module provides fixtures for:
- Example data files (water dimer, ethylene, benzene-cyclohexane)
- Molecular systems
- Test grids
- Reference outputs
"""

import pytest
import numpy as np
from pathlib import Path


# ==============================================================================
# Path Fixtures
# ==============================================================================

@pytest.fixture(scope="session")
def project_root():
    """Return path to project root directory."""
    return Path(__file__).parent.parent


@pytest.fixture(scope="session")
def examples_dir(project_root):
    """Return path to Examples directory."""
    return project_root / "Examples"


@pytest.fixture(scope="session")
def water_dimer_dir(examples_dir):
    """Return path to water dimer example."""
    return examples_dir


@pytest.fixture(scope="session")
def ethylene_dimer_dir(examples_dir):
    """Return path to ethylene dimer example."""
    return examples_dir / "C2H4_C2H4"


@pytest.fixture(scope="session")
def benzene_cyclo_dir(examples_dir):
    """Return path to benzene-cyclohexane example."""
    return examples_dir / "benzcyclo"


@pytest.fixture(scope="session")
def benzene_cyclo_excite_dir(examples_dir):
    """Return path to benzene-cyclohexane excited state example."""
    return examples_dir / "benzcyclo_excite"


# ==============================================================================
# Data File Fixtures
# ==============================================================================

@pytest.fixture(scope="session")
def water_dimer_output(water_dimer_dir):
    """Return path to H2O_H2O.out file."""
    output_file = water_dimer_dir / "H2O_H2O.out"
    if not output_file.exists():
        pytest.skip(f"Water dimer output file not found: {output_file}")
    return output_file


@pytest.fixture(scope="session")
def water_dimer_tarball(water_dimer_dir):
    """Return path to H2O_H2O.tar.gz file."""
    tar_file = water_dimer_dir / "H2O_H2O.tar.gz"
    if not tar_file.exists():
        pytest.skip(f"Water dimer tarball not found: {tar_file}")
    return tar_file


@pytest.fixture(scope="session")
def ethylene_output(ethylene_dimer_dir):
    """Return path to ethylene.out file."""
    output_file = ethylene_dimer_dir / "ethylene.out"
    if not output_file.exists():
        pytest.skip(f"Ethylene output file not found: {output_file}")
    return output_file


@pytest.fixture(scope="session")
def ethylene_tarball(ethylene_dimer_dir):
    """Return path to ethylene.tar.gz file."""
    tar_file = ethylene_dimer_dir / "ethylene.tar.gz"
    if not tar_file.exists():
        pytest.skip(f"Ethylene tarball not found: {tar_file}")
    return tar_file


@pytest.fixture(scope="session")
def benzene_output(benzene_cyclo_dir):
    """Return path to benzcyclo.out file."""
    output_file = benzene_cyclo_dir / "benzcyclo.out"
    if not output_file.exists():
        pytest.skip(f"Benzene output file not found: {output_file}")
    return output_file


@pytest.fixture(scope="session")
def benzene_sapt_output(benzene_cyclo_dir):
    """Return path to sapt.out file."""
    output_file = benzene_cyclo_dir / "sapt.out"
    if not output_file.exists():
        pytest.skip(f"SAPT output file not found: {output_file}")
    return output_file


@pytest.fixture(scope="session")
def benzene_saptvis_binary(benzene_cyclo_dir):
    """Return path to SAPTVIS binary file."""
    binary_file = benzene_cyclo_dir / "SAPTVIS"
    if not binary_file.exists():
        pytest.skip(f"SAPTVIS binary not found: {binary_file}")
    return binary_file


# ==============================================================================
# Numerical Tolerance Fixtures
# ==============================================================================

@pytest.fixture
def default_tolerance():
    """Default numerical tolerance for comparisons."""
    return 1e-10


@pytest.fixture
def relaxed_tolerance():
    """Relaxed tolerance for less precise comparisons."""
    return 1e-6


@pytest.fixture
def normalization_tolerance():
    """Tolerance for orbital normalization checks (∫ψ²dV = 1)."""
    return 1e-3  # 0.1% tolerance


@pytest.fixture
def orthogonality_tolerance():
    """Tolerance for orbital orthogonality checks (∫ψᵢψⱼdV = 0)."""
    return 1e-6


# ==============================================================================
# Grid Fixtures
# ==============================================================================

@pytest.fixture
def test_grid_small():
    """Small test grid for fast tests (10x10x10)."""
    from visualization.grid import Grid
    grid = Grid(R_max_multip=3.0, x_n=10, y_n=10, z_n=10)
    return grid


@pytest.fixture
def test_grid_medium():
    """Medium test grid (30x30x30)."""
    from visualization.grid import Grid
    grid = Grid(R_max_multip=3.0, x_n=30, y_n=30, z_n=30)
    return grid


# ==============================================================================
# Helper Functions
# ==============================================================================

@pytest.fixture
def assert_normalized():
    """Return function to check if orbital is normalized."""
    def _check(orbital, grid_spacing, tolerance=1e-3):
        """Check if ∫|ψ|²dV ≈ 1"""
        dv = grid_spacing[0] * grid_spacing[1] * grid_spacing[2]
        integral = np.sum(orbital**2) * dv
        assert abs(integral - 1.0) < tolerance, \
            f"Orbital not normalized: ∫ψ²dV = {integral:.6f} (expected 1.0)"
        return True
    return _check


@pytest.fixture
def assert_orthogonal():
    """Return function to check if two orbitals are orthogonal."""
    def _check(orbital1, orbital2, grid_spacing, tolerance=1e-6):
        """Check if ∫ψᵢψⱼdV ≈ 0"""
        dv = grid_spacing[0] * grid_spacing[1] * grid_spacing[2]
        integral = np.sum(orbital1 * orbital2) * dv
        assert abs(integral) < tolerance, \
            f"Orbitals not orthogonal: ∫ψᵢψⱼdV = {integral:.6e} (expected ~0)"
        return True
    return _check


# ==============================================================================
# Skip Markers
# ==============================================================================

def pytest_configure(config):
    """Add custom markers."""
    config.addinivalue_line(
        "markers", "requires_gpu: mark test as requiring GPU/CuPy"
    )


def pytest_collection_modifyitems(config, items):
    """Automatically skip GPU tests if CuPy not available."""
    try:
        import cupy
        gpu_available = True
    except ImportError:
        gpu_available = False

    if not gpu_available:
        skip_gpu = pytest.mark.skip(reason="CuPy not available - GPU tests skipped")
        for item in items:
            if "gpu" in item.keywords:
                item.add_marker(skip_gpu)
