"""
Unit tests for orbital calculation functionality.

Tests verify:
- Grid initialization and setup
- AO (Atomic Orbital) generation
- MO (Molecular Orbital) generation
- Shape consistency and basic properties
"""

import pytest
import numpy as np
from pathlib import Path


@pytest.mark.unit
class TestGridInitialization:
    """Test grid setup and initialization."""

    def test_grid_creation_default(self):
        """Test creating grid with default parameters."""
        from visualization.grid import Grid

        grid = Grid()

        # Verify default attributes exist
        assert hasattr(grid, 'x_n')
        assert hasattr(grid, 'y_n')
        assert hasattr(grid, 'z_n')
        assert hasattr(grid, 'R_max_multip')

    def test_grid_creation_custom(self):
        """Test creating grid with custom parameters."""
        from visualization.grid import Grid

        grid = Grid(x_n=15, y_n=20, z_n=25, R_max_multip=2.5)

        assert grid.x_n == 15
        assert grid.y_n == 20
        assert grid.z_n == 25
        assert grid.R_max_multip == 2.5

    @pytest.mark.requires_examples
    def test_grid_array_shapes(self, ethylene_dimer_dir, test_grid_small):
        """Test that grid arrays have correct shapes after initialization."""
        from visualization.visualization import Visualization

        # Need molecular system to initialize grid properly
        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        vis = Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()

        X, Y, Z = vis.orbital_generator.grid.return_grid_arrays()

        expected_shape = (test_grid_small.x_n, test_grid_small.y_n, test_grid_small.z_n)

        assert X.shape == expected_shape
        assert Y.shape == expected_shape
        assert Z.shape == expected_shape


@pytest.mark.unit
@pytest.mark.requires_examples
class TestAOGeneration:
    """Test Atomic Orbital (AO) generation."""

    def test_ao_generation_ethylene(self, ethylene_dimer_dir, test_grid_small):
        """Test AO generation for ethylene with small grid."""
        from visualization.visualization import Visualization

        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        vis = Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()

        # Generate AOs
        vis.generate_ao_orbitals()

        # Verify AOs were generated
        assert vis.molecular_system.AOs is not None
        assert len(vis.molecular_system.AOs) > 0

    def test_ao_shape_consistency(self, ethylene_dimer_dir, test_grid_small):
        """Test that all AOs have consistent shape."""
        from visualization.visualization import Visualization

        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        vis = Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()
        vis.generate_ao_orbitals()

        expected_shape = (test_grid_small.x_n, test_grid_small.y_n, test_grid_small.z_n)

        # Check all AOs have same shape
        for i, ao in enumerate(vis.molecular_system.AOs):
            assert ao.shape == expected_shape, \
                f"AO {i} has shape {ao.shape}, expected {expected_shape}"

    def test_ao_not_all_zero(self, ethylene_dimer_dir, test_grid_small):
        """Test that AOs contain non-zero values."""
        from visualization.visualization import Visualization

        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        vis = Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()
        vis.generate_ao_orbitals()

        # Check first few AOs are not all zero
        for i in range(min(5, len(vis.molecular_system.AOs))):
            ao = vis.molecular_system.AOs[i]
            assert not np.all(ao == 0), f"AO {i} is all zeros"
            assert np.any(np.abs(ao) > 1e-10), f"AO {i} has no significant values"

    def test_ao_count_matches_basis(self, ethylene_dimer_dir, test_grid_small):
        """Test that number of AOs matches basis set size."""
        from visualization.visualization import Visualization

        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        vis = Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()
        vis.generate_ao_orbitals()

        # Number of AOs should equal number of basis functions
        assert len(vis.molecular_system.AOs) == vis.molecular_system.nb, \
            f"Expected {vis.molecular_system.nb} AOs, got {len(vis.molecular_system.AOs)}"


@pytest.mark.unit
@pytest.mark.requires_examples
class TestMOGeneration:
    """Test Molecular Orbital (MO) generation."""

    def test_mo_generation_from_aos(self, ethylene_dimer_dir, test_grid_small):
        """Test MO generation from AOs."""
        from visualization.visualization import Visualization

        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        vis = Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()

        # Generate AOs then MOs
        vis.generate_ao_orbitals()
        vis.generate_mo_orbitals()

        # Verify MOs were generated
        assert vis.molecular_system.MOs is not None
        assert len(vis.molecular_system.MOs) > 0

    def test_mo_shape_consistency(self, ethylene_dimer_dir, test_grid_small):
        """Test that all MOs have consistent shape."""
        from visualization.visualization import Visualization

        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        vis = Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()
        vis.generate_ao_orbitals()
        vis.generate_mo_orbitals()

        expected_shape = (test_grid_small.x_n, test_grid_small.y_n, test_grid_small.z_n)

        # Check all MOs have same shape
        for i, mo in enumerate(vis.molecular_system.MOs):
            assert mo.shape == expected_shape, \
                f"MO {i} has shape {mo.shape}, expected {expected_shape}"

    def test_mo_not_all_zero(self, ethylene_dimer_dir, test_grid_small):
        """Test that MOs contain non-zero values."""
        from visualization.visualization import Visualization

        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        vis = Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()
        vis.generate_ao_orbitals()
        vis.generate_mo_orbitals()

        # Check first few MOs are not all zero
        for i in range(min(5, len(vis.molecular_system.MOs))):
            mo = vis.molecular_system.MOs[i]
            assert not np.all(mo == 0), f"MO {i} is all zeros"
            assert np.any(np.abs(mo) > 1e-10), f"MO {i} has no significant values"

    def test_mo_count_matches_ao_count(self, ethylene_dimer_dir, test_grid_small):
        """Test that number of MOs equals number of AOs."""
        from visualization.visualization import Visualization

        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        vis = Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()
        vis.generate_ao_orbitals()
        vis.generate_mo_orbitals()

        # Number of MOs should equal number of AOs/basis functions
        assert len(vis.molecular_system.MOs) == len(vis.molecular_system.AOs), \
            f"MO count ({len(vis.molecular_system.MOs)}) != AO count ({len(vis.molecular_system.AOs)})"


@pytest.mark.unit
@pytest.mark.requires_examples
class TestOrbitalCalculationWorkflow:
    """Test complete orbital calculation workflows."""

    def test_full_workflow_tiny_grid(self, ethylene_dimer_dir):
        """Test complete workflow with tiny grid (very fast)."""
        from visualization.visualization import Visualization
        from visualization.grid import Grid

        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        vis = Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        # Use tiny grid for speed
        tiny_grid = Grid(x_n=5, y_n=5, z_n=5, R_max_multip=2.0)

        vis.get_orbital_data()
        vis.orbital_generator.grid = tiny_grid
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()
        vis.generate_ao_orbitals()
        vis.generate_mo_orbitals()

        # Verify complete workflow succeeded
        assert vis.molecular_system.AOs is not None
        assert vis.molecular_system.MOs is not None
        assert len(vis.molecular_system.AOs) > 0
        assert len(vis.molecular_system.MOs) > 0

    def test_orbitals_method_shortcut(self, ethylene_dimer_dir):
        """Test the get_orbitals() convenience method."""
        from visualization.visualization import Visualization

        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        vis = Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        # Use convenience method (includes all steps)
        vis.get_orbitals(x_n=5, y_n=5, z_n=5, R_max_multip=2.0, gpu=False)

        # Verify orbitals generated
        assert vis.molecular_system.AOs is not None
        assert vis.molecular_system.MOs is not None
        assert len(vis.molecular_system.AOs) > 0
        assert len(vis.molecular_system.MOs) > 0


@pytest.mark.unit
@pytest.mark.requires_examples
@pytest.mark.slow
class TestOrbitalPerformance:
    """Test orbital calculation performance."""

    def test_calculation_completes_in_reasonable_time(self, ethylene_dimer_dir):
        """Test that orbital calculation completes within reasonable time."""
        import time
        from visualization.visualization import Visualization

        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        vis = Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        # Small grid should complete in < 10 seconds
        start_time = time.time()
        vis.get_orbitals(x_n=10, y_n=10, z_n=10, gpu=False)
        elapsed = time.time() - start_time

        assert elapsed < 10.0, \
            f"Orbital calculation took {elapsed:.2f}s (expected <10s for 10x10x10 grid)"

        print(f"Orbital calculation time: {elapsed:.2f}s")
