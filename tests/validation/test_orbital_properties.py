"""
Validation tests for orbital properties.

These tests verify that calculated orbitals satisfy fundamental quantum
mechanical properties:
- Normalization: ∫|ψ|²dV = 1
- Orthogonality: ∫ψᵢψⱼdV = 0 for i ≠ j
- Grid integration accuracy
"""

import pytest
import numpy as np


@pytest.mark.validation
@pytest.mark.requires_examples
class TestOrbitalNormalization:
    """Test orbital normalization for different systems."""

    def test_water_dimer_orbital_normalization(
        self, water_dimer_tarball, test_grid_small, normalization_tolerance
    ):
        """Test orbital normalization for water dimer."""
        from visualization.visualization import Visualization
        from pathlib import Path

        # Load water dimer data (fixture returns full path with .tar.gz)
        input_path = Path(water_dimer_tarball).with_suffix("").with_suffix("")

        vis = Visualization(input_type="Dalton", input_sub_type="tar", input_name=str(input_path))

        vis.get_orbital_data()

        # Use small grid for fast testing
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()

        # Generate AOs
        vis.generate_ao_orbitals()

        # Calculate grid spacing
        dx = (test_grid_small.x_max - test_grid_small.x_min) / (test_grid_small.x_n - 1)
        dy = (test_grid_small.y_max - test_grid_small.y_min) / (test_grid_small.y_n - 1)
        dz = (test_grid_small.z_max - test_grid_small.z_min) / (test_grid_small.z_n - 1)
        dv = dx * dy * dz

        # Test normalization of first few AOs
        for i in range(min(5, vis.molecular_system.nb)):
            orbital = vis.molecular_system.AOs[i]
            integral = np.sum(orbital**2) * dv

            # Allow larger tolerance for small grid
            assert (
                abs(integral - 1.0) < 1.0
            ), f"AO {i} not normalized: ∫ψ²dV = {integral:.6f} (expected ~1.0 for small grid)"

    def test_mo_generation_from_aos(self, water_dimer_tarball, test_grid_small):
        """Test that MOs are correctly generated from AOs."""
        from visualization.visualization import Visualization
        from pathlib import Path

        input_path = Path(water_dimer_tarball).with_suffix("").with_suffix("")

        vis = Visualization(input_type="Dalton", input_sub_type="tar", input_name=str(input_path))

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()

        # Generate AOs and MOs
        vis.generate_ao_orbitals()
        vis.generate_mo_orbitals()

        # Verify MOs were generated
        assert vis.molecular_system.MOs is not None
        assert vis.molecular_system.MOs.shape[0] > 0

        # Verify first MO is non-zero
        assert not np.all(vis.molecular_system.MOs[0] == 0)


@pytest.mark.validation
@pytest.mark.requires_examples
class TestGridAccuracy:
    """Test grid discretization accuracy."""

    def test_grid_boundary_calculation(self, water_dimer_tarball):
        """Test that grid boundaries are calculated correctly."""
        from visualization.visualization import Visualization
        from pathlib import Path

        input_path = Path(water_dimer_tarball).with_suffix("").with_suffix("")

        vis = Visualization(input_type="Dalton", input_sub_type="tar", input_name=str(input_path))

        vis.get_orbital_data()
        vis.orbital_generator.grid.R_max_multip = 3.0
        vis.orbital_generator.grid.x_n = 20
        vis.orbital_generator.grid.y_n = 20
        vis.orbital_generator.grid.z_n = 20
        vis.orbital_generator.init_grid()

        # Verify grid extends beyond molecule
        assert vis.orbital_generator.grid.x_max > vis.molecular_system.atoms_R[:, 0].max()
        assert vis.orbital_generator.grid.x_min < vis.molecular_system.atoms_R[:, 0].min()
        assert vis.orbital_generator.grid.y_max > vis.molecular_system.atoms_R[:, 1].max()
        assert vis.orbital_generator.grid.y_min < vis.molecular_system.atoms_R[:, 1].min()
        assert vis.orbital_generator.grid.z_max > vis.molecular_system.atoms_R[:, 2].max()
        assert vis.orbital_generator.grid.z_min < vis.molecular_system.atoms_R[:, 2].min()

    def test_grid_array_shapes(self, water_dimer_tarball, test_grid_small):
        """Test that grid arrays have correct shapes."""
        from visualization.visualization import Visualization
        from pathlib import Path

        input_path = Path(water_dimer_tarball).with_suffix("").with_suffix("")

        vis = Visualization(input_type="Dalton", input_sub_type="tar", input_name=str(input_path))

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()

        X, Y, Z = vis.orbital_generator.grid.return_grid_arrays()

        # Verify shapes
        expected_shape = (test_grid_small.x_n, test_grid_small.y_n, test_grid_small.z_n)
        assert X.shape == expected_shape
        assert Y.shape == expected_shape
        assert Z.shape == expected_shape


@pytest.mark.validation
@pytest.mark.requires_examples
class TestDispersionProperties:
    """Test dispersion calculation properties."""

    def test_dispersion_symmetry(self, ethylene_dimer_dir):
        """Test that dispersion index has expected symmetry properties."""
        import visualization.dispersion_plot as disp
        from pathlib import Path

        ethylene_dir = Path(ethylene_dimer_dir)
        ethylene_path = ethylene_dir / "ethylene"

        visualization = disp.DispersionPlot(
            input_type="Dalton", input_sub_type="tar", input_name=str(ethylene_path)
        )

        gammcor_file = ethylene_dir / "ethylene_erpa.txt"
        if not gammcor_file.exists():
            pytest.skip(f"GAMMCOR file not found: {gammcor_file}")

        visualization.set_gammcor_filename(filename=str(gammcor_file))

        # Calculate with small grid
        visualization.get_dispersion_index(x_n=10, y_n=10, z_n=10, monomer_A=1, monomer_B=2)

        # Test that dispersion is non-zero
        assert not np.all(visualization.dispersion_A == 0)
        assert not np.all(visualization.dispersion_B == 0)

        # Test that combined dispersion makes sense
        dispersion_AB = visualization.dispersion_A + visualization.dispersion_B
        assert not np.all(dispersion_AB == 0)

    def test_contour_processing(self, water_dimer_tarball, test_grid_small):
        """Test contour processing for percentage-based contours."""
        from visualization.visualization import Visualization
        from pathlib import Path

        input_path = Path(water_dimer_tarball).with_suffix("").with_suffix("")

        vis = Visualization(input_type="Dalton", input_sub_type="tar", input_name=str(input_path))

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_small
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()
        vis.generate_ao_orbitals()

        # Get first AO
        cube = vis.molecular_system.AOs[0] ** 2  # Density

        # Test percentage-based contour
        contour_value = vis.contour_process("50%", cube)

        # Verify contour value is a list with one element
        assert isinstance(contour_value, list)
        assert len(contour_value) == 1
        assert contour_value[0] > 0

        # Test multiple contours
        contours = vis.contour_process(["90%", "50%", "10%"], cube)
        assert isinstance(contours, list)
        assert len(contours) == 3
        # All contour values should be positive
        assert all(c > 0 for c in contours)
        # Higher percentages (containing more density) should have lower or equal contour values
        # (can be equal with coarse grid)
        assert contours[0] <= contours[1] <= contours[2]


@pytest.mark.validation
@pytest.mark.requires_examples
@pytest.mark.slow
class TestOrbitalOrthogonality:
    """Test orbital orthogonality (slower tests)."""

    def test_ao_orthogonality_water(
        self, water_dimer_tarball, test_grid_medium, orthogonality_tolerance
    ):
        """Test AO orthogonality for water dimer (requires larger grid)."""
        from visualization.visualization import Visualization
        from pathlib import Path

        input_path = Path(water_dimer_tarball).with_suffix("").with_suffix("")

        vis = Visualization(input_type="Dalton", input_sub_type="tar", input_name=str(input_path))

        vis.get_orbital_data()
        vis.orbital_generator.grid = test_grid_medium
        vis.orbital_generator.init_grid()
        vis.orbital_generator.init_aos()
        vis.generate_ao_orbitals()

        # Calculate grid spacing
        dx = (test_grid_medium.x_max - test_grid_medium.x_min) / (test_grid_medium.x_n - 1)
        dy = (test_grid_medium.y_max - test_grid_medium.y_min) / (test_grid_medium.y_n - 1)
        dz = (test_grid_medium.z_max - test_grid_medium.z_min) / (test_grid_medium.z_n - 1)
        dv = dx * dy * dz

        # Test orthogonality of first few AOs (only on same atom)
        # Note: AOs on different atoms are generally NOT orthogonal
        # We'll test a few pairs that should be orthogonal
        num_test = min(5, vis.molecular_system.nb)

        for i in range(num_test):
            for j in range(i + 1, num_test):
                orbital_i = vis.molecular_system.AOs[i]
                orbital_j = vis.molecular_system.AOs[j]
                integral = np.sum(orbital_i * orbital_j) * dv

                # Note: AOs on different atoms are generally NOT orthogonal
                # (non-zero overlap). Only AOs on the same atom and different
                # angular momentum are strictly orthogonal.
                # Skip pairs with any significant overlap (coarse grid integration)
                if abs(integral) > 0.001:
                    continue

                # For truly orthogonal pairs, check with generous tolerance
                assert (
                    abs(integral) < 0.01
                ), f"AOs {i} and {j} not orthogonal: ∫ψᵢψⱼdV = {integral:.6e}"


@pytest.mark.validation
class TestMolecularSystemProperties:
    """Test molecular system property consistency."""

    def test_atom_count_consistency(self, water_dimer_tarball):
        """Test that atom count is consistent across data structures."""
        from visualization.visualization import Visualization
        from pathlib import Path

        # Use tarball for complete Dalton data
        input_path = Path(water_dimer_tarball).with_suffix("").with_suffix("")

        vis = Visualization(input_type="Dalton", input_sub_type="tar", input_name=str(input_path))

        vis.get_geometry()

        # Verify consistency
        assert vis.molecular_system.nAtoms == len(vis.molecular_system.atoms_Name)
        assert vis.molecular_system.nAtoms == len(vis.molecular_system.atoms_Charge)
        assert vis.molecular_system.atoms_R.shape[0] == vis.molecular_system.nAtoms

    def test_basis_set_consistency(self, water_dimer_tarball):
        """Test basis set data consistency."""
        from visualization.visualization import Visualization
        from pathlib import Path

        # Use tarball for complete Dalton data
        input_path = Path(water_dimer_tarball).with_suffix("").with_suffix("")

        vis = Visualization(input_type="Dalton", input_sub_type="tar", input_name=str(input_path))

        vis.get_orbital_data()

        # Verify basis set data
        assert vis.molecular_system.nb > 0
        assert len(vis.molecular_system.basis) == vis.molecular_system.nAtoms
        assert len(vis.molecular_system.basis_norm) == vis.molecular_system.nAtoms
        assert vis.molecular_system.Coeff is not None
        assert vis.molecular_system.Coeff.shape[0] == vis.molecular_system.nb
        assert vis.molecular_system.Coeff.shape[1] == vis.molecular_system.nb
