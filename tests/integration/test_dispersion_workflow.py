"""
Integration tests for dispersion plotting workflow.

Based on example_pl.ipynb - tests the complete workflow:
1. Load Dalton tar files
2. Load GAMMCOR results
3. Calculate dispersion index
4. Visualize dispersion
"""

import pytest
import numpy as np
from pathlib import Path


@pytest.mark.integration
@pytest.mark.requires_examples
class TestDispersionWorkflow:
    """Test complete dispersion plotting workflow."""

    def test_load_ethylene_dispersion_data(self, ethylene_dimer_dir):
        """Test loading ethylene dimer data for dispersion calculation."""
        import visualization.dispersion_plot as disp

        # Use pathlib for path handling
        ethylene_path = Path(ethylene_dimer_dir) / 'ethylene'

        # Load Dalton data
        visualization = disp.DispersionPlot(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        # Verify visualization object created
        assert visualization is not None
        assert hasattr(visualization, 'molecular_system')

    def test_ethylene_dispersion_calculation_small_grid(self, ethylene_dimer_dir):
        """Test dispersion calculation with small grid (fast test)."""
        import visualization.dispersion_plot as disp

        # Use pathlib for all path operations
        ethylene_dir = Path(ethylene_dimer_dir)
        ethylene_path = ethylene_dir / 'ethylene'
        gammcor_file = ethylene_dir / 'ethylene_erpa.txt'

        visualization = disp.DispersionPlot(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        # Check if GAMMCOR file exists before trying to load
        if not gammcor_file.exists():
            pytest.skip(f"GAMMCOR file not found: {gammcor_file}")

        visualization.set_GAMMCOR_filename(filename=str(gammcor_file))

        # Calculate dispersion with small grid for speed
        visualization.get_dispersion_index(
            x_n=10,
            y_n=10,
            z_n=10,
            monomer_A=1,
            monomer_B=2
        )

        # Verify dispersion data was calculated
        assert hasattr(visualization, 'dispersion_A')
        assert hasattr(visualization, 'dispersion_B')
        assert visualization.dispersion_A is not None
        assert visualization.dispersion_B is not None

    @pytest.mark.slow
    def test_ethylene_dispersion_calculation_medium_grid(self, ethylene_dimer_dir):
        """Test dispersion calculation with medium grid (slower, more accurate)."""
        import visualization.dispersion_plot as disp

        ethylene_dir = Path(ethylene_dimer_dir)
        ethylene_path = ethylene_dir / "ethylene"

        visualization = disp.DispersionPlot(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        gammcor_file = Path(ethylene_dimer_dir) / "ethylene_erpa.txt"
        if not gammcor_file.exists():
            pytest.skip(f"GAMMCOR file not found: {gammcor_file}")

        visualization.set_GAMMCOR_filename(filename=str(gammcor_file))

        # Calculate with medium grid (as in notebook: 50x50x50)
        visualization.get_dispersion_index(
            x_n=50,
            y_n=50,
            z_n=50,
            monomer_A=1,
            monomer_B=2
        )

        # Verify grid dimensions
        assert visualization.dispersion_A.shape == (50, 50, 50)
        assert visualization.dispersion_B.shape == (50, 50, 50)


@pytest.mark.integration
@pytest.mark.requires_examples
class TestGeminalWorkflow:
    """Test geminal visualization workflow."""

    def test_geminal_data_loading(self, ethylene_dimer_dir):
        """Test loading data for geminal visualization."""
        import visualization.dispersion_plot as disp

        ethylene_dir = Path(ethylene_dimer_dir)
        ethylene_path = ethylene_dir / "ethylene"

        visualization = disp.DispersionPlot(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        # Verify molecular system initialized
        assert visualization.molecular_system is not None

    def test_geminal_calculation_small_grid(self, ethylene_dimer_dir):
        """Test geminal calculation with small grid."""
        import visualization.dispersion_plot as disp

        ethylene_dir = Path(ethylene_dimer_dir)
        ethylene_path = ethylene_dir / "ethylene"

        visualization = disp.DispersionPlot(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(ethylene_path)
        )

        gammcor_file = Path(ethylene_dimer_dir) / "ethylene_erpa.txt"
        if not gammcor_file.exists():
            pytest.skip(f"GAMMCOR file not found: {gammcor_file}")

        visualization.set_GAMMCOR_filename(filename=str(gammcor_file))

        # Small grid for fast testing
        visualization.get_dispersion_index(
            x_n=10,
            y_n=10,
            z_n=10,
            monomer_A=1,
            monomer_B=2,
            gpu=False
        )

        # Verify geminal data available
        assert hasattr(visualization.molecular_system, 'geminals')


@pytest.mark.integration
@pytest.mark.requires_examples
class TestSAPTWorkflow:
    """Test SAPT(CAS) workflow from notebook."""

    def test_load_molpro_sapt_data(self, benzene_cyclo_dir, benzene_output):
        """Test loading Molpro SAPT data."""
        import visualization.visualization as V

        # benzene_output fixture returns Path with .out extension
        # Remove .out since MolproSapt adds it
        benzene_path = Path(benzene_output).with_suffix('')

        vis_A = V.Visualization(
            input_type='MolproSapt',
            input_sub_type='output',
            input_name=str(benzene_path)
        )

        vis_B = V.Visualization(
            input_type='MolproSapt',
            input_sub_type='output',
            input_name=str(benzene_path)
        )
        vis_B.data_input.monomer = 1

        # Verify both monomers loaded
        assert vis_A is not None
        assert vis_B is not None
        assert vis_B.data_input.monomer == 1

    def test_sapt_geometry_loading(self, benzene_output):
        """Test SAPT geometry loading for both monomers."""
        import visualization.visualization as V

        benzene_path = Path(benzene_output).with_suffix('')

        vis_A = V.Visualization(
            input_type='MolproSapt',
            input_sub_type='output',
            input_name=str(benzene_path)
        )

        vis_A.get_geometry()

        # Verify geometry loaded
        assert vis_A.molecular_system.nAtoms > 0
        assert vis_A.molecular_system.atoms_R is not None
        assert vis_A.molecular_system.atoms_Charge is not None

    def test_sapt_orbital_data_loading(self, benzene_output):
        """Test SAPT orbital data loading."""
        import visualization.visualization as V

        benzene_path = Path(benzene_output).with_suffix('')

        vis_A = V.Visualization(
            input_type='MolproSapt',
            input_sub_type='output',
            input_name=str(benzene_path)
        )

        vis_A.get_geometry()
        vis_A.get_orbital_data()

        # Verify orbital data loaded
        assert vis_A.orbital_generator is not None
        assert vis_A.molecular_system.nb > 0  # Number of basis functions
        assert vis_A.molecular_system.Coeff is not None

    def test_sapt_grid_initialization_small(self, benzene_output):
        """Test SAPT grid initialization with small grid."""
        import visualization.visualization as V

        benzene_path = Path(benzene_output).with_suffix('')

        vis_A = V.Visualization(
            input_type='MolproSapt',
            input_sub_type='output',
            input_name=str(benzene_path)
        )

        vis_A.get_geometry()
        vis_A.get_orbital_data()

        # Set small grid for fast testing
        vis_A.orbital_generator.grid.R_max_multip = 1.5
        vis_A.orbital_generator.grid.x_n = 10
        vis_A.orbital_generator.grid.y_n = 10
        vis_A.orbital_generator.grid.z_n = 10
        vis_A.orbital_generator.init_grid()
        vis_A.orbital_generator.init_AOs()
        vis_A.orbital_generator.spherical = 1

        # Verify grid initialized
        assert vis_A.orbital_generator.grid.x_n == 10
        assert vis_A.orbital_generator.grid.y_n == 10
        assert vis_A.orbital_generator.grid.z_n == 10
        assert vis_A.orbital_generator.AOs is not None

    def test_sapt_saptvis_file_reading(self, benzene_cyclo_dir, benzene_saptvis_binary):
        """Test reading SAPTVIS binary file."""
        import visualization.utils as utils

        ANBasis, BNBasis, NOccupA, NOccupB, A_Occ, B_Occ, ACMO, BCMO, Qmat = \
            utils.read_SAPTVIS(str(benzene_saptvis_binary))

        # Verify SAPTVIS data loaded
        assert ANBasis > 0
        assert BNBasis > 0
        assert NOccupA > 0
        assert NOccupB > 0
        assert A_Occ is not None
        assert B_Occ is not None
        assert ACMO is not None
        assert BCMO is not None
        assert Qmat is not None

    @pytest.mark.slow
    def test_sapt_dispersion_calculation_small(self, benzene_output, benzene_saptvis_binary):
        """Test SAPT dispersion calculation with small grid."""
        import visualization.visualization as V
        import visualization.utils as utils
        import numpy as np

        benzene_path = Path(benzene_output).with_suffix('')

        # Load data for monomer A
        vis_A = V.Visualization(
            input_type='MolproSapt',
            input_sub_type='output',
            input_name=str(benzene_path)
        )

        # Load data for monomer B
        vis_B = V.Visualization(
            input_type='MolproSapt',
            input_sub_type='output',
            input_name=str(benzene_path)
        )
        vis_B.data_input.monomer = 1

        # Get geometry and orbital data
        vis_A.get_geometry()
        vis_B.get_geometry()
        vis_A.get_orbital_data()
        vis_B.get_orbital_data()

        # Initialize small grids
        for vis in [vis_A, vis_B]:
            vis.orbital_generator.grid.R_max_multip = 1.5
            vis.orbital_generator.grid.x_n = 10
            vis.orbital_generator.grid.y_n = 10
            vis.orbital_generator.grid.z_n = 10
            vis.orbital_generator.init_grid()
            vis.orbital_generator.init_AOs()
            vis.orbital_generator.spherical = 1

        # Read SAPTVIS data
        ANBasis, BNBasis, NOccupA, NOccupB, A_Occ, B_Occ, ACMO, BCMO, Qmat = \
            utils.read_SAPTVIS(str(benzene_saptvis_binary))

        # Set coefficients
        vis_A.molecular_system.Coeff[:, :] = ACMO
        vis_B.molecular_system.Coeff[:, :] = BCMO

        # Generate orbitals (CPU version for compatibility)
        vis_A.generate_AO_orbitals()
        vis_A.generate_MO_orbitals()
        vis_B.generate_AO_orbitals()
        vis_B.generate_MO_orbitals()

        # Calculate dispersion
        sapt_dispersion_A = np.zeros_like(vis_A.molecular_system.AOs[0])
        sapt_dispersion_B = np.zeros_like(vis_B.molecular_system.AOs[0])

        Qmat_tmp = np.reshape(Qmat, [NOccupB, NOccupA]).transpose()

        for i in range(NOccupA):
            for j in range(NOccupB):
                sapt_dispersion_A += Qmat_tmp[i, j] * A_Occ[i] * vis_A.molecular_system.MOs[i]**2
                sapt_dispersion_B += Qmat_tmp[i, j] * B_Occ[j] * vis_B.molecular_system.MOs[j]**2

        sapt_dispersion_AB = -0.5 * (sapt_dispersion_A + sapt_dispersion_B)

        # Verify dispersion calculated
        assert sapt_dispersion_AB is not None
        assert sapt_dispersion_AB.shape == (10, 10, 10)
        assert not np.all(sapt_dispersion_AB == 0)  # Should have non-zero values
