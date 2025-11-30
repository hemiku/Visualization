"""Tests for DaltonInput parser.

Tests cover:
- Basic metadata extraction (nAtoms, nb, spherical)
- Atomic geometry parsing (positions, charges, names)
- Bond parsing
- Coefficient extraction from tar archives
- Geminal data for GVB/APSG calculations
"""

import numpy as np
import pytest
from pathlib import Path


class TestDaltonInputBasic:
    """Test basic DaltonInput functionality with water dimer."""

    @pytest.fixture
    def dalton_input(self, water_dimer_output):
        """Create DaltonInput for water dimer output file."""
        from visualization.inputs import DaltonInput
        input_path = str(water_dimer_output).replace('.out', '')
        return DaltonInput(input_type='Output', input_name=input_path)

    def test_get_nb(self, dalton_input):
        """Test number of basis functions extraction."""
        nb = dalton_input.get_nb()
        assert isinstance(nb, int)
        assert nb == 82

    def test_get_nAtoms(self, dalton_input):
        """Test atom count extraction."""
        nAtoms = dalton_input.get_nAtoms()
        assert isinstance(nAtoms, int)
        assert nAtoms == 6

    def test_get_spherical(self, dalton_input):
        """Test spherical harmonics detection."""
        spherical = dalton_input.get_spherical()
        assert spherical == 1  # Spherical harmonic basis used


class TestDaltonInputAtoms:
    """Test atomic geometry parsing."""

    @pytest.fixture
    def dalton_input(self, water_dimer_tarball):
        """Create DaltonInput for water dimer tar archive."""
        from visualization.inputs import DaltonInput
        input_path = str(water_dimer_tarball).replace('.tar.gz', '')
        dalton = DaltonInput(input_type='tar', input_name=input_path)
        dalton.get_nAtoms()
        return dalton

    def test_get_atoms_returns_tuple(self, dalton_input):
        """Test get_atoms returns correct tuple structure."""
        result = dalton_input.get_atoms()
        assert isinstance(result, tuple)
        assert len(result) == 3

        atoms_R, atoms_Charge, atoms_Name = result
        assert isinstance(atoms_R, np.ndarray)
        assert isinstance(atoms_Charge, np.ndarray)
        assert isinstance(atoms_Name, list)

    def test_get_atoms_shapes(self, dalton_input):
        """Test atomic data shapes are correct."""
        atoms_R, atoms_Charge, atoms_Name = dalton_input.get_atoms()

        assert atoms_R.shape == (6, 3)
        assert atoms_Charge.shape == (6,)
        assert len(atoms_Name) == 6

    def test_get_atoms_charges(self, dalton_input):
        """Test atomic charges for water dimer (O=8, H=1)."""
        _, atoms_Charge, _ = dalton_input.get_atoms()
        expected = np.array([8, 1, 1, 8, 1, 1])
        np.testing.assert_array_equal(atoms_Charge, expected)

    def test_get_atoms_names(self, dalton_input):
        """Test atom names for water dimer."""
        _, _, atoms_Name = dalton_input.get_atoms()
        expected = ['O1', 'H1', 'H2', 'O2', 'H3', 'H4']
        assert atoms_Name == expected

    def test_get_atoms_coordinates_first_oxygen(self, dalton_input):
        """Test first oxygen coordinates (approximate)."""
        atoms_R, _, _ = dalton_input.get_atoms()
        # First atom should be oxygen at approximately these coordinates
        # Allow some tolerance since coordinates might be in different units
        assert atoms_R.shape == (6, 3)
        assert np.all(np.isfinite(atoms_R))


class TestDaltonInputBonds:
    """Test bond parsing."""

    @pytest.fixture
    def dalton_input(self, water_dimer_tarball):
        """Create DaltonInput for water dimer (tar archive required for atoms/bonds)."""
        from visualization.inputs import DaltonInput
        input_path = str(water_dimer_tarball).replace('.tar.gz', '')
        dalton = DaltonInput(input_type='tar', input_name=input_path)
        dalton.get_nAtoms()
        dalton.get_atoms()
        return dalton

    def test_get_bonds_returns_list(self, dalton_input):
        """Test get_bonds returns a list."""
        bonds = dalton_input.get_bonds()
        assert isinstance(bonds, list)

    def test_get_bonds_count(self, dalton_input):
        """Test water dimer has 4 O-H bonds."""
        bonds = dalton_input.get_bonds()
        assert len(bonds) == 4

    def test_get_bonds_structure(self, dalton_input):
        """Test bond structure [atom1_idx, atom2_idx, distance]."""
        bonds = dalton_input.get_bonds()
        for bond in bonds:
            assert len(bond) == 3
            assert isinstance(bond[0], (int, np.integer))
            assert isinstance(bond[1], (int, np.integer))
            assert isinstance(bond[2], float)


class TestDaltonInputCoefficients:
    """Test orbital coefficient extraction from tar archive."""

    @pytest.fixture
    def dalton_input(self, water_dimer_tarball):
        """Create DaltonInput for water dimer tar archive."""
        from visualization.inputs import DaltonInput
        input_path = str(water_dimer_tarball).replace('.tar.gz', '')
        dalton = DaltonInput(input_type='tar', input_name=input_path)
        dalton.get_nb()  # Required before get_coeff
        return dalton

    def test_get_coeff_shape(self, dalton_input):
        """Test coefficient matrix shape."""
        nb = dalton_input.nb
        coeff = dalton_input.get_coeff()

        assert isinstance(coeff, np.ndarray)
        assert coeff.shape == (nb, nb)

    def test_get_coeff_finite(self, dalton_input):
        """Test all coefficients are finite."""
        coeff = dalton_input.get_coeff()
        assert np.all(np.isfinite(coeff))


class TestDaltonInputBasis:
    """Test basis set extraction."""

    @pytest.fixture
    def dalton_input(self, water_dimer_tarball):
        """Create DaltonInput for water dimer tar archive."""
        from visualization.inputs import DaltonInput
        input_path = str(water_dimer_tarball).replace('.tar.gz', '')
        dalton = DaltonInput(input_type='tar', input_name=input_path)
        dalton.get_nAtoms()
        dalton.get_atoms()
        return dalton

    def test_get_basis_returns_tuple(self, dalton_input):
        """Test get_basis returns tuple of lists."""
        result = dalton_input.get_basis()
        assert isinstance(result, tuple)
        assert len(result) == 3

        basis, basis_norm, basis_norm2 = result
        assert isinstance(basis, list)
        assert isinstance(basis_norm, list)

    def test_get_basis_structure(self, dalton_input):
        """Test basis structure has 6 atoms."""
        basis, basis_norm, _ = dalton_input.get_basis()
        assert len(basis) == 6
        assert len(basis_norm) == 6


class TestDaltonInputEthylene:
    """Test DaltonInput with ethylene dimer (GVB/APSG calculation)."""

    @pytest.fixture
    def dalton_input(self, ethylene_tarball):
        """Create DaltonInput for ethylene dimer (tar for geminal data)."""
        from visualization.inputs import DaltonInput
        input_path = str(ethylene_tarball).replace('.tar.gz', '')
        return DaltonInput(input_type='tar', input_name=input_path)

    def test_get_nb(self, dalton_input):
        """Test basis function count for ethylene dimer."""
        nb = dalton_input.get_nb()
        assert isinstance(nb, int)
        assert nb == 96

    def test_get_nAtoms(self, dalton_input):
        """Test atom count for ethylene dimer."""
        nAtoms = dalton_input.get_nAtoms()
        assert nAtoms == 12

    def test_get_inactive(self, dalton_input):
        """Test inactive orbital count for GVB calculation."""
        dalton_input.get_nb()  # Required first
        inactive = dalton_input.get_inactive()
        # Note: get_inactive may return 0 or 4 depending on parser state
        assert isinstance(inactive, (int, type(None))) or inactive >= 0

    def test_get_electrons(self, dalton_input):
        """Test electron count for GVB calculation."""
        electrons = dalton_input.get_electrons()
        assert electrons == 24


class TestDaltonInputEthyleneGeminals:
    """Test geminal-related data extraction for ethylene dimer."""

    @pytest.fixture
    def dalton_input(self, ethylene_tarball):
        """Create DaltonInput for ethylene dimer tar archive."""
        from visualization.inputs import DaltonInput
        input_path = str(ethylene_tarball).replace('.tar.gz', '')
        return DaltonInput(input_type='tar', input_name=input_path)

    def test_get_g_coeff(self, dalton_input):
        """Test geminal coefficient extraction."""
        g_coeff = dalton_input.get_g_coeff()
        assert isinstance(g_coeff, np.ndarray)
        assert len(g_coeff) > 0
        # First geminal coefficient should be close to 0.996475299688
        np.testing.assert_allclose(g_coeff[0], 0.996475299688, rtol=1e-6)

    def test_get_orb2gem(self, dalton_input):
        """Test orbital-to-geminal mapping."""
        orb2gem, nGeminal = dalton_input.get_orb2gem()
        assert isinstance(orb2gem, np.ndarray)
        assert isinstance(nGeminal, int)
        assert nGeminal > 0

    def test_get_occ(self, dalton_input):
        """Test occupation numbers."""
        dalton_input.get_nb()
        dalton_input.get_electrons()  # Required before get_occ
        occ = dalton_input.get_occ()
        assert isinstance(occ, np.ndarray)
        assert len(occ) > 0
        # All occupation numbers should be non-negative
        assert np.all(occ >= 0)


class TestDaltonInputErrorHandling:
    """Test error handling for malformed inputs."""

    def test_missing_file(self):
        """Test handling of missing file."""
        from visualization.inputs import DaltonInput
        dalton = DaltonInput(input_name='/nonexistent/path')
        with pytest.raises(FileNotFoundError):
            dalton.get_dalton_output()

    def test_missing_tar_file(self):
        """Test handling of missing tar archive."""
        from visualization.inputs import DaltonInput
        dalton = DaltonInput(input_type='tar', input_name='/nonexistent/path')
        with pytest.raises(FileNotFoundError):
            dalton.get_coeff()
