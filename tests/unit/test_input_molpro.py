"""Tests for MolproInput and MolproSaptInput parsers.

Tests cover:
- Basic metadata extraction (nAtoms, nb, spherical)
- Atomic geometry parsing (positions, charges, names)
- Bond parsing
- SAPT monomer selection
"""

import numpy as np
import pytest
from pathlib import Path


class TestMolproInputBasic:
    """Test basic MolproInput functionality with benzene-cyclohexane."""

    @pytest.fixture
    def molpro_input(self, benzene_output):
        """Create MolproInput for benzene-cyclohexane."""
        from visualization.input_molpro import MolproInput
        input_path = str(benzene_output).replace('.out', '')
        return MolproInput(input_name=input_path)

    def test_get_nb(self, molpro_input):
        """Test number of basis functions extraction."""
        nb = molpro_input.get_nb()
        assert isinstance(nb, int)
        assert nb == 874

    def test_get_nAtoms(self, molpro_input):
        """Test atom count extraction."""
        nAtoms = molpro_input.get_nAtoms()
        assert isinstance(nAtoms, int)
        assert nAtoms == 27

    def test_get_spherical(self, molpro_input):
        """Test spherical harmonics detection."""
        spherical = molpro_input.get_spherical()
        # May return 0 or 1 depending on parser implementation
        assert spherical in [0, 1]


class TestMolproInputAtoms:
    """Test atomic geometry parsing."""

    @pytest.fixture
    def molpro_input(self, benzene_output):
        """Create MolproInput for benzene-cyclohexane."""
        from visualization.input_molpro import MolproInput
        input_path = str(benzene_output).replace('.out', '')
        molpro = MolproInput(input_name=input_path)
        molpro.get_nAtoms()
        return molpro

    def test_get_atoms_returns_tuple(self, molpro_input):
        """Test get_atoms returns correct tuple structure."""
        result = molpro_input.get_atoms()
        assert isinstance(result, tuple)
        assert len(result) == 3

        atoms_R, atoms_Charge, atoms_Name = result
        assert isinstance(atoms_R, np.ndarray)
        assert isinstance(atoms_Charge, np.ndarray)
        assert isinstance(atoms_Name, list)

    def test_get_atoms_shapes(self, molpro_input):
        """Test atomic data shapes are correct."""
        atoms_R, atoms_Charge, atoms_Name = molpro_input.get_atoms()

        assert atoms_R.shape == (27, 3)
        assert atoms_Charge.shape == (27,)
        assert len(atoms_Name) == 27

    def test_get_atoms_charges_contain_carbon(self, molpro_input):
        """Test atomic charges include carbon (Z=6)."""
        _, atoms_Charge, _ = molpro_input.get_atoms()
        assert 6 in atoms_Charge  # Carbon

    def test_get_atoms_charges_contain_hydrogen(self, molpro_input):
        """Test atomic charges include hydrogen (Z=1)."""
        _, atoms_Charge, _ = molpro_input.get_atoms()
        assert 1 in atoms_Charge  # Hydrogen

    def test_get_atoms_coordinates_finite(self, molpro_input):
        """Test all coordinates are finite."""
        atoms_R, _, _ = molpro_input.get_atoms()
        assert np.all(np.isfinite(atoms_R))


class TestMolproInputBonds:
    """Test bond parsing."""

    @pytest.fixture
    def molpro_input(self, benzene_output):
        """Create MolproInput for benzene-cyclohexane."""
        from visualization.input_molpro import MolproInput
        input_path = str(benzene_output).replace('.out', '')
        molpro = MolproInput(input_name=input_path)
        molpro.get_nAtoms()
        molpro.get_atoms()
        return molpro

    def test_get_bonds_returns_list(self, molpro_input):
        """Test get_bonds returns a list."""
        bonds = molpro_input.get_bonds()
        assert isinstance(bonds, list)

    def test_get_bonds_not_empty(self, molpro_input):
        """Test benzene-cyclohexane has bonds."""
        bonds = molpro_input.get_bonds()
        assert len(bonds) > 0

    def test_get_bonds_structure(self, molpro_input):
        """Test bond structure [atom1_idx, atom2_idx, distance]."""
        bonds = molpro_input.get_bonds()
        for bond in bonds:
            assert len(bond) == 3
            assert isinstance(bond[0], (int, np.integer))
            assert isinstance(bond[1], (int, np.integer))
            assert isinstance(bond[2], float)


class TestMolproInputBasis:
    """Test basis set extraction."""

    @pytest.fixture
    def molpro_input(self, benzene_output):
        """Create MolproInput for benzene-cyclohexane."""
        from visualization.input_molpro import MolproInput
        input_path = str(benzene_output).replace('.out', '')
        molpro = MolproInput(input_name=input_path)
        molpro.get_nAtoms()
        molpro.get_atoms()
        return molpro

    def test_get_basis_returns_tuple(self, molpro_input):
        """Test get_basis returns tuple of lists."""
        result = molpro_input.get_basis()
        assert isinstance(result, tuple)
        assert len(result) == 3

        basis, basis_norm, basis_norm2 = result
        assert isinstance(basis, list)
        assert isinstance(basis_norm, list)

    def test_get_basis_not_empty(self, molpro_input):
        """Test basis has entries."""
        basis, _, _ = molpro_input.get_basis()
        assert len(basis) > 0


class TestMolproInputElectrons:
    """Test electron count extraction."""

    @pytest.fixture
    def molpro_input(self, benzene_output):
        """Create MolproInput for benzene-cyclohexane."""
        from visualization.input_molpro import MolproInput
        input_path = str(benzene_output).replace('.out', '')
        return MolproInput(input_name=input_path)

    def test_get_electrons(self, molpro_input):
        """Test electron count (nuclear charge) extraction."""
        electrons = molpro_input.get_electrons()
        assert isinstance(electrons, int)
        assert electrons == 42  # Benzene nuclear charge


class TestMolproSaptInput:
    """Test MolproSaptInput for SAPT calculations."""

    @pytest.fixture
    def molpro_sapt_input(self, benzene_output):
        """Create MolproSaptInput for benzene-cyclohexane SAPT."""
        from visualization.input_molpro import MolproSaptInput
        input_path = str(benzene_output).replace('.out', '')
        return MolproSaptInput(input_name=input_path)

    def test_default_monomer_is_zero(self, molpro_sapt_input):
        """Test default monomer selection is 0 (monomer A)."""
        assert molpro_sapt_input.monomer == 0

    def test_monomer_a_electrons(self, molpro_sapt_input):
        """Test monomer A electrons (nuclear charge = 42)."""
        molpro_sapt_input.monomer = 0
        electrons = molpro_sapt_input.get_electrons()
        assert electrons == 42

    def test_monomer_b_electrons(self, molpro_sapt_input):
        """Test monomer B electrons (nuclear charge = 40)."""
        molpro_sapt_input.monomer = 1
        molpro_sapt_input.output = None  # Reset cached output
        electrons = molpro_sapt_input.get_electrons()
        assert electrons == 40

    def test_monomer_a_nb(self, molpro_sapt_input):
        """Test monomer A basis function count."""
        molpro_sapt_input.monomer = 0
        nb = molpro_sapt_input.get_nb()
        assert nb == 874

    def test_monomer_b_nb(self, molpro_sapt_input):
        """Test monomer B basis function count."""
        molpro_sapt_input.monomer = 1
        molpro_sapt_input.output = None
        nb = molpro_sapt_input.get_nb()
        assert nb == 874


class TestMolproInputErrorHandling:
    """Test error handling for malformed inputs."""

    def test_missing_file(self):
        """Test handling of missing file."""
        from visualization.input_molpro import MolproInput
        molpro = MolproInput(input_name='/nonexistent/path')
        with pytest.raises(FileNotFoundError):
            molpro.get_output()


class TestMolproInputOutput:
    """Test output retrieval."""

    @pytest.fixture
    def molpro_input(self, benzene_output):
        """Create MolproInput for benzene-cyclohexane."""
        from visualization.input_molpro import MolproInput
        input_path = str(benzene_output).replace('.out', '')
        return MolproInput(input_name=input_path)

    def test_get_output(self, molpro_input):
        """Test get_output returns string."""
        output = molpro_input.get_output()
        assert isinstance(output, str)
        assert len(output) > 0

    def test_get_output_contains_seward(self, molpro_input):
        """Test output contains SEWARD section."""
        output = molpro_input.get_output()
        # The output should be split at SEWARD
        assert 'ATOMIC COORDINATES' in output or 'BASIS DATA' in output
