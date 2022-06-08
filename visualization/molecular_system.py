
class Molecular_System:

    import numpy as np

    spherical = None
    nb:int 
    nAtoms = None
    inactive = None
    electrons = None
    Occ = None
    Coeff:np.ndarray
    atoms_R:np.ndarray
    atoms_Charge: list 
    atoms_Name: list 
    bonds = None

    basis = None
    basis_norm = None
    basis_norm2 = None

    G_coeff:np.ndarray
    Orb2Gem:np.ndarray
    n_geminals:int

    AOs:np.ndarray
    MOs:np.ndarray
    geminals:np.ndarray

    def __init__(self, spherical=None, nb=None, nAtoms=None, inactive=None, electrons=None, Occ=None, Coeff=None, atoms_R=None, atoms_Charge=None, atoms_Name=None, bonds=None):

        if spherical is not None:
            self.spherical = spherical

        if nb is not None:
            self.nb = nb

        if nAtoms is not None:
            self.nAtoms = nAtoms

        if inactive is not None:
            self.inactive = inactive

        if electrons is not None:
            self.electrons = electrons

        if Occ is not None:
            self.Occ = Occ

        if Coeff is not None:
            self.Coeff = Coeff

        if atoms_R is not None:
            self.atoms_R = atoms_R

        if atoms_Charge is not None:
            self.atoms_Charge = atoms_Charge

        if atoms_Name is not None:
            self.atoms_Name = atoms_Name

        if bonds  is not None:
            self.bonds  = bonds


    def set_spherical(self, spherical):
        self.spherical = spherical

    def set_nb(self, nb):
        self.nb = nb

    def set_nAtoms(self, nAtoms):
        self.nAtoms = nAtoms

    def set_inactive(self, inactive):
        self.inactive = inactive

    def set_electrons(self, electrons):
        self.electrons = electrons

    def set_Occ(self, Occ):
        self.Occ = Occ

    def set_Coeff(self, Coeff):
        self.Coeff = Coeff

    def set_atoms_R(self, atoms_R):
        self.atoms_R = atoms_R

    def set_atoms_Charge(self, atoms_Charge):
        self.atoms_Charge = atoms_Charge

    def set_atoms_Name(self, atoms_Name):
        self.atoms_Name = atoms_Name

    def set_atoms(self, atoms_R, atoms_Charge, atoms_Name):

        self.set_atoms_R(atoms_R)
        self.set_atoms_Charge(atoms_Charge)
        self.set_atoms_Name(atoms_Name)

    def set_bonds(self, bonds):
        self.bonds = bonds

    def get_spherical(self):
        return self.spherical

    def get_nb(self):
        return self.nb

    def get_nAtoms(self):
        return self.nAtoms

    def get_inactive(self):
        return self.inactive

    def get_electrons(self):
        return self.electrons

    def get_Occ(self):
        return self.Occ

    def get_coeff(self):
        return self.Coeff

    def get_atoms_R(self):
        return self.atoms_R

    def get_atoms_Charge(self):
        return self.atoms_Charge

    def get_atoms_Name(self):
        return self.atoms_Name

    def get_bonds(self):
        return self.bonds

    def get_basis(self):
        return self.basis

    def get_basis_norm(self):
        return self.basis_norm

    def set_basis_and_norms(self, input):

        self.basis = input[0]
        self.basis_norm = input[1]
        self.basis_norm2 = input[2]
