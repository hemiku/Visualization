
class Molecular_System:

    Spherical = None
    nb = None
    nAtoms = None
    inactive = None
    electrons = None
    Occ = None
    Coeff = None
    atoms_R = None
    atoms_Charge = None
    atoms_Name = None
    Bonds = None

    basis = None
    basis_norm = None
    basis_norm2 = None

    AOs = None
    MOs = None
    geminals = None

    def __init__(self, Spherical=None, nb=None, nAtoms=None, inactive=None, electrons=None, Occ=None, Coeff=None, atoms_R=None, atoms_Charge=None, atoms_Name=None, Bonds=None):

        if Spherical is not None:
            self.Spherical = Spherical

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

        if Bonds  is not None:
            self.Bonds  = Bonds


    def set_Spherical(self, Spherical):
        self.Spherical = Spherical

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

    def set_Bonds(self, Bonds):
        self.Bonds = Bonds

    def get_Spherical(self):
        return self.Spherical

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

    def get_Bonds(self):
        return self.Bonds

    def get_basis(self):
        return self.basis

    def get_basis_norm(self):
        return self.basis_norm

    def set_basis_and_norms(self, input):

        self.basis = input[0]
        self.basis_norm = input[1]
        self.basis_norm2 = input[2]
