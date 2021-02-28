
class Molecular_System:

    Spherical = None
    nb = None
    nAtoms = None
    inactive = None
    electrons = None
    Occ = None
    Coeff = None
    Atoms_R = None
    Atoms_Charge = None
    Atoms_Name = None
    Bonds = None



    AOs = None
    MOs = None
    geminals = None

    def __init__(self, Spherical=None, nb=None, nAtoms=None, inactive=None, electrons=None, Occ=None, Coeff=None, Atoms_R=None, Atoms_Charge=None, Atoms_Name=None, Bonds=None):

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

        if Atoms_R is not None:
            self.Atoms_R = Atoms_R

        if Atoms_Charge is not None:
            self.Atoms_Charge = Atoms_Charge

        if Atoms_Name is not None:
            self.Atoms_Name = Atoms_Name

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

    def set_Atoms_R(self, Atoms_R):
        self.Atoms_R = Atoms_R

    def set_Atoms_Charge(self, Atoms_Charge):
        self.Atoms_Charge = Atoms_Charge

    def set_Atoms_Name(self, Atoms_Name):
        self.Atoms_Name = Atoms_Name

    def set_Atoms(self, Atoms_R, Atoms_Charge, Atoms_Name):

        self.set_Atoms_R(Atoms_R)
        self.set_Atoms_Charge(Atoms_Charge)
        self.set_Atoms_Name(Atoms_Name)

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

    def get_Coeff(self):
        return self.Coeff

    def get_Atoms_R(self):
        return self.Atoms_R

    def get_Atoms_Charge(self):
        return self.Atoms_Charge

    def get_Atoms_Name(self):
        return self.Atoms_Name

    def get_Bonds(self):
        return self.Bonds



    def set_basis_and_norms(self, input):

        self.basis = input[0]
        self.basis_norm = input[1]
        self.basis_norm2 = input[2]
