#!/usr/bin/env python
# coding: utf-8


class DaltonPlt:
    import numpy as np

    def __init__(self, path):
        self.path = path

        self.Atoms_R = None
        self.Atoms_Charge = None
        self.Atoms_Name = None

        self.basis = None
        self.basis_norm = None
        self.basis_norm2 = None

        self.X = None
        self.Y = None
        self.Z = None

        f = open(path, 'r')
        self.Out = f.read()
        f.close()

    def read_data(self):
        self.record1 = self.np.fromstring(self.Out[:4], dtype=self.np.int32)[0]
        self.record2 = self.np.fromstring(self.Out[4:8], dtype=self.np.int32)[0]

        self.x_n = self.np.fromstring(self.Out[8:12], dtype=self.np.int32)[0]
        self.y_n = self.np.fromstring(self.Out[12:16], dtype=self.np.int32)[0]
        self.z_n = self.np.fromstring(self.Out[16:120], dtype=self.np.int32)[0]

        self.z_min = self.np.fromstring(self.Out[20:24], dtype=self.np.float32)[0]
        self.z_max = self.np.fromstring(self.Out[24:28], dtype=self.np.float32)[0]

        self.y_min = self.np.fromstring(self.Out[28:32], dtype=self.np.float32)[0]
        self.y_max = self.np.fromstring(self.Out[32:36], dtype=self.np.float32)[0]

        self.x_min = self.np.fromstring(self.Out[36:40], dtype=self.np.float32)[0]
        self.x_max = self.np.fromstring(self.Out[40:44], dtype=self.np.float32)[0]

        self.values = self.np.reshape(self.np.fromstring(self.Out[44:], dtype=self.np.float32, count=self.x_n * self.y_n * self.z_n), [self.x_n, self.y_n, self.z_n])

        self.XYZ_str = "  self.np.mgrid[ " + str(self.x_min) + ":" + str(self.x_max) + ":" + str(self.x_n) + "j, " + \
                       str(self.y_min) + ":" + str(self.y_max) + ":" + str(self.y_n) + "j, " + \
                       str(self.z_min) + ":" + str(self.z_max) + ":" + str(self.z_n) + "j ]"

    def get_grid(self):
        exec("X, Y ,Z = " + self.XYZ_str)

        return X, Y, Z

    def gen_grid(self):
        exec("self.X, self.Y ,self.Z = " + self.XYZ_str)

    def get_self_grid(self):
        if not ('X' in dir(self) and 'Y' in dir(self) and 'Z' in dir(self)):
            self.gen_grid()

        return self.X, self.Y, self.Z


class Grid():
    def __init__(self, x_min, x_max, x_n, y_min, y_max, y_n, z_min, z_max, z_n):
        self.x_min = x_min
        self.x_max = x_max
        self.x_n = x_n
        self.y_min = y_min
        self.y_max = y_max
        self.y_n = y_n
        self.z_min = z_min
        self.z_max = z_max
        self.z_n = z_n


class Geminal():
    def __init__(self, Grid, values):
        self.g = Grid
        self.v = values


def letters(input):
    return ''.join(filter(str.isalpha, input))


class GemialGenerator:
    import tarfile

    import numpy as np
    from mayavi import mlab

    def __init__(self, Filename):
        if Filename[-4:] == '.out':
            self.File_Name = Filename[:-4]
        else:
            self.File_Name = Filename

    def set_source(self, source):

        self.source = source

    def get_Dalton_Output(self):

        f = open(self.File_Name + ".out", 'r')
        Out = f.read()
        f.close()

        return Out

    def get_harmonic(self):

        Out = self.get_Dalton_Output()

        if (Out.find("Spherical harmonic basis used.") > 0):
            self.Spherical = 1
        else:
            self.Spherical = 0

    def get_nb(self):

        Out = self.get_Dalton_Output()

        self.nb = int(Out[Out.find("Number of basis functions"):Out.find("Number of basis functions") + 38].split()[-1])

    def get_nAtoms(self):

        Out = self.get_Dalton_Output()

        self.nAtoms = int(Out[Out.find("Total number of atoms:"):Out.find("Total number of atoms:") + 30].split()[-1])

    def get_inactive(self):

        Out = self.get_Dalton_Output()

        self.inactive = int(Out[Out.find(".INACTIVE"):Out.find("Number of basis functions") + 100].split()[1])

    def get_electrons(self):

        Out = self.get_Dalton_Output()

        self.electrons = int(Out[Out.find("@    Number of electrons in active shells"):Out.find(
            "@    Total charge of the molecule") + 100].split()[7])

    def get_Occ(self):

        self.Occ = self.np.zeros([self.nb], dtype=self.np.float64)

        Out = self.get_Dalton_Output()

        Geminal_buf = Out[Out.find("APSG geminal coefficients and natural orbital occupations:"):Out.rfind(
            'NEWORB " orbitals punched.')].split(
            "====================================================================")[1].split("\n")[1:-1]

        for i in range(self.electrons):
            self.Occ[self.inactive + i] = float(Geminal_buf[i].split()[6])

        self.Occ[:self.inactive] = 1

    def get_G_coeff(self):

        self.G_coeff = self.np.zeros([self.electrons], dtype=self.np.float64)

        Out = self.get_Dalton_Output()
        Geminal_buf = Out[Out.find("APSG geminal coefficients and natural orbital occupations:"):
                          Out.rfind('NEWORB " orbitals punched.')].split(
            "====================================================================")[1].split("\n")[1:-1]

        for i in range(self.electrons):
            self.G_coeff[i] = float(Geminal_buf[i].split()[5])

    def get_Orb2Gem(self):
        self.Orb2Gem = self.np.zeros([self.electrons], dtype=self.np.int64)

        Out = self.get_Dalton_Output()
        Geminal_buf = Out[Out.find("APSG geminal coefficients and natural orbital occupations:"):
                          Out.rfind('NEWORB " orbitals punched.')].split(
            "====================================================================")[1].split("\n")[1:-1]

        for i in range(self.electrons):
            self.Orb2Gem[i] = int(Geminal_buf[i].split()[3]) - 1

        self.nGeminal = int(len(self.Orb2Gem) / 2)

    def get_Coeff(self):

        if self.source == 'MOPUN':
            f = open("DALTON.MOPUN", 'r')
            F_MOPUN = f.read()
            f.close()

        if self.source == 'tar':
            tar = self.tarfile.open(self.File_Name + ".tar.gz")
            f = tar.extractfile(tar.getmember("DALTON.MOPUN"))
            F_MOPUN = f.read().decode(encoding='utf-8')
            tar.close()

        F_MOPUN = F_MOPUN.replace('-', ' -')

        self.Coeff = self.np.reshape(self.np.fromstring(" ".join(F_MOPUN[F_MOPUN.find("\n"):].split()),
                                                        dtype=self.np.float64, sep=' '), [self.nb, self.nb])

    def get_Basis_File(self):

        if self.source == 'MOPUN':
            f = open("DALTON.BAS", 'r')
            F_BAS = f.read()
            f.close()

        if self.source == 'tar':
            tar = self.tarfile.open(self.File_Name + ".tar.gz")
            f = tar.extractfile(tar.getmember("DALTON.BAS"))
            F_BAS = f.read().decode(encoding='utf-8')
            tar.close()

        return F_BAS

    def get_Atoms(self):

        self.Atoms_R = self.np.zeros([self.nAtoms, 3], dtype=self.np.float64)
        self.Atoms_Charge = self.np.zeros(self.nAtoms, dtype=self.np.int64)
        self.Atoms_Name = []

        F_BAS = self.get_Basis_File()

        BASIS = F_BAS[F_BAS[F_BAS.find("ATOMBASIS"):].find("\n") + F_BAS.find("ATOMBASIS") + 1:]
        BAS_split = BASIS.split("\n")[:-1]
        n = 0
        for i in range(len(BAS_split)):

            if (BAS_split[i][0] != " " and len(BAS_split[i].split()) == 4):
                self.Atoms_Name += [BAS_split[i].split()[0]]
                self.Atoms_Charge[n] = float(BAS_split[i - 1].split()[0])
                self.Atoms_R[n] = self.np.array(
                    [float(BAS_split[i].split()[1]), float(BAS_split[i].split()[2]), float(BAS_split[i].split()[3])])
                n += 1

    def get_Bonds(self):

        self.Bonds = []

        Out = self.get_Dalton_Output()

        Bonds_str = Out[Out.find("Bond distances (Angstrom):"):Out.find("Bond angles (degrees)")]
        Bonds_str = Bonds_str[Bonds_str.find("bond"):]
        if 1:
            for i in range(len(Bonds_str.split("\n")) - 3):
                self.Bonds.append([(Bonds_str.split("\n")[i]).split()[2], (Bonds_str.split("\n")[i]).split()[3],
                                   float((Bonds_str.split("\n")[i]).split()[4])])
        else:
            for i in range(len(Bonds_str.split("\n")) - 3):
                self.Bonds.append([(Bonds_str.split("\n")[i]).split()[2], (Bonds_str.split("\n")[i]).split()[4],
                                   float((Bonds_str.split("\n")[i]).split()[6])])

    def set_init_visual_data(self):

        self.Atoms_Color = {
            "H": (1, 1, 1),
            "He": (0.5, 0.5, 0.5),
            "O": (1, 0, 0),
            "B": (233.0 / 255, 165.0 / 255, 165.0 / 255),
            "N": (0, 0, 1),
            "C": (0, 0, 0),
            "F": (0, 1, 0),
            "Ne": (1.0 * 178 / 256, 1.0 * 229 / 256, 1.0 * 247 / 256)
        }

        self.Atoms_Scale = {
            "H": 1.0,
            "He": 1.0,
            "O": 1.4,
            "B": 1.4,
            "N": 1.4,
            "C": 1.4,
            "F": 1.4,
            "Ne": 1.4
        }

        self.Atoms_Mass = {
            "H": 1.00794,
            "He": 4.002602,
            "Li": 6.941,
            "Be": 9.012182,
            "B": 10.811,
            "C": 12.011,
            "N": 14.00674,
            "O": 15.9994,
            "F": 18.9984032,
            "Ne": 20.1797,
            "Na": 22.989768,
            "Mg": 24.305,
            "Al": 26.981539,
            "Si": 28.0855,
            "P": 30.973762,
            "S": 32.066,
            "Cl": 35.4527,
            "Ar": 39.948,
            "K": 39.0983,
            "Ca": 40.078,
            "Sc": 44.95591,
            "Ti": 47.88,
            "V": 50.9415,
            "Cr": 51.9961,
            "Mn": 54.93805,
            "Fe": 55.847,
            "Co": 58.9332,
            "Ni": 58.69,
            "Cu": 63.546,
            "Zn": 65.39,
            "Ga": 69.723,
            "Ge": 72.61,
            "As": 74.92159,
            "Se": 78.96,
            "Br": 79.904,
            "Kr": 83.8,
            "Rb": 85.4678,
            "Sr": 87.62,
            "Y": 88.90585,
            "Zr": 91.224,
            "Nb": 92.90638,
            "Mo": 95.94,
        }

    def Norm_S(self, C):

        NOrb = self.np.shape(C)[-1] - 1
        NExpans = self.np.shape(C)[0]

        Norm = self.np.zeros(NOrb, dtype=self.np.float64)

        for n in range(NOrb):
            for i in range(NExpans):
                Norm[n] += self.np.sqrt(2.0) * self.np.pi ** (3.0 / 2.0) * C[i, n + 1] ** 2 / (
                            4.0 * C[i, 0] ** (3.0 / 2.0))
                for j in range(i):
                    Norm[n] += 2.0 * self.np.pi ** (3.0 / 2.0) * C[i, n + 1] * C[j, n + 1] / (
                                C[i, 0] * self.np.sqrt(C[i, 0] + C[j, 0]) + C[j, 0] * self.np.sqrt(C[i, 0] + C[j, 0]))

        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_P(self, Dane):

        if (len(self.np.shape(Dane)) == 2):

            NOrb = self.np.shape(Dane)[-1] - 1
            NExpans = self.np.shape(Dane)[0]
            Norm = self.np.zeros(NOrb, dtype=self.np.float64)

            for n in range(NOrb):
                for i in range(NExpans):

                    Norm[n] += self.np.sqrt(2.0) * self.np.pi ** (3.0 / 2.0) * Dane[i, n + 1] ** 2.0 / (
                                16.0 * Dane[i, 0] ** (5.0 / 2.0))
                    for j in range(i):
                        Norm[n] += self.np.pi ** (3.0 / 2.0) * Dane[i, n + 1] * Dane[j, n + 1] / (
                                    Dane[i, 0] ** 2.0 * self.np.sqrt(Dane[i, 0] + Dane[j, 0]) + 2.0 * Dane[i, 0] * Dane[
                                j, 0] * self.np.sqrt(Dane[i, 0] + Dane[j, 0]) + Dane[j, 0] ** 2.0 * self.np.sqrt(
                                Dane[i, 0] + Dane[j, 0]))

        else:

            Norm += self.np.sqrt(2.0) * self.np.pi ** (3.0 / 2.0) * Dane[1] ** 2.0 / (16.0 * Dane[0] ** (5.0 / 2.0))

        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_D(self, Dane):
        if (len(self.np.shape(Dane)) == 2):

            NOrb = self.np.shape(Dane)[-1] - 1
            Norm = self.np.zeros(NOrb, dtype=self.np.float64)

            for n in range(NOrb):
                Norm[n] = 2 ** (3.0 / 4.0) * self.np.sqrt(3.0) * Dane[0, n] ** (7.0 / 4.0) * self.np.sqrt(
                    Dane[n, n + 1] ** (-2)) / (3.0 * self.np.pi ** (3.0 / 4.0))
        else:
            Norm = 2 ** (3.0 / 4.0) * self.np.sqrt(3.0) * Dane[0] ** (7.0 / 4.0) * self.np.sqrt(Dane[1] ** (-2)) / (
                        3.0 * self.np.pi ** (3.0 / 4.0))

        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_F(self, Dane):
        if (len(self.np.shape(Dane)) == 2):

            NOrb = self.np.shape(Dane)[-1] - 1
            Norm = self.np.zeros(NOrb, dtype=self.np.float64)

            for n in range(NOrb):
                Norm[n] = 2 ** (3.0 / 4.0) * self.np.sqrt(3.0) * Dane[0, n] ** (7.0 / 4.0) * self.np.sqrt(
                    Dane[n, n + 1] ** (-2)) / (3.0 * self.np.pi ** (3.0 / 4.0))
        else:
            Norm = 2 ** (3.0 / 4.0) * self.np.sqrt(3.0) * Dane[0] ** (7.0 / 4.0) * self.np.sqrt(Dane[1] ** (-2)) / (
                        3.0 * self.np.pi ** (3.0 / 4.0))

        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm(self, Dane, Angular):

        pow = (3.0 + 2.0 * Angular) / 4.0

        if (len(self.np.shape(Dane)) == 2):

            NOrb = self.np.shape(Dane)[-1] - 1
            NExpans = self.np.shape(Dane)[0]
            Norm = self.np.zeros(NOrb, dtype=self.np.float64)

            for n in range(NOrb):
                Norm[n] = 2 ** (3.0 / 4.0) * self.np.sqrt(3.0) * Dane[0, n] ** (7.0 / 4.0) * self.np.sqrt(
                    Dane[n, n + 1] ** (-2)) / (3.0 * self.np.pi ** (3.0 / 4.0))

        # else:


    def get_Basis_new(self):

        self.basis = []
        self.basis_norm = []
        self.basis_norm2 = []

        F_BAS = self.get_Basis_File()

        import re

        F_BAS_split_lines = F_BAS.splitlines()

        Atom_header_pattern    = ' {1,9}\d{1,9}\. '
        Atom_Geometrie_pattern = '^\w{1,6} \D \D \D'
        Basis_header_pattern   = '^H  {1,3}\d {1,4}\d$'

        Atoms_Gropus = []

        for  line in F_BAS_split_lines[7:]:

            if re.match(Atom_header_pattern, line):

                Atoms_Gropu = {'headerLine':line, 'Geometries':[], 'Basis':[] }
                Atoms_Gropus.append(Atoms_Gropu )

            elif re.match(Atom_Geometrie_pattern, line):

                Atoms_Gropu['Geometries'].append(line)

            elif re.match(Basis_header_pattern, line):
                
                basis_part = {'header':line, 'data':[] }
                Atoms_Gropu['Basis'].append(basis_part)

            else:
                basis_part['data'].append(line)

        self.map_atoms_from_BAS_file(Atoms_Gropus )


    def map_atoms_from_BAS_file(self, Atoms_Gropus ):
        """
        docstring
        """

        self.Atoms_R = self.np.zeros([self.nAtoms, 3], dtype=self.np.float64)
        self.Atoms_Charge = self.np.zeros(self.nAtoms, dtype=self.np.int64)
        self.Atoms_Name = []

        i = 0

        for Atoms_Gropu in Atoms_Gropus:

            #TODO basis 
            Orbitals = []

            for Orbital in Atoms_Gropu['Basis']:
                
                dim = self.np.fromstring(Orbital['header'][1:], dtype=self.np.int64, sep=' ')  
                dim[1] += 1

                Orbital =  self.np.reshape( self.np.fromstring( ''.join(Orbital['data']), dtype=self.np.float64, sep=' '), dim)
                Orbitals.append(Orbital)

            Atom_Charge = int(Atoms_Gropu['headerLine'][:Atoms_Gropu['headerLine'].find('.')] )

            for Atom in Atoms_Gropu['Geometries']:
                
                Atom_split = Atom.split()

                Atom_name = Atom_split[0]
                Atom_R    = self.np.array([float(Atom_split[1]), float(Atom_split[2]), float(Atom_split[3])])


                self.Atoms_Name.append( Atom_name )
                self.Atoms_Charge[i] = Atom_Charge
                self.Atoms_R[i] = Atom_R

                self.basis.append(Orbitals)

                i += 1

        for n in range(self.nAtoms):

            self.basis_norm.append([])
            self.basis_norm2.append([])
            # print(n, self.basis[n])
            for j in range(len(self.basis[n])):
                # self.basis_norm2[n].append(self.Norm_2(self.basis[n][j] , j))
                if (j == 0):
                    # self.basis_norm[n].append(self.Norm_1(self.basis[n][j] , j))
                    self.basis_norm[n].append(self.Norm_S(self.basis[n][j]))
                    self.basis_norm2[n].append(self.Norm_S2(self.basis[n][j]))
                if (j == 1):
                    self.basis_norm[n].append(self.Norm_P(self.basis[n][j]))
                    self.basis_norm2[n].append(self.Norm_P2(self.basis[n][j]))
                if (j == 2):
                    self.basis_norm[n].append(self.Norm_D(self.basis[n][j]))
                    self.basis_norm2[n].append(self.Norm_D2(self.basis[n][j]))
                if (j == 3):
                    # self.basis_norm[n].append(self.Norm_F(self.basis[n][j]))
                    self.basis_norm2[n].append(self.Norm_F2(self.basis[n][j]))
                if (j == 4):
                    # self.basis_norm[n].append(self.Norm_G(self.basis[n][j]))
                    self.basis_norm2[n].append(self.Norm_G2(self.basis[n][j]))

    def get_Basis(self):

        self.basis = []
        self.basis_norm = []
        self.basis_norm2 = []

        F_BAS = self.get_Basis_File()

        BASIS = F_BAS[F_BAS[F_BAS.find("ATOMBASIS"):].find("\n") + F_BAS.find("ATOMBASIS") + 1:]
        BAS_split = BASIS.split("\n")[:-1]

        n = 0
        for i in range(len(self.Atoms_Name)):
            self.basis.append([])

            if (i != len(self.Atoms_Name) - 1):

                buf = "\n".join(
                    BASIS[BASIS.find(self.Atoms_Name[i]):BASIS.find(self.Atoms_Name[i + 1])].split("\n")[1:-1])
                for j in range(len(buf.split("H"))):
                    if (buf.split("H")[j] != ''):
                        dim = self.np.fromstring(buf.split("H")[j].split("\n")[0], dtype=self.np.int64, sep=' ',
                                                 count=2)
                        dim[1] += 1

                        data = self.np.reshape(self.np.array(
                            self.np.fromstring(" ".join(buf.split("H")[j].split("\n")[1:]), dtype=self.np.float64,
                                               sep=' ', count=self.np.prod(dim))), dim)

                        self.basis[i].append(data)

            else:

                buf = "\n".join(BASIS[BASIS.find(self.Atoms_Name[i]):].split("\n")[1:-1])
                for j in range(len(buf.split("H"))):
                    if (buf.split("H")[j] != ''):
                        dim = self.np.fromstring(buf.split("H")[j].split("\n")[0], dtype=self.np.int64, sep=' ',
                                                 count=2)
                        dim[1] += 1

                        data = self.np.reshape(self.np.array(
                            self.np.fromstring(" ".join(buf.split("H")[j].split("\n")[1:]), dtype=self.np.float64,
                                               sep=' ', count=self.np.prod(dim))), dim)
                        self.basis[i].append(data)

        for n in range(self.nAtoms):

            self.basis_norm.append([])
            self.basis_norm2.append([])
            # print(n, self.basis[n])
            for j in range(len(self.basis[n])):
                # self.basis_norm2[n].append(self.Norm_2(self.basis[n][j] , j))
                if (j == 0):
                    # self.basis_norm[n].append(self.Norm_1(self.basis[n][j] , j))
                    self.basis_norm[n].append(self.Norm_S(self.basis[n][j]))
                    self.basis_norm2[n].append(self.Norm_S2(self.basis[n][j]))
                if (j == 1):
                    self.basis_norm[n].append(self.Norm_P(self.basis[n][j]))
                    self.basis_norm2[n].append(self.Norm_P2(self.basis[n][j]))
                if (j == 2):
                    self.basis_norm[n].append(self.Norm_D(self.basis[n][j]))
                    self.basis_norm2[n].append(self.Norm_D2(self.basis[n][j]))
                if (j == 3):
                    # self.basis_norm[n].append(self.Norm_F(self.basis[n][j]))
                    self.basis_norm2[n].append(self.Norm_F2(self.basis[n][j]))
                if (j == 4):
                    # self.basis_norm[n].append(self.Norm_G(self.basis[n][j]))
                    self.basis_norm2[n].append(self.Norm_G2(self.basis[n][j]))

    def get_data(self):

        self.set_init_visual_data()

        self.get_nb()
        self.get_harmonic()
        self.get_nAtoms()
        self.get_inactive()
        self.get_electrons()
        self.get_Occ()
        self.get_G_coeff()
        self.get_Coeff()
        self.get_Basis_new()        
        self.get_Bonds()
        self.get_Orb2Gem()

        self.set_Init_Grid_param()
        self.Prepe_Grid()
        self.Gen_Grid()

    def set_Init_Grid_param(self, R_max_multip=12.0, x_min='Initial', x_max='Initial', x_n=100,
                            y_min='Initial', y_max='Initial', y_n=100,
                            z_min='Initial', z_max='Initial', z_n=100):

        self.R_max_multip = R_max_multip

        self.x_min = x_min
        self.x_max = x_max
        self.x_n = x_n
        self.y_min = y_min
        self.y_max = y_max
        self.y_n = y_n
        self.z_min = z_min
        self.z_max = z_max
        self.z_n = z_n

    def Prepe_Grid(self):

        # import pdb; pdb.set_trace()

        Basis_max = []
        for i in range(len(self.basis)):
            for j in range(len(self.basis[i])):
                Basis_max.append(1.64526336574595 / self.np.sqrt((self.basis[i][j]).max()))

        R_max = self.R_max_multip * self.np.max(Basis_max)

        if self.x_max == 'Initial':
            self.x_max = self.np.max(self.Atoms_R[:, 0]) + R_max
        if self.x_min == 'Initial':
            self.x_min = self.np.min(self.Atoms_R[:, 0]) - R_max

        if self.y_max == 'Initial':
            self.y_max = self.np.max(self.Atoms_R[:, 1]) + R_max
        if self.y_min == 'Initial':
            self.y_min = self.np.min(self.Atoms_R[:, 1]) - R_max

        if self.z_max == 'Initial':
            self.z_max = self.np.max(self.Atoms_R[:, 2]) + R_max
        if self.z_min == 'Initial':
            self.z_min = self.np.min(self.Atoms_R[:, 2]) - R_max

        self.XYZ_str = "self.X,self.Y,self.Z =  self.np.mgrid[" + str(self.x_min) + ":" + str(self.x_max) + ":" + str(
            self.x_n) + "j, " \
                       + str(self.y_min) + ":" + str(self.y_max) + ":" + str(self.y_n) + "j, " \
                       + str(self.z_min) + ":" + str(self.z_max) + ":" + str(self.z_n) + "j ]"

    def Gen_Grid(self):

        exec(self.XYZ_str)

    def Calc_Mas_Inertia(self):

        Center_of_mass = self.np.zeros([3], dtype=self.np.float64)
        Moment_of_inertia = self.np.zeros([3, 3], dtype=self.np.float64)
        Mass = 0.0
        for n in range(self.nAtoms):
            Center_of_mass += self.Atoms_R[n, :] * self.Atoms_Mass[letters(self.Atoms_Name[n])]
            Mass += self.Atoms_Mass[letters(self.Atoms_Name[n])]

        Center_of_mass = Center_of_mass / Mass

        for n in range(self.nAtoms):
            Moment_of_inertia += self.np.outer((self.Atoms_R[n, :] - Center_of_mass),
                                               (self.Atoms_R[n, :] - Center_of_mass)) * self.Atoms_Mass[
                                     letters(self.Atoms_Name[n])]

        print("Center_of_mass", Center_of_mass)
        print("Moment_of_inertia", Moment_of_inertia)

        v, w = self.np.linalg.eig(Moment_of_inertia)
        print(v)
        print(w)

    def Gen_AO_init_array(self):

        self.AO_init = self.np.zeros(
            [self.nb, self.np.shape(self.X)[0], self.np.shape(self.Y)[1], self.np.shape(self.Z)[2]],
            dtype=self.np.float64)

    def Calc_AO_init(self):

        self.Gen_AO_init_array()

        self.Calc_AO(self.AO_init, self.X, self.Y, self.Z)

    def Calc_AO_init2(self):

        self.Gen_AO_init_array()

        self.Calc_AO_2(self.AO_init, self.X, self.Y, self.Z)

    def Calc_AO(self, AO, X, Y, Z):

        if self.Spherical:
            dx = X[1, 0, 0] - X[0, 0, 0]
            dy = Y[0, 1, 0] - Y[0, 0, 0]
            dz = Z[0, 0, 1] - Z[0, 0, 0]

            dv = dx * dy * dz

        nOrb = 0

        for n in range(self.nAtoms):

            x = X - self.Atoms_R[n, 0]
            y = Y - self.Atoms_R[n, 1]
            z = Z - self.Atoms_R[n, 2]

            R = x * x + y * y + z * z

            r = self.np.sqrt(R)
            theta = self.np.arctan(y / x)
            phi = self.np.arccos(z / r)

            for j in range(len(self.basis[n])):

                if (j == 0):

                    for k in range(self.np.shape(self.basis[n][j])[1] - 1):

                        for l in range(self.np.shape(self.basis[n][j])[0]):
                            AO[nOrb, :, :, :] += ((self.basis[n][j])[l, k + 1] * self.np.exp(
                                -(self.basis[n][j])[l, 0] * R)) ** 1
                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                        nOrb += 1

                if (j == 1):

                    for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                        for m in ['x', 'y', 'z']:
                            if (m == 'x'):
                                # print "nOrb: ", nOrb, " typ p_" +m
                                for l in range(self.np.shape(self.basis[n][j])[0]):
                                    AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp(
                                        -(self.basis[n][j])[l, 0] * R) * (x)
                                AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                nOrb += 1

                            if (m == 'y'):
                                # print "nOrb: ", nOrb, " typ p_" +m
                                for l in range(self.np.shape(self.basis[n][j])[0]):
                                    AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp(
                                        -(self.basis[n][j])[l, 0] * R) * (y)
                                AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                nOrb += 1

                            if (m == 'z'):
                                # print "nOrb: ", nOrb, " typ p_" +m
                                for l in range(self.np.shape(self.basis[n][j])[0]):
                                    AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp(
                                        -(self.basis[n][j])[l, 0] * R) * (z)
                                AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                nOrb += 1
                if (j == 2):
                    if (self.Spherical == 1):

                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in ['zz', 'xz', 'yz', 'xx-yy', 'xy']:
                                    if (m == 'zz'):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp(
                                                -(self.basis[n][j])[l, 0] * R) * (3 * z * z - R)
                                        AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt(
                                            self.np.sum(AO[nOrb, :, :, :] * AO[nOrb, :, :, :]) * dv)
                                        nOrb += 1

                                    if (m == 'xz'):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp(
                                                -(self.basis[n][j])[l, 0] * R) * x * z
                                        AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt(
                                            self.np.sum(AO[nOrb, :, :, :] * AO[nOrb, :, :, :]) * dv)
                                        nOrb += 1

                                    if (m == 'yz'):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp(
                                                -(self.basis[n][j])[l, 0] * R) * y * z
                                        AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt(
                                            self.np.sum(AO[nOrb, :, :, :] * AO[nOrb, :, :, :]) * dv)
                                        nOrb += 1

                                    if (m == 'xx-yy'):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp(
                                                -(self.basis[n][j])[l, 0] * R) * (x * x - y * y)
                                        AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt(
                                            self.np.sum(AO[nOrb, :, :, :] * AO[nOrb, :, :, :]) * dv)
                                        nOrb += 1

                                    if (m == 'xy'):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp(
                                                -(self.basis[n][j])[l, 0] * R) * (x * y)
                                        AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt(
                                            self.np.sum(AO[nOrb, :, :, :] * AO[nOrb, :, :, :]) * dv)
                                        nOrb += 1

                if (j == 3):
                    if (self.Spherical == 1):

                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in ['-3', '-2', '-1', '0', '1', '2', '3']:
                                    nOrb += 1

    def Plot_Geometry(self, Plot_Atoms=1, Atom_Names=1, Plot_Bonds=1):

        # from mayavi import mlab

        self.mlab.figure("Geometria", bgcolor=(.5, .5, .75), size=(1000, 1000))
        self.mlab.clf()

        if Plot_Atoms:
            for i in range(self.nAtoms):
                self.mlab.points3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], scale_factor=self.Atoms_Scale[letters(self.Atoms_Name[i])], resolution=20, color=self.Atoms_Color[letters(self.Atoms_Name[i])], scale_mode='none')

        if Atom_Names:
            for i in range(self.nAtoms):
                self.mlab.text3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], self.Atoms_Name[i], scale=(.9, .9, .9))

        if Plot_Bonds:
            for i in range(len(self.Bonds)):
                self.mlab.plot3d(
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 0], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 0]]),
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 1], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 1]]),
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 2], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 2]]),
                    tube_radius=0.2, color=(233.0 / 255, 165.0 / 255, 165.0 / 255))

        self.mlab.show()

    def Calc_MO_init(self):

        self.MO_init = self.np.zeros( [self.nb, self.np.shape(self.X)[0], self.np.shape(self.Y)[1], self.np.shape(self.Z)[2]], dtype=self.np.float64)

        for i in range(self.inactive + 2 * self.nGeminal):
            for j in range(self.nb):
                self.MO_init[i, :, :, :] += self.Coeff[i, j] * self.AO_init[j, :, :, :]

    def Calc_Geminals_init(self):

        self.GEMINALS_init = self.np.zeros( [self.nGeminal, self.np.shape(self.X)[0], self.np.shape(self.Y)[1], self.np.shape(self.Z)[2]], dtype=self.np.float64)

        for i in range(len(self.Orb2Gem)):
            self.GEMINALS_init[self.Orb2Gem[i], :, :, :] += ( self.G_coeff[i] * self.MO_init[i + self.inactive, :, :, :] ) ** 2

    def Plot_AO_Orbital0(self, Atom_Names=1):

        # from mayavi import mlab

        self.mlab.figure("Geometria", bgcolor=(.5, .5, .75), size=(1000, 1000))
        self.mlab.clf()

        for i in range(self.nAtoms):
            if Atom_Names:
                self.mlab.text3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], self.Atoms_Name[i], scale=(.9, .9, .9))

        self.mlab.contour3d(self.X, self.Y, self.Z, (self.AO[0]), contours=12, opacity=0.5)

        self.mlab.show()

    def Plot_MO_Orbital0(self, Atom_Names=1):

        # from mayavi import mlab

        self.mlab.figure("Geometria", bgcolor=(.5, .5, .75), size=(1000, 1000))
        self.mlab.clf()

        for i in range(self.nAtoms):
            self.mlab.points3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], scale_factor=self.Atoms_Scale[letters(self.Atoms_Name[i])], resolution=20, color=self.Atoms_Color[letters(self.Atoms_Name[i])], scale_mode='none')
            if Atom_Names:
                self.mlab.text3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], self.Atoms_Name[i],    scale=(.9, .9, .9))

        self.mlab.contour3d(self.X, self.Y, self.Z, (self.MO[0]), contours=12, opacity=0.5)

        self.mlab.show()

    def Plot_Geminal0_init(self, Atom_Names=1):

        self.Plot_Geminal_init(Atom_Names=Atom_Names)

    def Plot_Geminal_init(self, Plot_Atoms=1, Atom_Names=1, Plot_Bonds=1, Geminal_number=0, contours=12, opacity=0.5, Add_atom_number=False):

        # from mayavi import mlab

        self.mlab.figure("Geometria", bgcolor=(.5, .5, .75), size=(1000, 1000))
        self.mlab.clf()

        if Plot_Atoms:
            for i in range(self.nAtoms):
                self.mlab.points3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], scale_factor=self.Atoms_Scale[letters(self.Atoms_Name[i])], resolution=20, color=self.Atoms_Color[letters(self.Atoms_Name[i])], scale_mode='none')

        if Atom_Names:
            for i in range(self.nAtoms):                
                if Add_atom_number:
                    self.mlab.text3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], self.Atoms_Name[i]+'_'+str(i), scale=(.75, .75, .75))
                else:
                    self.mlab.text3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], self.Atoms_Name[i], scale=(.75, .75, .75))



        if Plot_Bonds:
            for i in range(len(self.Bonds)):
                self.mlab.plot3d(
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 0], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 0]]),
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 1], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 1]]),
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 2], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 2]]),
                    tube_radius=0.2, color=(233.0 / 255, 165.0 / 255, 165.0 / 255))
        for i in Geminal_number:
            self.mlab.contour3d(self.X, self.Y, self.Z, (self.GEMINALS_init[i]), contours=contours, opacity=opacity)

        self.mlab.show()

    def Plot_MO_Orbital_init(self, Plot_Atoms=1, Atom_Names=1, Plot_Bonds=1, Orbital_number=0):

        # from mayavi import mlab

        self.mlab.figure("Geometria", bgcolor=(.5, .5, .75), size=(1000, 1000))
        self.mlab.clf()

        if Plot_Atoms:
            for i in range(self.nAtoms):
                self.mlab.points3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2],
                                   scale_factor=self.Atoms_Scale[letters(self.Atoms_Name[i])], resolution=20,
                                   color=self.Atoms_Color[letters(self.Atoms_Name[i])], scale_mode='none')

        if Atom_Names:
            for i in range(self.nAtoms):
                self.mlab.text3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], self.Atoms_Name[i], scale=(.9, .9, .9))

        if Plot_Bonds:
            for i in range(len(self.Bonds)):
                self.mlab.plot3d(
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 0], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 0]]),
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 1], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 1]]),
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 2], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 2]]),
                    tube_radius=0.2, color=(233.0 / 255, 165.0 / 255, 165.0 / 255))

        self.mlab.contour3d(self.X, self.Y, self.Z, (self.MO_init[Orbital_number]), contours=12, opacity=0.5)

        self.mlab.show()

    def Plot_AO_Orbital_init(self, Plot_Atoms=1, Atom_Names=1, Plot_Bonds=1, Orbital_number=0):

        # from mayavi import mlab

        self.mlab.figure("AO_" + str(Orbital_number), bgcolor=(.5, .5, .75), size=(1000, 1000))
        self.mlab.clf()

        if Plot_Atoms:
            for i in range(self.nAtoms):
                self.mlab.points3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], scale_factor=self.Atoms_Scale[letters(self.Atoms_Name[i])], resolution=20, color=self.Atoms_Color[letters(self.Atoms_Name[i])], scale_mode='none')

        if Atom_Names:
            for i in range(self.nAtoms):
                self.mlab.text3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], self.Atoms_Name[i], scale=(.2, .2, .2))

        if Plot_Bonds:
            for i in range(len(self.Bonds)):
                self.mlab.plot3d(
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 0], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 0]]),
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 1], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 1]]),
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 2], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 2]]),
                    tube_radius=0.2, color=(233.0 / 255, 165.0 / 255, 165.0 / 255))

        self.mlab.contour3d(self.X, self.Y, self.Z, (self.AO_init[Orbital_number]), contours=12, opacity=0.5)

        self.mlab.show()

    def Geminal_Generator(self, x_n=100, y_n=100, z_n=100, fast = False, fast_threshold = 10**-4):

        self.Geminnals = []

        self.CONTUR = self.np.zeros([self.nGeminal], dtype=self.np.float64)
        Center_Geminal = self.np.zeros([self.nGeminal, 3], dtype=self.np.float64)

        x_max = self.np.zeros([self.nGeminal], dtype=self.np.float64)
        x_min = self.np.zeros([self.nGeminal], dtype=self.np.float64)

        y_max = self.np.zeros([self.nGeminal], dtype=self.np.float64)
        y_min = self.np.zeros([self.nGeminal], dtype=self.np.float64)

        z_max = self.np.zeros([self.nGeminal], dtype=self.np.float64)
        z_min = self.np.zeros([self.nGeminal], dtype=self.np.float64)

        for i in range(self.nGeminal):
            Sort_geminal = self.np.sort(self.GEMINALS_init[i, :, :, :], axis=None)
            CSum_geminal = self.np.cumsum(Sort_geminal)
            Sum_geminal = self.np.sum(Sort_geminal)
            Contur_sum = CSum_geminal < Sum_geminal * 0.1
            self.CONTUR[i] = Sort_geminal[Contur_sum][-1]
            maska = self.GEMINALS_init[i, :, :, :] > self.CONTUR[i]

            Center_Geminal[i, 0] = self.np.sum(self.GEMINALS_init[i, :, :, :] * self.X) / self.np.sum( self.GEMINALS_init[i, :, :, :])
            Center_Geminal[i, 1] = self.np.sum(self.GEMINALS_init[i, :, :, :] * self.Y) / self.np.sum( self.GEMINALS_init[i, :, :, :])
            Center_Geminal[i, 2] = self.np.sum(self.GEMINALS_init[i, :, :, :] * self.Z) / self.np.sum( self.GEMINALS_init[i, :, :, :])

            x_max[i] = 1.5 * (self.np.max(self.X[maska]) - Center_Geminal[i, 0]) + Center_Geminal[i, 0]
            x_min[i] = 1.5 * (self.np.min(self.X[maska]) - Center_Geminal[i, 0]) + Center_Geminal[i, 0]

            y_max[i] = 1.5 * (self.np.max(self.Y[maska]) - Center_Geminal[i, 1]) + Center_Geminal[i, 1]
            y_min[i] = 1.5 * (self.np.min(self.Y[maska]) - Center_Geminal[i, 1]) + Center_Geminal[i, 1]

            z_max[i] = 1.5 * (self.np.max(self.Z[maska]) - Center_Geminal[i, 2]) + Center_Geminal[i, 2]
            z_min[i] = 1.5 * (self.np.min(self.Z[maska]) - Center_Geminal[i, 2]) + Center_Geminal[i, 2]

        for i in range(self.nGeminal):

            XYZ_str = "X,Y,Z =  self.np.mgrid[" + str(x_min[i]) + ":" + str(x_max[i]) + ":" + str(x_n) + "j, " + str(
                y_min[i]) + ":" + str(y_max[i]) + ":" + str(y_n) + "j, " + str(z_min[i]) + ":" + str(
                z_max[i]) + ":" + str(z_n) + "j ]"
        

            grid = Grid(x_min[i], x_max[i], x_n, y_min[i], y_max[i], y_n, z_min[i], z_max[i], z_n)

            X, Y, Z = self.np.mgrid[x_min[i]:x_max[i]:x_n * 1j, y_min[i]:y_max[i]:y_n * 1j, z_min[i]:z_max[i]:z_n * 1j]
            AO_buf = self.np.zeros([self.nb, self.np.shape(X)[0], self.np.shape(Y)[1], self.np.shape(Z)[2]], dtype=self.np.float64)
            MO_buf = self.np.zeros( [2 * self.nGeminal + self.inactive, self.np.shape(X)[0], self.np.shape(Y)[1], self.np.shape(Z)[2]], dtype=self.np.float64)
            GEMINAL_buf = self.np.zeros([self.np.shape(X)[0], self.np.shape(Y)[1], self.np.shape(Z)[2]], dtype=self.np.float64)            
            
            self.Calc_AO_2(AO_buf, X, Y, Z)

            for j in range( 2 * self.nGeminal + self.inactive):

                for k in range(self.nb):
                    MO_buf[j, :, :, :] += self.Coeff[j, k] * AO_buf[k, :, :, :]

            for j in range(len(self.Orb2Gem)):
                if (i == self.Orb2Gem[j]):
                    GEMINAL_buf[:, :, :] += (self.G_coeff[j] * MO_buf[j + self.inactive, :, :, :]) ** 2

            self.Geminnals.append(Geminal(grid, GEMINAL_buf))

    def Plot_Geminal(self, Plot_Atoms=1, Plot_Atom_Names=1, Plot_Bonds=1, Geminal_number=0, Add_atom_number=False):
        #
        # from mayavi import mlab

        fig1 = self.mlab.figure("Geminal_" + str(Geminal_number), bgcolor=(.5, .5, .75), size=(1000, 1000))
        fig1.scene.parallel_projection = False

        self.mlab.clf()

        for i in Geminal_number:
            X, Y, Z = self.np.mgrid[self.Geminnals[i].g.x_min:self.Geminnals[i].g.x_max:self.Geminnals[i].g.x_n * 1j,
                      self.Geminnals[i].g.y_min:self.Geminnals[i].g.y_max:self.Geminnals[i].g.y_n * 1j,
                      self.Geminnals[i].g.z_min:self.Geminnals[i].g.z_max:self.Geminnals[i].g.z_n * 1j]
            #          exec(self.Geminnals[i].g)
            self.mlab.contour3d(X, Y, Z, (self.Geminnals[i].v), contours=12, opacity=0.5)

        if Plot_Atoms:
            for i in range(self.nAtoms):
                self.mlab.points3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2],
                                   scale_factor=self.Atoms_Scale[letters(self.Atoms_Name[i])], resolution=20,
                                   color=self.Atoms_Color[letters(self.Atoms_Name[i])], scale_mode='none')

        if Plot_Atom_Names:
            for i in range(self.nAtoms):
                if Add_atom_number:
                    self.mlab.text3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], self.Atoms_Name[i]+'_'+str(i), scale=(.75, .75, .75))
                else:
                    self.mlab.text3d(self.Atoms_R[i, 0], self.Atoms_R[i, 1], self.Atoms_R[i, 2], self.Atoms_Name[i], scale=(.75, .75, .75))


        if Plot_Bonds:
            for i in range(len(self.Bonds)):
                self.mlab.plot3d(
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 0], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 0]]),
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 1], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 1]]),
                    self.np.array([self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][0]), 2], self.Atoms_R[self.Atoms_Name.index(self.Bonds[i][1]), 2]]),
                    tube_radius=0.2, color=(233.0 / 255, 165.0 / 255, 165.0 / 255))

        self.mlab.show()

    def normalization_sumation(self, Data, pow_val):

        if (len(self.np.shape(Data)) == 2):

            NOrb = self.np.shape(Data)[-1] - 1
            NExpans = self.np.shape(Data)[0]
            Norm = self.np.zeros(NOrb, dtype=self.np.float64)

            for i in range(NExpans):
                for j in range(i + 1):
                    Norm += Data[i, 1:] * Data[j, 1:] / (Data[i, 0] + Data[j, 0]) ** pow_val

        else:

            Norm += Data[1] * Data[1] / (Data[0] + Data[0]) ** pow_val

        return Norm


    def Norm_S2(self, Data):

        pow_val = 3.0 / 2.0
        fact = 1.0

        Norm =  self.normalization_sumation( Data, pow_val) * self.np.pi ** (3.0 / 2.0) * fact
        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_P2(self, Data):

        pow_val = 5.0 / 2.0
        fact = 1.0 / 2.0

        Norm =  self.normalization_sumation( Data, pow_val) * self.np.pi ** (3.0 / 2.0) * fact
        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_D2(self, Data):

        pow_val = 7.0 / 2.0
        fact = 3.0 / 4.0

        Norm =  self.normalization_sumation( Data, pow_val) * self.np.pi ** (3.0 / 2.0) * fact

        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_F2(self, Data):

        pow_val = 9.0 / 2.0
        fact = 15.0 / 8.0

        Norm =  self.normalization_sumation( Data, pow_val) * self.np.pi ** (3.0 / 2.0) * fact
        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_G2(self, Data):

        pow_val = 11.0 / 2.0
        fact = 105.0 / 16.0

        Norm =  self.normalization_sumation( Data, pow_val) * self.np.pi ** (3.0 / 2.0) * fact       
        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_H2(self, Data):

        pow_val = 13.0 / 2.0
        fact = 945.0 / 32.0

        Norm =  self.normalization_sumation( Data, pow_val) * self.np.pi ** (3.0 / 2.0) * fact
        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Calc_AO_2(self, AO, X, Y, Z):

        sqrt3 = self.np.sqrt( 3.0 )
        sqrt3p4 = self.np.sqrt( 3.0/4.0) 
        sqrt5 = self.np.sqrt( 5.0 )
        sqrt15 = self.np.sqrt( 15.0 )
        sqrt10p16 = self.np.sqrt( 10.0 / 16.0 )
        sqrt45p100 = self.np.sqrt( 45.0 / 100 )
        sqrt40_25 = self.np.sqrt( 40.0 / 25 )
        sqrt125p108_1440 = self.np.sqrt( 125.0 / ( 108 - self.np.sqrt( 1440 ) ) )
        sqrt30p25 = self.np.sqrt( 30.0 / 25.0 )
        sqrt50p16 = self.np.sqrt( 50.0 / 16.0 )

        if self.Spherical:
            dx = X[1, 0, 0] - X[0, 0, 0]
            dy = Y[0, 1, 0] - Y[0, 0, 0]
            dz = Z[0, 0, 1] - Z[0, 0, 0]

            dv = dx * dy * dz

        nOrb = 0

        for n in range(self.nAtoms):

            x = X - self.Atoms_R[n, 0]
            y = Y - self.Atoms_R[n, 1]
            z = Z - self.Atoms_R[n, 2]

            RR = x * x + y * y + z * z

            r = self.np.sqrt(RR)
            theta = self.np.arctan(y / x)
            phi = self.np.arccos(z / r)

            for j in range(len(self.basis[n])):

                if (j == 0):

                    for k in range(self.np.shape(self.basis[n][j])[1] - 1):

                        for l in range(self.np.shape(self.basis[n][j])[0]):
                            # if (self.basis[n][j])[l, k+1]  > 10**-16:
                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR)  # / ( (self.basis[n][j])[l,0]  **( 3.0/2.0))

                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                        nOrb += 1

                if (j == 1):

                    for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                        for m in ['x', 'y', 'z']:
                            if (m == 'x'):

                                for l in range(self.np.shape(self.basis[n][j])[0]):
                                    AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x)
                                AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                nOrb += 1

                            if (m == 'y'):

                                for l in range(self.np.shape(self.basis[n][j])[0]):
                                    AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y)
                                AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                nOrb += 1

                            if (m == 'z'):

                                for l in range(self.np.shape(self.basis[n][j])[0]):
                                    AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (z)
                                AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                nOrb += 1
                if (j == 2):
                    if (self.Spherical == 1):

                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in range(-j, j + 1, 1):

                                    if (m == -2):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (sqrt3 * x * y)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == -1):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (sqrt3 * y * z)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 0):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (-0.5 * (x * x + y * y) + z * z)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 1):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (sqrt3 * x * z)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 2):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (
                                                                             sqrt3p4 * (x * x - y * y))
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1
                    else:
                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']:

                                    if (m == 'xx'):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x ** 2)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'yy'):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y ** 2)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'zz'):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (z ** 2)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xy'):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * y * sqrt3)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xz'):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -( self.basis[n][j])[l, 0] * RR) * (x * z * sqrt3)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'yz'):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y * z * sqrt3)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                if (j == 3):
                    if (self.Spherical == 1):

                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in range(-j, j + 1, 1):

                                    if (m == -3):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (sqrt10p16 * ( 3 * x ** 2 * y - y ** 3))
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == -2):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( sqrt15 * (x * y * z))
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == -1):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( sqrt125p108_1440 * sqrt45p100 * ( -x ** 2 * y - y ** 3 + sqrt40_25 * y * z ** 2))
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 0):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( z ** 3 - sqrt45p100 * ( sqrt5 * x ** 2 * z + sqrt5 * y ** 2 * z))
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 1):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( sqrt125p108_1440 * sqrt45p100 * ( -x ** 3 - x * y ** 2 + sqrt40_25 * x * z ** 2))
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 2):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( sqrt30p25 * sqrt50p16 * (x ** 2 * z - y ** 2 * z))
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 3):
                                        # print "nOrb: ", nOrb, " typ d_" +m
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( sqrt10p16 * ( x ** 3 - 3 * x * y ** 2))
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                    else:
                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in ['xxx', 'yyy', 'zzz', 'xyy', 'xxy', 'xxz', 'xzz', 'yzz', 'yyz', 'xyz']:

                                    if (m == 'xxx'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * x * x)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'yyy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y * y * y)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'zzz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (z * z * z)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xyy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * y * y * sqrt5)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xxy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * x * y * sqrt5)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xxz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * x * z * sqrt5)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xzz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * z * z * sqrt5)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'yzz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y * z * z * sqrt5)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'yyz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y * y * z * sqrt5)
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xyz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * y * z * sqrt15 )
                                        AO[nOrb, :, :, :] = (self.basis_norm2[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1
