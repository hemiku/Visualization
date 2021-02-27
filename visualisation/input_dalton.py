
import visualisation.inputs as input

class DaltonInput(input.Input):

    Dalton_Output = None

    def set_source(self, source):

        self.input_type = source

    def get_data_source(self, data_source):

        self.data_source = data_source

    def get_Dalton_Output(self):

        if self.Dalton_Output is not None:

            return self.Dalton_Output

        else:
            with open(self.input_name + ".out", 'r') as f:
                self.Dalton_Output = f.read()

            return self.Dalton_Output

    def get_harmonic(self):

        Out = self.get_Dalton_Output()

        if (Out.find("Spherical harmonic basis used.") > 0):
            self.Spherical = 1
        else:
            self.Spherical = 0

    def get_nb(self):

        if self.nb is not None:
            return self.nb

        Out = self.get_Dalton_Output()

        self.nb = int(Out[Out.find("Number of basis functions"):Out.find("Number of basis functions") + 38].split()[-1])

        return self.nb

    def get_nAtoms(self):

        if self.nAtoms is not None:
            return self.nAtoms

        Out = self.get_Dalton_Output()

        self.nAtoms = int(Out[Out.find("Total number of atoms:"):Out.find("Total number of atoms:") + 30].split()[-1])

        return self.nAtoms

    def get_inactive(self):

        if self.inactive is not None:
            return self.inactive

        Out = self.get_Dalton_Output()

        self.inactive = int(Out[Out.find(".INACTIVE"):Out.find("Number of basis functions") + 100].split()[1])

        return self.inactive

    def get_electrons(self):

        if self.electrons is not None:
            return self.electrons

        Out = self.get_Dalton_Output()

        self.electrons = int(Out[Out.find("@    Number of electrons in active shells"):Out.find("@    Total charge of the molecule") + 100].split()[7])

        return self.electrons

    def get_Occ(self):

        if self.Occ is not None:
            return self.self.Occ

        self.Occ = self.np.zeros([self.nb], dtype=self.np.float64)

        Out = self.get_Dalton_Output()

        Geminal_buf = Out[Out.find("APSG geminal coefficients and natural orbital occupations:"):Out.rfind(
            'NEWORB " orbitals punched.')].split(
            "====================================================================")[1].split("\n")[1:-1]

        for i in range(self.electrons):
            self.Occ[self.inactive + i] = float(Geminal_buf[i].split()[6])

        self.Occ[:self.inactive] = 1

        return self.Occ

    def get_Coeff(self):

        source_file_name = 'DALTON.MOPUN'

        if self.Coeff is not None:
            return self.Coeff

        if self.input_type == 'MOPUN':

            if self.input_name.count( '/' ):
                with open('/'.join(self.input_namesplit(sep='/')[:-1]) + '/'+ source_file_name, 'r') as f:
                    F_MOPUN = f.read()
            else:
                with open(source_file_name, 'r') as f:
                    F_MOPUN = f.read()

        if self.input_type == 'tar':
            tar = self.tarfile.open(self.input_name + ".tar.gz")
            f = tar.extractfile(tar.getmember(source_file_name))
            F_MOPUN = f.read().decode(encoding='utf-8')
            tar.close()

        F_MOPUN = F_MOPUN.replace('-', ' -')

        self.Coeff = self.np.reshape(self.np.fromstring(" ".join(F_MOPUN[F_MOPUN.find("\n"):].split()),
                                                        dtype=self.np.float64, sep=' '), [self.nb, self.nb])

        return self.Coeff

    def get_Basis_File(self):

        source_file_name = 'DALTON.BAS'

        if self.input_type == 'MOPUN':

            with open(source_file_name, 'r') as f:
                F_BAS = f.read()

            if self.input_name.count( '/' ):
                with open('/'.join(self.input_namesplit(sep='/')[:-1]) + '/' + source_file_name, 'r') as f:
                    F_BAS = f.read()
            else:
                with open( source_file_name , 'r') as f:
                    F_BAS = f.read()

        if self.input_type == 'tar':
            tar = self.tarfile.open(self.input_name + ".tar.gz")
            f = tar.extractfile(tar.getmember(source_file_name ))
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

        return self.Atoms_R, self.Atoms_Charge, self.Atoms_Name

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

        return self.Bonds

    def get_Basis(self):

        self.basis = []
        self.basis_norm = []
        self.basis_norm2 = []

        F_BAS = self.get_Basis_File()

        import re

        F_BAS_split_lines = F_BAS.splitlines()

        Atom_header_pattern = ' {1,9}\d{1,9}\. '
        Atom_Geometrie_pattern = '^\w{1,6} \D \D \D'
        Basis_header_pattern = '^H  {1,3}\d {1,4}\d$'

        Atoms_Gropus = []

        for line in F_BAS_split_lines[7:]:

            if re.match(Atom_header_pattern, line):

                Atoms_Gropu = {'headerLine': line, 'Geometries': [], 'Basis': []}
                Atoms_Gropus.append(Atoms_Gropu)

            elif re.match(Atom_Geometrie_pattern, line):

                Atoms_Gropu['Geometries'].append(line)

            elif re.match(Basis_header_pattern, line):

                basis_part = {'header': line, 'data': []}
                Atoms_Gropu['Basis'].append(basis_part)

            else:
                basis_part['data'].append(line)

        self.map_atoms_from_BAS_file(Atoms_Gropus)


        return self.basis, self.basis_norm, self.basis_norm2

    def map_atoms_from_BAS_file(self, Atoms_Gropus):
        """
        docstring
        """

        self.Atoms_R = self.np.zeros([self.nAtoms, 3], dtype=self.np.float64)
        self.Atoms_Charge = self.np.zeros(self.nAtoms, dtype=self.np.int64)
        self.Atoms_Name = []

        i = 0

        for Atoms_Gropu in Atoms_Gropus:

            # TODO basis
            Orbitals = []

            for Orbital in Atoms_Gropu['Basis']:
                dim = self.np.fromstring(Orbital['header'][1:], dtype=self.np.int64, sep=' ')
                dim[1] += 1

                Orbital = self.np.reshape(self.np.fromstring(''.join(Orbital['data']), dtype=self.np.float64, sep=' '),
                                          dim)
                Orbitals.append(Orbital)

            Atom_Charge = int(Atoms_Gropu['headerLine'][:Atoms_Gropu['headerLine'].find('.')])

            for Atom in Atoms_Gropu['Geometries']:
                Atom_split = Atom.split()

                Atom_name = Atom_split[0]
                Atom_R = self.np.array([float(Atom_split[1]), float(Atom_split[2]), float(Atom_split[3])])

                self.Atoms_Name.append(Atom_name)
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

    def get_G_coeff(self):

        self.G_coeff = self.np.zeros([self.electrons], dtype=self.np.float64)

        Out = self.get_Dalton_Output()
        Geminal_buf = Out[Out.find("APSG geminal coefficients and natural orbital occupations:"):
                          Out.rfind('NEWORB " orbitals punched.')].split(
            "====================================================================")[1].split("\n")[1:-1]

        for i in range(self.electrons):
            self.G_coeff[i] = float(Geminal_buf[i].split()[5])

        return self.G_coeff

    def get_Orb2Gem(self):

        self.Orb2Gem = self.np.zeros([self.electrons], dtype=self.np.int64)

        Out = self.get_Dalton_Output()
        Geminal_buf = Out[Out.find("APSG geminal coefficients and natural orbital occupations:"):
                          Out.rfind('NEWORB " orbitals punched.')].split(
            "====================================================================")[1].split("\n")[1:-1]

        for i in range(self.electrons):
            self.Orb2Gem[i] = int(Geminal_buf[i].split()[3]) - 1

        self.nGeminal = int(len(self.Orb2Gem) / 2)

        return self.Orb2Gem, self.nGeminal

