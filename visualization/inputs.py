


#dddd = DaltonInput( )
#INPUT_TYPES = {'Dalton': DaltonInput() }

class Input( ):

    import tarfile
    import numpy as np

    input_type = None
    input_type = None
    input_name = None
    file_string = None
    BAS_filename = None
    data_source = None
    
    spherical = None
    nb = None
    nAtoms = None
    inactive = None
    electrons = None
    Occ = None

    Coeff= None

    Atoms_R = None
    Atoms_Charge = None
    Atoms_Name = None

    Bonds = None

    def __init__(self, input_type=None, input_sub_type=None, input_name=None, file_string=None, BAS_filename=None, data_source=None ):


        if input_type is not None:
            self.input_type = input_type

        if input_sub_type is not None:
            self.input_sub_type = input_sub_type

        if input_name is not None:
            self.input_name = input_name

        if file_string is not None:
            self.file_string = file_string

        if BAS_filename is not None:
            self.BAS_filename = BAS_filename

        if data_source is not None:
            self.data_source = data_source
    
#    @staticmethod
#    def get_input( input_type=None, input_sub_type=None, input_name=None, file_string=None, BAS_filename=None, data_source=None ):#
#
#        try:
#            #data_input = Input.input_types[input_type]
#            data_input = Input( input_type=input_type, input_sub_type=input_sub_type, input_name=input_name, file_string=file_string, BAS_filename=BAS_filename, data_source=data_source )
#            data_input.__class__ = Input.input_types[input_type]
#        except KeyError:
#            Exception(f"Input_type: {input_type} not found")

#        print(data_input)
#        return data_input

        #return type(data_input, (data_input,), {'input_type' : input_type, 'input_sub_type':input_sub_type, 'input_name':input_name, 
#'file_string':file_string, 'BAS_filename':BAS_filename, 'data_source':data_source})
        #return data_input( input_type = input_type, input_sub_type=input_sub_type, input_name=input_name, file_string=file_string, BAS_filename=BAS_filename, data_source=data_source)
        #if input_type in Input.input_types:
        #    pass

        #else:


    

    def set_input_type(self, input_type):
        pass

    def set_input_sub_type(self, input_sub_type):
        pass

    def set_input_name(self, source):
        pass

    def get_data_source(self, data_source):
        pass

    def get_harmonic(self):
        pass

    def get_nb(self):
        pass

    def get_nAtoms(self):
        pass

    def get_inactive(self):
        pass

    def get_electrons(self):
        pass

    def get_Occ(self):
        pass

    def get_Coeff(self):
        pass

    def get_Atoms(self):
        pass

    def get_Bonds(self):
        pass

class DaltonInput(Input):

    input_name = 'Dalton'

    dalton_output = None
    F_BAS = None

    def set_source(self, source):

        self.input_type = source

    def get_data_source(self, data_source):

        self.data_source = data_source

    def get_Dalton_Output(self):

        if self.dalton_output is not None:

            return self.dalton_output

        else:
            with open(self.input_name + ".out", 'r') as f:
                self.dalton_output = f.read()

            return self.dalton_output

    def get_spherical(self):

        Out = self.get_Dalton_Output()

        if (Out.find("Spherical harmonic basis used.") > 0):
            self.spherical = 1
        else:
            self.spherical = 0

        return self.spherical

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
            return self.Occ

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

        F_MOPUN = None

        if self.Coeff is not None:
            return self.Coeff

        if self.input_type == 'MOPUN':
            with open("DALTON.MOPUN", 'r') as f:
                F_MOPUN = f.read()

        if self.input_type == 'tar':
            tar = self.tarfile.open(self.input_name + ".tar.gz")
            f = tar.extractfile(tar.getmember("DALTON.MOPUN"))
            F_MOPUN = f.read().decode(encoding='utf-8')
            tar.close()

        F_MOPUN = F_MOPUN.replace('-', ' -')

        a = " ".join( F_MOPUN[F_MOPUN.find("\n"):].split() )
        b = self.np.fromstring( a, dtype=self.np.float64, sep=' ')

        self.Coeff = self.np.reshape(self.np.fromstring(" ".join(F_MOPUN[F_MOPUN.find("\n"):].split()),
                                                        dtype=self.np.float64, sep=' '), [self.nb, self.nb])

        return self.Coeff

    def get_Basis_File(self):

        if self.F_BAS is not None:
            return self.F_BAS

        if self.input_type == 'MOPUN':
            with open("DALTON.BAS", 'r') as f:
                self.F_BAS = f.read()

        if self.input_type == 'tar':
            tar = self.tarfile.open(self.input_name + ".tar.gz")
            f = tar.extractfile(tar.getmember("DALTON.BAS"))
            self.F_BAS = f.read().decode(encoding='utf-8')
            tar.close()

        return self.F_BAS

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
                self.Bonds.append(  [(Bonds_str.split("\n")[i]).split()[2], (Bonds_str.split("\n")[i]).split()[3],
                                    float((Bonds_str.split("\n")[i]).split()[4])])
        else:
            for i in range(len(Bonds_str.split("\n")) - 3):
                self.Bonds.append(  [(Bonds_str.split("\n")[i]).split()[2], (Bonds_str.split("\n")[i]).split()[4],
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
        Atom_Geometries_pattern = '^\w{1,6} \D \D \D'
        Basis_header_pattern = '^H {1,3}\d{1,3} {1,4}\d{1,3}$'

        Atoms_Gropus = []

        for line in F_BAS_split_lines[7:]:

            if re.match(Atom_header_pattern, line):

                Atoms_Group = {'headerLine': line, 'Geometries': [], 'Basis': []}
                Atoms_Gropus.append(Atoms_Group)

            elif re.match(Atom_Geometries_pattern, line):

                Atoms_Group['Geometries'].append(line)

            elif re.match(Basis_header_pattern, line):

                basis_part = {'header': line, 'data': []}
                Atoms_Group['Basis'].append(basis_part)

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

        for Atoms_Group in Atoms_Gropus:

            # TODO basis
            Orbitals = []

            for Orbital in Atoms_Group['Basis']:
                dim = self.np.fromstring(Orbital['header'][1:], dtype=self.np.int64, sep=' ')
                dim[1] += 1

                Orbital = self.np.reshape(self.np.fromstring(''.join(Orbital['data']), dtype=self.np.float64, sep=' '), dim)
                Orbitals.append(Orbital)

            Atom_Charge = int(Atoms_Group['headerLine'][:Atoms_Group['headerLine'].find('.')])

            for Atom in Atoms_Group['Geometries']:
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



    def normalization_summation(self, Data, pow_val):

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

        Norm =  self.normalization_summation( Data, pow_val) * self.np.pi ** (3.0 / 2.0) * fact
        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_P2(self, Data):

        pow_val = 5.0 / 2.0
        fact = 1.0 / 2.0

        Norm =  self.normalization_summation( Data, pow_val) * self.np.pi ** (3.0 / 2.0) * fact
        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_D2(self, Data):

        pow_val = 7.0 / 2.0
        fact = 3.0 / 4.0

        Norm =  self.normalization_summation( Data, pow_val) * self.np.pi ** (3.0 / 2.0) * fact

        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_F2(self, Data):

        pow_val = 9.0 / 2.0
        fact = 15.0 / 8.0

        Norm =  self.normalization_summation( Data, pow_val) * self.np.pi ** (3.0 / 2.0) * fact
        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_G2(self, Data):

        pow_val = 11.0 / 2.0
        fact = 105.0 / 16.0

        Norm =  self.normalization_summation( Data, pow_val) * self.np.pi ** (3.0 / 2.0) * fact       
        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_H2(self, Data):

        pow_val = 13.0 / 2.0
        fact = 945.0 / 32.0

        Norm =  self.normalization_summation( Data, pow_val) * self.np.pi ** (3.0 / 2.0) * fact
        Norm = 1 / self.np.sqrt(Norm)

        return Norm

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

INPUT_TYPES = {'Dalton': DaltonInput}

def get_input( input_type=None, input_sub_type=None, input_name=None, file_string=None, BAS_filename=None, data_source=None ):

    try:
        data_input = INPUT_TYPES[input_type](   input_type=input_sub_type, 
                                                input_name=input_name, 
                                                file_string=file_string, 
                                                BAS_filename=BAS_filename, 
                                                data_source=data_source )
    except KeyError:
        Exception(f"Input_type: {input_type} not found")

    return data_input

