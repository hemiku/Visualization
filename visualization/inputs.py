


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

    basis = None
    basis_norm = None
    basis_norm2 = None

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
            with open(self.input_name + ".out", 'r', encoding="utf-8") as f:
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

        self.inactive = int(Out[Out.find("@    Inactive orbitals"):Out.find("@    Inactive orbitals") + 100].split()[3] )

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

    def find_header_end_line(self, lines , Atom_header_pattern = ' {1,9}\d{1,9}\. ' ):

        import re

        for i, line in enumerate(lines):
            if re.match(Atom_header_pattern, line):
                return i

    def get_Atoms(self):

        import re

        self.Atomsa_R = self.np.zeros([self.nAtoms, 3], dtype=self.np.float64)
        self.Atoms_Charge = self.np.zeros(self.nAtoms, dtype=self.np.int64)
        self.Atoms_Name = []

        F_BAS = self.get_Basis_File()
        F_BAS_split_lines = F_BAS.splitlines()

        Atom_header_pattern = ' {1,9}\d{1,9}\. '
        Atom_Geometries_pattern = '^\w{1,6} \D \D \D'
        Basis_header_pattern = '^H {1,3}\d{1,3} {1,4}\d{1,3}$'

        Atoms_Gropus = []

        header_end_line = self.find_header_end_line( lines = F_BAS_split_lines, Atom_header_pattern = Atom_header_pattern )

        for line in F_BAS_split_lines[header_end_line:]:

            if re.match(Atom_header_pattern, line):

                Atoms_Group = {'headerLine': line, 'Geometries': [], 'Basis': []}
                Atoms_Gropus.append(Atoms_Group)

            elif re.match(Atom_Geometries_pattern, line):

                Atoms_Group['Geometries'].append(line)

            elif re.match(Basis_header_pattern, line):

                basis_part = {'header': line, 'data': []}
                Atoms_Group['Basis'].append(basis_part)

            else:
                pass
                #basis_part['data'].append(line)

        self.map_atoms_geometry_from_BAS_file( Atoms_Gropus )

        return self.Atoms_R, self.Atoms_Charge, self.Atoms_Name

    def get_Bonds(self):

        self.Bonds = []

        Out = self.get_Dalton_Output()

        start_bond_section = "Bond distances "
        start_next_section = "| Starting in Integral Section (HERMIT) |"

        Bonds_str = Out[Out.find( start_bond_section ):Out.find(start_next_section)]
        Bonds_str = Bonds_str[Bonds_str.find("bond"):]
        Bonds_str_lines = Bonds_str.splitlines()

        for bond_str in Bonds_str_lines:
            bond_str_split = bond_str.split()

            try:
                self.Bonds.append(  [ bond_str_split[2], bond_str_split[3], float(bond_str_split[4]) ] )
                continue
            except:
                pass

            try:
                self.Bonds.append( [  bond_str_split[2], bond_str_split[4], float(bond_str_split[6])  ] )
                continue
            except:
                pass

            if bond_str == "" :
                break

        return self.Bonds

    def get_Basis(self):

        import re

        self.basis = []
        self.basis_norm = []
        self.basis_norm2 = []

        F_BAS = self.get_Basis_File()
        F_BAS_split_lines = F_BAS.splitlines()

        Atom_header_pattern = ' {1,9}\d{1,9}\. '
        Atom_Geometries_pattern = '^\w{1,6} \D \D \D'
        Basis_header_pattern = '^H {1,3}\d{1,3} {1,4}\d{1,3}$'

        Atoms_Gropus = []

        header_end_line = self.find_header_end_line( lines = F_BAS_split_lines, Atom_header_pattern = Atom_header_pattern )

        for line in F_BAS_split_lines[header_end_line:]:

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


    def map_atoms_geometry_from_BAS_file(self, Atoms_Gropus):
        """
        docstring
        """

        self.Atoms_R = self.np.zeros([self.nAtoms, 3], dtype=self.np.float64)
        self.Atoms_Charge = self.np.zeros(self.nAtoms, dtype=self.np.int64)
        self.Atoms_Name = []

        i = 0

        for Atoms_Group in Atoms_Gropus:

            Atom_Charge = int(Atoms_Group['headerLine'][:Atoms_Group['headerLine'].find('.')])

            for Atom in Atoms_Group['Geometries']:
                Atom_split = Atom.split()

                Atom_name = Atom_split[0]
                Atom_R = self.np.array([float(Atom_split[1]), float(Atom_split[2]), float(Atom_split[3])])

                self.Atoms_Name.append(Atom_name)
                self.Atoms_Charge[i] = Atom_Charge
                self.Atoms_R[i] = Atom_R

                i += 1

    def map_atoms_from_BAS_file(self, Atoms_Gropus):
        """
        docstring
        """

        self.Atoms_R = self.np.zeros([self.nAtoms, 3], dtype=self.np.float64)
        self.Atoms_Charge = self.np.zeros(self.nAtoms, dtype=self.np.int64)
        self.Atoms_Name = []

        self.basis = []
        self.basis_norm = []
        self.basis_norm2 = []

        i = 0

        for Atoms_Group in Atoms_Gropus:

            # TODO basis
            Orbitals = []

            for Orbital in Atoms_Group['Basis']:
                dim = self.np.fromstring(Orbital['header'][1:], dtype=self.np.int64, sep=' ')
                dim[1] += 1

                if len( Orbital['data'] ):
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

                if len( Orbitals ):
                    self.basis.append(Orbitals)

                i += 1

        if( len( self.basis )):
            for n in range(self.nAtoms):
                
                self.basis_norm.append([])
                
                for j in range(len(self.basis[n])):
                    if (j == 0):
                        self.basis_norm[n].append( self.Norm_S2(self.basis[n][j]) )
                    if (j == 1):
                        self.basis_norm[n].append(self.Norm_P2(self.basis[n][j]))
                    if (j == 2):
                        self.basis_norm[n].append(self.Norm_D2(self.basis[n][j]))
                    if (j == 3):
                        self.basis_norm[n].append(self.Norm_F2(self.basis[n][j]))
                    if (j == 4):
                        self.basis_norm[n].append(self.Norm_G2(self.basis[n][j]))

    def get_Printout_of_final_geminals(self, dalton_output ):

        GEMINAL_PART_START = "Printout of final geminals"
        WAWE_FUNCTION_SECTION_END = "| End of Wave Function Section (SIRIUS) |"

        _geminal_part = dalton_output[  dalton_output.find( GEMINAL_PART_START ) : 
                                        dalton_output.find( WAWE_FUNCTION_SECTION_END )  ]

        return _geminal_part

    def get_G_coeff(self):

        Coeff_separator = "===================================================================="
        Out = self.get_Dalton_Output()
        Geminal_buf = self.get_Printout_of_final_geminals( dalton_output = Out )
        coeff_buf = Geminal_buf.split( Coeff_separator )[1].split("\n")[1:-1]

        G_coeff = []

        for coeff in coeff_buf:
            G_coeff.append( float(coeff.split()[5]) )

        self.G_coeff = self.np.array( G_coeff, dtype=self.np.float32)

        return self.G_coeff

    def get_Orb2Gem(self):

        Coeff_separator = "===================================================================="
        Out = self.get_Dalton_Output()
        Geminal_buf = self.get_Printout_of_final_geminals( dalton_output = Out )
        coeff_buf = Geminal_buf.split( Coeff_separator )[1].split("\n")[1:-1]

        Orb2Gem = []

        for coeff in coeff_buf:
            Orb2Gem.append( int(coeff.split()[3]) - 1)

        
        self.Orb2Gem = self.np.array( Orb2Gem, dtype=self.np.int32)
        self.nGeminal = int( self.np.max(self.Orb2Gem) ) + 1 

        # self.Orb2Gem = self.np.zeros([self.electrons], dtype=self.np.int64)

        # Out = self.get_Dalton_Output()
        # Geminal_buf = Out[Out.find("APSG geminal coefficients and natural orbital occupations:"):
        #                   Out.rfind('NEWORB " orbitals punched.')].split(
        #     "====================================================================")[1].split("\n")[1:-1]

        # for i in range(self.electrons):
        #     self.Orb2Gem[i] = int(Geminal_buf[i].split()[3]) - 1

        # self.nGeminal = int(len(self.Orb2Gem) / 2)

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
                for j in range(NExpans):
                    Norm += Data[i, 1:] * Data[j, 1:] / (Data[i, 0] + Data[j, 0]) ** pow_val

        else:

            Norm = Data[1] * Data[1] / (Data[0] + Data[0]) ** pow_val

        return Norm


    def normalization_summation_2(self, Data, pow_val):

        if (len(self.np.shape(Data)) == 2):

            NOrb = self.np.shape(Data)[-1] - 1
            NExpans = self.np.shape(Data)[0]
            Norm = self.np.zeros(NOrb, dtype=self.np.float64)

            exponents = Data[:,0]
            coefficents =  Data[:,1:]

            for i in range(NExpans):
                for j in range(NExpans):
                    Norm += Data[i, 1:] *Data[j, 1:]  / (Data[i, 0] + Data[j, 0]) ** pow_val

        else:

            Norm + Data[1] * Data[1] / (Data[0] + Data[0]) ** pow_val

        return Norm



    def Norm_S2(self, Data):

        pow_val = 3.0 / 2.0
        fact = self.np.sqrt(2)/4
        fact = 1.0 
        
        Norm = fact * self.np.pi ** (3.0 / 2.0) * self.normalization_summation( Data, pow_val) 
        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_P2(self, Data):

        pow_val = 5.0 / 2.0
        fact = 1.0 / 2.0

        Norm = fact * self.np.pi ** (3.0 / 2.0) * self.normalization_summation( Data, pow_val) 
        Norm = 1 / self.np.sqrt(Norm)

        return Norm

    def Norm_D2(self, Data):

        pow_val = 7.0 / 2.0
        fact = 1.0 / 4.0

        Norm = fact * self.np.pi ** (3.0 / 2.0) * self.normalization_summation( Data, pow_val)

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

class MoldenInput(Input):

    input_name = 'molden'

    output = None

    def get_Output(self, Filename):
        file_type = '.molden'
        if Filename[-len(file_type):] == file_type:
            self.File_Name = Filename[:-len(file_type)]
        else:
            self.File_Name = Filename

    def get_file(self):

        if self.output is not None:

            return self.output

        else:
            with open(self.input_name, 'r', encoding="utf-8") as f:
                self.output = f.read()

            return self.output

    def get_data(self):

        Atoms_begin = 'Atoms]'
        GTO_begin   = 'GTO]'
        MO_begin   = 'MO]'

        Out = self.get_file()

        if Out.find('[Molden Format]') >= 0:
            Out_sections = Out.split('[')

            sections = {}

            for section  in Out_sections:

                section_key = section[:section.find(']')]
                section_data = section[section.find(']')+1:]
                
                sections[section_key] = section_data
            
            self.get_atoms(sections['Atoms'])
            self.get_GTOs(sections['GTO'])
            self.get_MOs(sections['MO'])

        self.Spherical = True
        self.calc_nb()

    def calc_nb(self):
        self.nb = len(self.MOs)

    def get_MOs(self, section):

        MO_sym = 'Sym'
        MO_Ene = 'Ene'
        MO_Spin = 'Spin'
        MO_Occup = 'Occup'

        self.MOs_list = []
        section_lines = section.splitlines()
        for j, line in enumerate(section_lines[1:]):
            if line.count(MO_Occup):
                #MO = []
                #self.MOs_list.append(MO)
                MO_dict = {}
            
            elif ( line.count(MO_sym) + line.count(MO_Ene) + line.count(MO_Spin) + line.count(MO_Occup) ) == 0:
                #MO.append( float( line.split()[1] ) )

                MO_dict[ int( line.split()[0] )  ] = float( line.split()[1] ) 

        self.MOs = self.np.zeros([len(self.MOs_list),len(self.MOs_list)], dtype=self.np.float64)

        for i, MO_dict in enumerate(self.MOs_list):
            for j in MO_dict:
                self.MOs[i,j] = MO_dict[j]
        
        #self.MOs = self.np.asarray( self.MOs_list )

    def get_GTOs(self, section):
        self.GTOs = []
        section_lines = section.splitlines()
        for j, line in enumerate(section_lines[1:]):
            line_split = line.split()
            if len(line_split) == 2 and len(line) < 18:
                atom_GTOs = [[], [], [], [], [], []]
                self.GTOs.append(atom_GTOs)
            elif len(line_split) == 3 and  line_split[0] == 's':
                GTO = []
                harmonic_GTOs = atom_GTOs[0]
                harmonic_GTOs.append(GTO)
#            elif line[:2] == ' p':
            elif len(line_split) == 3 and line_split[0] == 'p':
                GTO = []
                harmonic_GTOs = atom_GTOs[1]
                harmonic_GTOs.append(GTO)
            #elif line[:2] == ' d':
            elif len(line_split) == 3 and line_split[0] == 'd':
                GTO = []
                harmonic_GTOs = atom_GTOs[2]
                harmonic_GTOs.append(GTO)
            #elif line[:2] == ' f':
            elif len(line_split) == 3 and line_split[0] == 'f':
                GTO = []
                harmonic_GTOs = atom_GTOs[3]
                harmonic_GTOs.append(GTO)
#            elif line[:2] == '  ' and len(line) > 2:
            elif len(line_split) == 2 and line[18] == ' ':

                GTO.append( [float(line_split[0].replace('D', 'E')),   float(line_split[1].replace('D', 'E'))]  )

        self.basis = []

        for atom_GTOs in self.GTOs:
            atom_basis = []
            self.basis.append( atom_basis)
            for harmonic_GTOs in atom_GTOs:
                if len(harmonic_GTOs) > 0 :
                    atom_basis.append( self.transform_molden_basis( harmonic_GTOs ) )

    def transform_molden_basis(self, harmonic_GTOs ):

        Sigmas_list = []
        number_of_orbital = len(harmonic_GTOs)

        for GTO in harmonic_GTOs:
            for GTO_line in GTO:
                Sigmas_list.append(GTO_line[0])

        Sigmas = self.np.array(sorted(set(Sigmas_list), reverse=True))

        basis = self.np.zeros([Sigmas.shape[0], number_of_orbital + 1], dtype=self.np.float64)
        basis[:, 0] = Sigmas
        for i, GTO in enumerate( harmonic_GTOs):
            for GTO_line in GTO:
                mask = Sigmas == GTO_line[0]
                basis[mask, i + 1] = GTO_line[1]

        return basis

    def get_atoms(self, section):
        self.Atoms = []
        section_lines = section.splitlines()
        for j, line in enumerate(section_lines[1:]):
            self.Atoms.append(line)

        self.nAtoms = len(self.Atoms)

        self.Atoms_R = self.np.zeros([self.nAtoms,3],dtype=self.np.float64)
        self.Atoms_Charge = self.np.zeros(self.nAtoms,dtype=self.np.int64)
        self.Atoms_Name=[]

        for i, atom in enumerate(self.Atoms):

            atom_split = atom.split()
            self.Atoms_Name.append( atom_split[0] )
            self.Atoms_Charge[i] = int(atom_split[2])
            self.Atoms_R[i] = self.np.array([ float(atom_split[-3]),   float(atom_split[-2]) , float(atom_split[-1])   ] )

    def get_nb(self):
        
        if self.nb is None:
            self.get_data()
        return self.nb

    def get_AO_number(self):
        return self.MOs.shape[1]

    def get_Atoms_R(self):
        return self.Atoms_R

    def get_Atoms_Name(self):
        return self.Atoms_Name

    def get_Atoms_Charge(self):
        return self.Atoms_Charge

    def get_Basis(self):
        #return self.Basis
        return self.basis, self.basis_norm, self.basis_norm2

    def get_spherical(self):
        return self.spherical

    def get_nAtoms(self):
        return self.nAtoms

    def get_Coeff(self):
        return self.Coeff

    def get_Atoms(self):

        self.get_data()
        return self.Atoms_R, self.Atoms_Charge, self.Atoms_Name


INPUT_TYPES = { 'Dalton': DaltonInput,
                'molden': MoldenInput}

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
