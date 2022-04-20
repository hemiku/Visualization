
from visualization.inputs import Input


class MolproInput(Input):

    input_name = 'Molpro'

    output = None
    F_BAS = None

    def get_output(self):

        if self.output is not None:

            return self.output

        else:
            with open(self.input_name + ".out", 'r', encoding="utf-8") as f:
                _output = f.read()

            geometry_block_begining_sentence = 'Geometry written to block  1 of record 700' 

            self.output = _output.split(geometry_block_begining_sentence)[1]

            return self.output


    def get_spherical(self):

        _output = self.get_output()

        if (_output.find("Using spherical harmonics") > 0):
            self.spherical = 1
        else:
            self.spherical = 0

        return self.spherical

    def get_nb(self):

        if self.nb is not None:
            return self.nb

        _output = self.get_output()

        search_string = 'NUMBER OF CONTRACTIONS:'
        search_offset = 38

        self.nb = int(_output[_output.find(search_string):_output.find(search_string) + search_offset].split()[3] )

        return self.nb

    def get_nAtoms(self):

        if self.nAtoms is not None:
            return self.nAtoms

        _output = self.get_output()

        nAtoms = 0 
        _atomic_data = _output[_output.find("ATOMIC COORDINATES"):_output.find("BASIS DATA") + 30].splitlines()

        _sepparator = ['']
        _sepparator_countdown = 3
        for i, atomic_data_line in enumerate( _atomic_data ): 

            if atomic_data_line in _sepparator:
                _sepparator_countdown -= 1       
                continue    
            
            if not _sepparator_countdown:
                break
            elif _sepparator_countdown == 1:
                nAtoms += 1

        self.nAtoms = nAtoms


        return self.nAtoms

    def get_inactive(self):

        if self.inactive is not None:
            return self.inactive

        _output = self.get_output()
        
        search_string = 'Number of closed-shell orbitals:'
        search_offset = 38

        self.inactive = int(_output[_output.find(search_string):_output.find(search_string) + search_offset].split()[2])

        return self.inactive

    def get_electrons(self):

        if self.electrons is not None:
            return self.electrons

        _output = self.get_output()

        search_string = 'NUCLEAR CHARGE:'
        search_offset = 38

        self.electrons = int(_output[_output.find(search_string):_output.find(search_string) + search_offset].split()[2])

        return self.electrons

    def get_Occ(self):

        if self.Occ is not None:
            return self.Occ

        self.Occ = self.np.zeros([self.nb], dtype=self.np.float64)

        _output = self.get_output()

        beginning_orbital_part = 'Orb     Occ        Energy       Coefficients'
        ending_orbital_part = 'orbital dump at molpro section'
        empty_line = '\n\n'

        orbital_buf = _output[_output.find( beginning_orbital_part ):_output.find( ending_orbital_part )].split( empty_line )[2:-1]

        for i , orbital in enumerate(orbital_buf):
            orbital_splti = orbital.split()
            #number = orbital_splti[0]
            Occ = orbital_splti[1]
            #energy = orbital_splti[2]
            #coeff = orbital_splti[2:]
            #print(coeff)

            self.Occ[i] = 0.5 *float(  Occ ) 

        #for i in range(self.electrons):
        #    self.Occ[self.inactive + i] = float(Geminal_buf[i].split()[6])

        #self.Occ[:self.inactive] = 1

        return self.Occ


    def get_Coeff(self):

        self.Coeff = self.np.zeros([self.nb, self.nb], dtype=self.np.float64)

        _output = self.get_output()

        beginning_orbital_data = 'Orb     Occ        Energy       Coefficients'
        ending_orbital_data = 'orbital dump at molpro section'
        empty_line = '\n\n'

        _orbital_data = _output[_output.find( beginning_orbital_data ):_output.find( ending_orbital_data )].split( empty_line )[2:-1]

        for i , orbital in enumerate(_orbital_data):
            orbital_splti = orbital.split()
            #number = orbital_splti[0]
            Occ = orbital_splti[1]
            #energy = orbital_splti[2]
            coeff = orbital_splti[3:]
            #print(coeff)
            coeff_line = self.np.fromstring(' '.join(coeff), dtype=self.np.float64, sep=' ')
            self.Coeff[i,:] = coeff_line

        #F_MOPUN = None

        #if self.Coeff is not None:
        #    return self.Coeff

        #if self.input_type == 'MOPUN':
        #    with open("DALTON.MOPUN", 'r') as f:
        #        F_MOPUN = f.read()

        #if self.input_type == 'tar':
        #    tar = self.tarfile.open(self.input_name + ".tar.gz")
        #    f = tar.extractfile(tar.getmember("DALTON.MOPUN"))
        #    F_MOPUN = f.read().decode(encoding='utf-8')
        #    tar.close()

        #F_MOPUN = F_MOPUN.replace('-', ' -')

        #a = " ".join( F_MOPUN[F_MOPUN.find("\n"):].split() )
        #b = self.np.fromstring( a, dtype=self.np.float64, sep=' ')

        #self.Coeff = self.np.reshape(self.np.fromstring(" ".join(F_MOPUN[F_MOPUN.find("\n"):].split()),
        #                                                dtype=self.np.float64, sep=' '), [self.nb, self.nb])

        return self.Coeff


    def get_Atoms(self):

        if self.Atoms_R is not None and self.Atoms_Charge is not None and self.Atoms_Name is not None:
            return self.Atoms_R, self.Atoms_Charge, self.Atoms_Name

        self.Atoms_R = self.np.zeros([self.nAtoms, 3], dtype=self.np.float64)
        self.Atoms_Charge = self.np.zeros(self.nAtoms, dtype=self.np.int64)
        self.Atoms_Name = []

        _output = self.get_output()

        beginning_atomic_data = ' ATOMIC COORDINATES'
        ending_atomic_data = ' BASIS DATA'
        empty_line = '\n\n'

        _atomic_data = _output[_output.find( beginning_atomic_data ):_output.find( ending_atomic_data )].split( empty_line )[2].splitlines()

        for i , atom in enumerate( _atomic_data):
            atom_splti = atom.split()
            self.Atoms_Name.append( atom_splti[1] )
            self.Atoms_Charge[i] = self.Atoms_Charge[i] = int( float(atom_splti[2]) )
            self.Atoms_R[i,:] = self.np.fromstring(' '.join(atom_splti[3:]), dtype=self.np.float64, sep=' ')

        return self.Atoms_R, self.Atoms_Charge, self.Atoms_Name

    
    def get_Bonds(self):

        return None

    def get_Basis(self):

        _output = self.get_output()

        beginning_basis_data = ' BASIS DATA'
        ending_basis_data = ' NUCLEAR CHARGE:'
        empty_line = '\n\n'

        _basis_data = _output[_output.find( beginning_basis_data ):_output.find( ending_basis_data )].split( empty_line )[2].splitlines()




        return self.basis, self.basis_norm, self.basis_norm2
