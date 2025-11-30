import numpy as np

from visualization.inputs import Input


class MolproInput(Input):
	"""

		Input from molpro

	"""
	input_name = 'Molpro'

	output = None
	F_BAS = None

	def get_output(self):

		if self.output is not None:

			return self.output

		else:
			with open(self.input_name + ".out", 'r', encoding="utf-8") as f:
				_output = f.read()

			geometry_block_beginning_sentence = '1PROGRAM * SEWARD (Integral evaluation for generally contracted gaussian basis sets)     Author: Roland Lindh, 1990'

			self.output = _output.split(geometry_block_beginning_sentence)[1]

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

		_separator = ['']
		_separator_countdown = 3
		for i, atomic_data_line in enumerate( _atomic_data ):

			if atomic_data_line in _separator:
				_separator_countdown -= 1
				continue

			if not _separator_countdown:
				break
			elif _separator_countdown == 1:
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

	def get_occ(self):

		if self.Occ is not None:
			return self.Occ

		self.Occ = np.zeros([self.nb], dtype=np.float64)

		_output = self.get_output()

		beginning_orbital_part = 'Orb     Occ        Energy       Coefficients'
		ending_orbital_part = 'orbital dump at molpro section'
		empty_line = '\n\n'

		orbital_buf = _output[_output.find( beginning_orbital_part ):_output.find( ending_orbital_part )].split( empty_line )[2:-1]

		for i , orbital in enumerate(orbital_buf):
			orbital_split = orbital.split()
			#number = orbital_split[0]
			Occ = orbital_split[1]
			#energy = orbital_split[2]
			#coeff = orbital_split[2:]
			#print(coeff)

			self.Occ[i] = 0.5 *float(  Occ )

		#for i in range(self.electrons):
		#    self.Occ[self.inactive + i] = float(Geminal_buf[i].split()[6])

		#self.Occ[:self.inactive] = 1

		return self.Occ


	def get_coeff(self):

		self.Coeff = np.zeros([self.nb, self.nb], dtype=np.float64)

		_output = self.get_output()


		program_split_str = "1PROGRAM *"

		_output_last_program = _output.split(program_split_str)[-1]

		# output_programs[-1]

		beginning_orbital_data = 'Orb     Occ        Energy       Coefficients'
		ending_orbital_data = '*****************************************************************************************'
		empty_line = '\n\n'

		_orbital_data = _output_last_program[   _output_last_program.find( beginning_orbital_data   ):
												_output_last_program.find( ending_orbital_data      )].split( empty_line )[2:-3]


		for i , orbital in enumerate(_orbital_data):
			orbital_split = orbital.split()
			Occ = orbital_split[1]
			coeff = orbital_split[3:]
			coeff_line = np.fromstring(' '.join(coeff), dtype=np.float64, sep=' ')
			self.Coeff[i,:] = coeff_line


		return self.Coeff


	def get_atoms(self):

		if self.Atoms_R is not None and self.Atoms_Charge is not None and self.Atoms_Name is not None:
			return self.Atoms_R, self.Atoms_Charge, self.Atoms_Name

		self.Atoms_R = np.zeros([self.nAtoms, 3], dtype=np.float64)
		self.Atoms_Charge = np.zeros(self.nAtoms, dtype=np.int64)
		self.Atoms_Name = []

		_output = self.get_output()

		beginning_atomic_data = ' ATOMIC COORDINATES'
		ending_atomic_data = ' BASIS DATA'
		empty_line = '\n\n'

		_atomic_data = _output[_output.find( beginning_atomic_data ):_output.find( ending_atomic_data )].split( empty_line )[2].splitlines()

		for i , atom in enumerate( _atomic_data):
			atom_split = atom.split()
			self.Atoms_Name.append( atom_split[1] )
			self.Atoms_Charge[i] = self.Atoms_Charge[i] = int( float(atom_split[2]) )
			self.Atoms_R[i,:] = np.fromstring(' '.join(atom_split[3:]), dtype=np.float64, sep=' ')

		return self.Atoms_R, self.Atoms_Charge, self.Atoms_Name


	def get_bonds(self):

		_bonds = []

		beginning_atomic_data = ' ATOMIC COORDINATES'
		ending_atomic_data = ' BASIS DATA'
		empty_line = '\n\n'

		bond_length_start = 'Bond lengths'
		bond_angles_start = 'Bond angles'

		_output = self.get_output()

		_atomic_data = _output[_output.find( beginning_atomic_data ):_output.find( ending_atomic_data )]
		_bond_length_data_lines = _atomic_data[ _atomic_data.find( bond_length_start):_atomic_data.find( bond_angles_start) ].split(empty_line)[1:-1]

		for bonds_twoline in _bond_length_data_lines:
			_bonds_line = bonds_twoline.splitlines()[0].replace('-', ' ' ).split( )

			for i in range( int( len(_bonds_line)/3 ) ):
				_bond = [ int(_bonds_line[0 + i*3])-1, int(_bonds_line[1 + i*3])-1, float(_bonds_line[2 + i*3 ]) ]
				_bonds.append(_bond)

		self.Bonds = _bonds

		return self.Bonds

	def _read_basis_to_list(self, basis_data:str ):

		_basis:list = []

		_atom_basis:list = []
		_atom_number: int = 0

		_orbital_type_basis:list = []
		_orbitals_type:int = 0

		_orbital_basis:list = []

		_orbitals_coefficient_count:int = 0
		_orbitals_exponent_count:int = 0


		_header_length:int = 4

		for i, basis_line in enumerate ( basis_data ):

			beginning = basis_line[:21].split()
			data = basis_line[21:].split()

			if len(beginning) ==_header_length:

				if _atom_number != int( beginning[2] ):

					_orbitals_coefficient_count = 0
					_orbitals_exponent_count = 0

					_atom_basis = []
					_basis.append(_atom_basis)
					_atom_number = int( beginning[2] )

				if _orbitals_type != int( beginning[3][0]  ):

					_orbital_type_basis = []
					_atom_basis.append( _orbital_type_basis )
					_orbitals_type = int( beginning[3][0] )

				_orbital_basis = []
				_orbital_type_basis.append( _orbital_basis )

			if beginning:
				_orbitals_coefficient_count += 1

			_orbitals_exponent_count += 1

			_orbital_basis.append(data)

		return _basis


	def get_basis(self):
		""" get basis for the molpro input """

		_output = self.get_output()

		beginning_basis_data = ' BASIS DATA'
		ending_basis_data = ' NUCLEAR CHARGE:'
		empty_line = '\n\n'

		_basis_data = _output[_output.find( beginning_basis_data ):_output.find( ending_basis_data )].split( empty_line )[2].splitlines()

		_basis:list = []


		_basis_read = self._read_basis_to_list( basis_data = _basis_data )   # type: ignore

		_orbitals_count:int = 0
		_processed_orbitals:int = 0

		_orbital_exponent_set: type( set())    # type: ignore

		for _atom_basis in _basis_read:

			_atom_basis_process = []
			_basis.append(_atom_basis_process)


			for i, _orbital_type_basis in enumerate( _atom_basis ):

				_orbital_exponent_set = set([])
				_orbitals_count = 0

				for _orbital_group in  _orbital_type_basis[0::(2*i+1) ]:

					for j, _orbital_entry in enumerate(_orbital_group):

						if not j:
							_orbitals_count += len(_orbital_entry) -1

						_orbital_exponent_set = _orbital_exponent_set.union( [ float( _orbital_entry[0] )] )

				_orbital_type_basis_array = np.zeros( [len(_orbital_exponent_set) , _orbitals_count+1 ] ,dtype=np.float32)
				_atom_basis_process.append(_orbital_type_basis_array)

				_orbital_type_basis_array[:,0] = np.sort(np.array( list( _orbital_exponent_set ) ))[::-1]

				_processed_orbitals = 0
				for _orbital_group in  _orbital_type_basis[0::(2*i+1) ]:

					for j, _orbital_entry in enumerate(_orbital_group):

						_exponent_mask =_orbital_type_basis_array[:,0] == float( _orbital_entry[0] )
						_number_of_coefficients = len(_orbital_entry) - 1
						_exponent_line = np.array( [ float( coefficient ) for coefficient in _orbital_entry[1:] ]  ,dtype=np.float32)

						_orbital_type_basis_array[_exponent_mask, _processed_orbitals + 1:(_processed_orbitals + _number_of_coefficients + 1) ] = _exponent_line

						if j == len(_orbital_group)-1:
							_processed_orbitals += _number_of_coefficients


		basis_norm =  self.calc_norm_from_basis( basis = _basis)
		basis_norm2 = basis_norm

		self.basis = _basis
		self.basis_norm = basis_norm
		self.basis_norm2 = basis_norm2

		return self.basis, self.basis_norm, self.basis_norm2

class MolproSaptInput(MolproInput):
	"""

		Input from molpro

	"""
	input_name = 'Molpro'
	monomer:int = 0


	def __init__(self, *args, **kwargs):
		super(MolproSaptInput, self).__init__(*args, **kwargs)

		if 'monomer' in kwargs:
			self.monomer = kwargs['monomer']


	def get_output(self):

		if self.output is not None:

			return self.output

		else:
			with open(self.input_name + ".out", 'r', encoding="utf-8") as f:
				_output = f.read()

			geometry_block_beginning_sentence = 'Geometry written to block  1 of record 700'

			self.output = _output.split(geometry_block_beginning_sentence)[ 1+ self.monomer ]

			return self.output