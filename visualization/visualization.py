from visualization.geminals import GeminalGenerator

class VisualizationData( ):
    
    background_colors = {
        'Gaussian': (.5, .5, .75),
        'White': (1., 1., 1.)
    }

    Atoms_Color = {
            "H": (1, 1, 1),
            "He": (0.5, 0.5, 0.5),
            "Ar": (0.75, 0.75, 0.75),
            "O": (1, 0, 0),
            "B": (233.0 / 255, 165.0 / 255, 165.0 / 255),
            "N": (0, 0, 1),
            "C": (0, 0, 0),
            "F": (0, 1, 0),
            "Ne": (1.0 * 178 / 256, 1.0 * 229 / 256, 1.0 * 247 / 256)
        }

    Atoms_Scale = {
            "H": 1.0,
            "He": 1.0,
            "O": 1.4,
            "B": 1.4,
            "N": 1.4,
            "C": 1.4,
            "F": 1.4,
            "Ne": 1.4,
            "Ar": 1.4
        }

    Atoms_Mass = {
            "H": 1.00794, "He": 4.002602, "Li": 6.941, "Be": 9.012182, "B": 10.811, "C": 12.011, "N": 14.00674, "O": 15.9994, "F": 18.9984032, "Ne": 20.1797, 
            "Na": 22.989768, "Mg": 24.305, "Al": 26.981539, "Si": 28.0855, "P": 30.973762, "S": 32.066, "Cl": 35.4527, "Ar": 39.948, 
            "K": 39.0983, "Ca": 40.078, "Sc": 44.95591, "Ti": 47.88, "V": 50.9415, "Cr": 51.9961, "Mn": 54.93805, "Fe": 55.847, "Co": 58.9332, "Ni": 58.69, 
            "Cu": 63.546, "Zn": 65.39, "Ga": 69.723, "Ge": 72.61, "As": 74.92159, "Se": 78.96, "Br": 79.904, "Kr": 83.8, "Rb": 85.4678, "Sr": 87.62, "Y": 88.90585, "Zr": 91.224, "Nb": 92.90638, "Mo": 95.94,
        }
        
    Name_to_Atom_Number = {
            'H':1, 'He':2, 'Li':3, 'Be':4, 'B':5, 'C':6, 'N':7, 'O':8, 'F':9, 'Ne':10, 'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P':15, 'S':16, 'Cl':17, 'Ar':18, 
            'K':19, 'Ca':20, 'Sc':21, 'Ti':22, 'V':23, 'Cr':24, 'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 'Ga':31, 'Ge':32, 'As':33, 'Se':34, 'Br':35, 
            'Kr':36, 'Rb':37, 'Sr':38, 'Y':39, 'Zr':40, 'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48, 'In':49, 'Sn':50, 'Sb':51, 'Te':52, 
            'I':53, 'Xe':54, 'Cs':55, 'Ba':56, 'La':57, 'Ce':58, 'Pr':59, 'Nd':60, 'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64, 'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 
            'Yb':70, 'Lu':71, 'Hf':72, 'Ta':73, 'W':74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80, 'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 
            'Fr':87, 'Ra':88, 'Ac':89, 'Th':90, 'Pa':91, 'U':92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96, 'Bk':97, 'Cf':98, 'Es':99, 'Fm':100, 'Md':101, 'No':102, 'Lr':103, 
            'Rf':104, 'Db':105, 'Sg':106, 'Bh':107, 'Hs':108, 'Mt':109, 'Ds':110, 'Rg':111, 'Cn':112, 'Nh':113, 'Fl':114, 'Mc':115, 'Lv':116, 'Ts':117, 'Og':118
        }


class Visualization():

    from mayavi import mlab
    import numpy as np

    #import visualization.visualization_data as vd
    import visualization.utils as u

    input_type: str = None
    input_sub_type: str = None
    input_name: str = None
    file_string = None
    BAS_filename = None
    data_source = None

    data_input = None
    molecular_system = None
    orbital_generator = None

    visualization_data = VisualizationData()


    def __init__(self, input_type=None, input_sub_type=None, input_name=None, file_string=None, BAS_filename=None, data_source=None):

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

        #if self.input_type is not None and self.input_sub_type and self.input_name is not None:
        if self.input_type is not None and self.input_name is not None:
            
            self.initialize_data_input()

        self.initialize_molecular_system()

    def set_input_type(self, input_type):

        self.input_type = input_type

    def set_input_sub_type(self, input_sub_type):

        self.input_sub_type = input_sub_type

    def set_input_name(self, input_name):

        self.input_name = input_name

    def set_source(self, source):

        self.source = source

    def set_data_source(self, data_source):

        self.data_source = data_source

    def initialize_data_input(self):

        import visualization.inputs

        self.data_input = visualization.inputs.get_input( input_type=self.input_type, input_sub_type=self.input_sub_type, input_name=self.input_name)


    def plot_Geometry(self, plot_atoms = 1, atom_names = 1, plot_bonds = 1 , atom_scaling = 1.0, bond_scaling = 1.0, contours = 6, background_color = None, sclalarbar = False, auto_show = True, figure = None ):

        self.mlab.figure("Geometria", bgcolor=(.5, .5, .75), size=(1000, 1000))
        self.mlab.clf()

        if plot_atoms:
            for i in range(self.molecular_system.nAtoms):
                self.mlab.points3d(self.molecular_system.atoms_R[i, 0],
                                   self.molecular_system.atoms_R[i, 1],
                                   self.molecular_system.atoms_R[i, 2],
                                   scale_factor=self.visualization_data.Atoms_Scale[
                                       self.u.letters(self.molecular_system.atoms_Name[i])],
                                   resolution=20,
                                   color=self.visualization_data.Atoms_Color[
                                       self.u.letters(self.molecular_system.atoms_Name[i])],
                                   scale_mode='none')

        if atom_names:
            for i in range(self.molecular_system.nAtoms):
                self.mlab.text3d(self.molecular_system.atoms_R[i, 0],
                                 self.molecular_system.atoms_R[i, 1],
                                 self.molecular_system.atoms_R[i, 2],
                                 self.molecular_system.atoms_Name[i], scale=(.9, .9, .9))

        if plot_bonds:
            for i, bond in enumerate( self.molecular_system.bonds ) :

                bond_begin = self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(bond[0]) ]
                bond_end   = self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(bond[1]) ]

                bond_half = 0.5 * ( self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(bond[0]) ] + self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(bond[1]) ] )

                self.mlab.plot3d(
                    self.np.array([bond_begin[0] , bond_half[0]]),
                    self.np.array([bond_begin[1] , bond_half[1]]),
                    self.np.array([bond_begin[2] , bond_half[2]]),
                    tube_radius=0.2 * bond_scaling , 
                    color= self.visualization_data.Atoms_Color[  self.u.letters( bond[0] ) ])

                self.mlab.plot3d(
                    self.np.array([bond_end[0] , bond_half[0]]),
                    self.np.array([bond_end[1] , bond_half[1]]),
                    self.np.array([bond_end[2] , bond_half[2]]),
                    tube_radius=0.2 * bond_scaling , 
                    color= self.visualization_data.Atoms_Color[  self.u.letters( bond[1] ) ])

        self.mlab.show()

    def initialize_molecular_system(self):

        import visualization.molecular_system

        self.molecular_system = visualization.molecular_system.Molecular_System()

    def get_geometry(self, get_bonds = True):

        if self.data_input is None or self.molecular_system is None:

            print('Either self.data_input or self.Molecular_System')

        else:
            
            nAtoms=self.data_input.get_nAtoms()
            atoms_R, atoms_Charge, atoms_Name =self.data_input.get_Atoms()

            self.molecular_system.set_nAtoms(nAtoms=nAtoms)
            self.molecular_system.set_atoms_R(atoms_R = atoms_R )
            self.molecular_system.set_atoms_Charge(atoms_Charge=atoms_Charge)
            self.molecular_system.set_atoms_Name(atoms_Name=atoms_Name)   
            
            if get_bonds:
                bonds=self.data_input.get_Bonds()
                self.molecular_system.set_bonds(bonds=bonds )


    def get_orbital_data(self, get_bonds = True ):

        if self.data_input is None or self.molecular_system is None:

            print('Either self.data_input or self.Molecular_System')

        else:

            from visualization.orbitals import OrbitalsGenerator

            self.get_geometry( get_bonds = get_bonds)

            self.molecular_system.set_spherical( spherical=self.data_input.get_spherical())

            self.molecular_system.set_nb( nb=self.data_input.get_nb())

            self.molecular_system.set_Coeff(Coeff=self.data_input.get_Coeff())

            self.molecular_system.set_basis_and_norms(input=self.data_input.get_Basis())

            self.orbital_generator = OrbitalsGenerator( nAtoms = self.molecular_system.get_nAtoms(), 
                                                        atoms_R = self.molecular_system.get_atoms_R(), 
                                                        spherical=self.molecular_system.get_spherical(), 
                                                        nb =  self.molecular_system.get_nb(),
                                                        coeff=self.molecular_system.get_coeff(),
                                                        basis = self.molecular_system.get_basis(), 
                                                        basis_norm = self.molecular_system.get_basis_norm() )

    def generate_AO_orbitals(self):

        self.orbital_generator.calc_AOs( AO = self.orbital_generator.AOs )
        self.molecular_system.AOs = self.orbital_generator.AOs


    def generate_AO_orbitals_gpu(self):

        self.orbital_generator.calc_AOs_gpu( AO = self.orbital_generator.AOs )
        self.molecular_system.AOs = self.orbital_generator.AOs



#    @numba.cuda.jit
    def generate_AO_orbitals_numba(self):


        self.orbital_generator.calc_AOs( AO = self.orbital_generator.AOs )
        self.molecular_system.AOs = self.orbital_generator.AOs

    def generate_MO_orbitals(self, gpu=False):

        self.orbital_generator.calc_MOs( )
        self.molecular_system.MOs = self.orbital_generator.MOs

    def generate_MO_orbitals_gpu(self):

        self.orbital_generator.calc_MOs_gpu( )
        self.molecular_system.MOs = self.orbital_generator.MOs

    def generate_MO_orbitals_gpu_low_memory(self):

        self.orbital_generator.calc_MOs_gpu_low_memory( )
        self.molecular_system.MOs = self.orbital_generator.MOs


    def get_geminal_data(self):

        if self.data_input is None or self.molecular_system is None:

            print('Either self.data_input or self.Molecular_System')

        else:


            self.data_input.get_inactive()
            self.data_input.get_electrons()
            self.data_input.get_Occ()

            self.data_input.get_G_coeff()
            self.data_input.get_Orb2Gem()

            self.molecular_system.inactive = self.data_input.inactive
            self.molecular_system.G_coeff = self.data_input.G_coeff            
            self.molecular_system.Orb2Gem = self.data_input.Orb2Gem
            self.molecular_system.n_geminals = self.data_input.nGeminal

            self.geminal_generator = GeminalGenerator(  MOs = self.molecular_system.MOs, 
                                                        n_geminal = self.molecular_system.n_geminals, 
                                                        inactive = self.molecular_system.inactive, 
                                                        Orb2Gem = self.molecular_system.Orb2Gem, 
                                                        G_coeff = self.molecular_system.G_coeff, 
                                                        geminals = self.molecular_system.geminals  ) 

    def contur_process(self, contour, cube ):
        
        if type( contour ) == str:
            if contour[-1] == '%':

                contour_input_value = float(contour[:-1])

                cube_processed = self.np.sort(cube.flatten())[::-1]
                cube_csum = self.np.cumsum( cube_processed ) 
                cube_csum = cube_csum / cube_csum[-1]
            
                #print("cube_csum",  cube_csum )
                
                limit_value = cube_processed[cube_csum > ( (contour_input_value / 100 ) )][0]
                print(f"for contur {contour} limit value is {limit_value:.2e}")

                return [limit_value]

        if type( contour ) == list:

            return_list = []

            for contour_element in contour:
                return_list.append( self.contur_process(contour_element, cube) )

            return list(self.np.array(return_list).flat)

        else: 
            return contour

        #if type( contour ) == int:
        #    return contour

        #if type( contour ) == list:
        #    return contour


        

    def get_generate_geminals(self):

        self.molecular_system.geminals = self.geminal_generator.Calc_Geminals( )
        self.molecular_system.geminals = self.geminal_generator.geminals


    def plot_Orbitals(self, orbital_numbers, plot_atoms=1, atom_names=1, plot_bonds=1):

        fig1 = self.mlab.figure("Orbitals_" + str(orbital_numbers), bgcolor=(.5, .5, .75), size=(1000, 1000))
        fig1.scene.parallel_projection = False

        self.mlab.clf()

        if plot_atoms:
            for i in range(self.molecular_system.nAtoms):
                self.mlab.points3d(self.molecular_system.atoms_R[i, 0],
                                   self.molecular_system.atoms_R[i, 1],
                                   self.molecular_system.atoms_R[i, 2],
                                   scale_factor=self.visualization_data.Atoms_Scale[ self.u.letters(self.molecular_system.atoms_Name[i])],
                                   resolution=20,
                                   color=self.visualization_data.Atoms_Color[ self.u.letters(self.molecular_system.atoms_Name[i])],
                                   scale_mode='none')

        if atom_names:
            for i in range(self.molecular_system.nAtoms):
                self.mlab.text3d(self.molecular_system.atoms_R[i, 0],
                                 self.molecular_system.atoms_R[i, 1],
                                 self.molecular_system.atoms_R[i, 2],
                                 self.molecular_system.atoms_Name[i], scale=(.9, .9, .9))

        if plot_bonds:
            for i in range(len(self.molecular_system.bonds)):
                self.mlab.plot3d(
                    self.np.array([ self.molecular_system.atoms_R[
                                    self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][0]), 0],
                                    self.molecular_system.atoms_R[
                                    self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][1]), 0]]),
                    self.np.array([ self.molecular_system.atoms_R[
                                    self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][0]), 1],
                                    self.molecular_system.atoms_R[
                                    self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][1]), 1]]),
                    self.np.array([ self.molecular_system.atoms_R[
                                    self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][0]), 2],
                                    self.molecular_system.atoms_R[
                                    self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][1]), 2]]),
                    tube_radius=0.2, color=(233.0 / 255, 165.0 / 255, 165.0 / 255))

        for number in orbital_numbers:
            
            #X, Y, Z = self.molecular_system.grid.return_grid_arrays()
            X, Y, Z = self.orbital_generator.grid.return_grid_arrays()
            self.mlab.contour3d( X, Y, Z, (self.molecular_system.MOs[number]), contours=12, opacity=0.5)


        self.mlab.show()

    def plot_Orbitals_AO(self, orbital_numbers, plot_atoms=1, atom_names=1, plot_bonds=1):

        fig1 = self.mlab.figure("Orbitals_AO_" + str(orbital_numbers), bgcolor=(.5, .5, .75), size=(1000, 1000))
        fig1.scene.parallel_projection = False

        self.mlab.clf()

        if plot_atoms:
            for i in range(self.molecular_system.nAtoms):
                self.mlab.points3d(self.molecular_system.atoms_R[i, 0],
                                   self.molecular_system.atoms_R[i, 1],
                                   self.molecular_system.atoms_R[i, 2],
                                   scale_factor=self.visualization_data.Atoms_Scale[ self.u.letters(self.molecular_system.atoms_Name[i])],
                                   resolution=20,
                                   color=self.visualization_data.Atoms_Color[ self.u.letters(self.molecular_system.atoms_Name[i])],
                                   scale_mode='none')

        if atom_names:
            for i in range(self.molecular_system.nAtoms):
                self.mlab.text3d(self.molecular_system.atoms_R[i, 0],
                                 self.molecular_system.atoms_R[i, 1],
                                 self.molecular_system.atoms_R[i, 2],
                                 self.molecular_system.atoms_Name[i], scale=(.9, .9, .9))

        if plot_bonds:
            for i in range(len(self.molecular_system.bonds)):
                self.mlab.plot3d(
                    self.np.array([ self.molecular_system.atoms_R[
                                    self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][0]), 0],
                                    self.molecular_system.atoms_R[
                                    self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][1]), 0]]),
                    self.np.array([ self.molecular_system.atoms_R[
                                    self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][0]), 1],
                                    self.molecular_system.atoms_R[
                                    self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][1]), 1]]),
                    self.np.array([ self.molecular_system.atoms_R[
                                    self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][0]), 2],
                                    self.molecular_system.atoms_R[
                                    self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][1]), 2]]),
                    tube_radius=0.2, color=(233.0 / 255, 165.0 / 255, 165.0 / 255))

        for number in orbital_numbers:
            
            #X, Y, Z = self.molecular_system.grid.return_grid_arrays()
            X, Y, Z = self.orbital_generator.grid.return_grid_arrays()
            
            self.mlab.contour3d( X, Y, Z, (self.molecular_system.AOs[number]), contours=12, opacity=0.5)


        self.mlab.show()

    def plot_Geminals(self, geminal_numbers = [0],  
                            plot_atoms = True, 
                            atom_names =True, 
                            plot_bonds = True , 
                            atom_scaling = 1.0, 
                            bond_scaling = 1.0, 
                            contours = 6, 
                            background_color = None, 
                            sclalarbar = False, 
                            auto_show = True, 
                            figure = None ):  
        
        if background_color is None:
            _background_color = self.visualization_data.background_colors['White']
        else:
            _background_color = background_color

        if figure is None:
            _figure = self.mlab.figure(   "Dispersion", 
                            bgcolor=_background_color,
                            size=(600, 600) )

            self.mlab.clf()

        else:

            _figure = figure

        if plot_atoms:
            for i in range(self.molecular_system.nAtoms):
                self.mlab.points3d(self.molecular_system.atoms_R[i, 0],
                                   self.molecular_system.atoms_R[i, 1],
                                   self.molecular_system.atoms_R[i, 2],
                                   scale_factor= atom_scaling * self.visualization_data.Atoms_Scale[ self.u.letters(self.molecular_system.atoms_Name[i])],
                                   resolution=20,
                                   color=self.visualization_data.Atoms_Color[ self.u.letters(self.molecular_system.atoms_Name[i])],
                                   scale_mode='none')

        if atom_names:
            for i in range(self.molecular_system.nAtoms):
                self.mlab.text3d(self.molecular_system.atoms_R[i, 0],
                                 self.molecular_system.atoms_R[i, 1],
                                 self.molecular_system.atoms_R[i, 2],
                                 self.molecular_system.atoms_Name[i], scale=(.9, .9, .9))

        if plot_bonds:
            for i, bond in enumerate( self.molecular_system.bonds ) :

                bond_begin = self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(bond[0]) ]
                bond_end   = self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(bond[1]) ]

                bond_half = 0.5 * ( self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(bond[0]) ] + self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(bond[1]) ] )

                self.mlab.plot3d(
                    self.np.array([bond_begin[0] , bond_half[0]]),
                    self.np.array([bond_begin[1] , bond_half[1]]),
                    self.np.array([bond_begin[2] , bond_half[2]]),
                    tube_radius=0.2 * bond_scaling , 
                    color= self.visualization_data.Atoms_Color[  self.u.letters( bond[0] ) ])

                self.mlab.plot3d(
                    self.np.array([bond_end[0] , bond_half[0]]),
                    self.np.array([bond_end[1] , bond_half[1]]),
                    self.np.array([bond_end[2] , bond_half[2]]),
                    tube_radius=0.2 * bond_scaling , 
                    color= self.visualization_data.Atoms_Color[  self.u.letters( bond[1] ) ])

        for number in geminal_numbers:

            #X, Y, Z = self.molecular_system.grid.return_grid_arrays()
            X, Y, Z = self.orbital_generator.grid.return_grid_arrays()

            #          exec(self.Geminnals[i].g)
            #self.mlab.contour3d(X, Y, Z, (self.Geminnals[number].v), contours=12, opacity=0.5)
            self.mlab.contour3d( X, Y, Z, (self.molecular_system.geminals[number]), 
                                contours= self.contur_process( contours, self.molecular_system.geminals[number] ) , opacity=0.5)

        if auto_show:
            self.mlab.show()



    def _plot_bonds(self, bond_scaling):

        for i, bond in enumerate( self.molecular_system.bonds ) :

            bond_begin = self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(bond[0]) ]
            bond_end   = self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(bond[1]) ]

            bond_half = 0.5 * ( self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(bond[0]) ] + self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(bond[1]) ] )

            self.mlab.plot3d(   self.np.array([bond_begin[0] , bond_half[0]]),
                                self.np.array([bond_begin[1] , bond_half[1]]),
                                self.np.array([bond_begin[2] , bond_half[2]]),
                                tube_radius=0.2 * bond_scaling ,                     
                                color= self.visualization_data.Atoms_Color[  self.u.letters( bond[0] ) ])

            self.mlab.plot3d(   self.np.array([bond_end[0] , bond_half[0]]),
                                self.np.array([bond_end[1] , bond_half[1]]),
                                self.np.array([bond_end[2] , bond_half[2]]),
                                tube_radius=0.2 * bond_scaling , 
                                color= self.visualization_data.Atoms_Color[  self.u.letters( bond[1] ) ])

    def _plot_atoms(self, atom_scaling):
            
        for i in range(self.molecular_system.nAtoms):
            self.mlab.points3d( self.molecular_system.atoms_R[i, 0],
                                self.molecular_system.atoms_R[i, 1],
                                self.molecular_system.atoms_R[i, 2],
                                scale_factor= atom_scaling * self.visualization_data.Atoms_Scale[ self.u.letters(self.molecular_system.atoms_Name[i])],
                                resolution=20,
                                color=self.visualization_data.Atoms_Color[ self.u.letters(self.molecular_system.atoms_Name[i])],
                                scale_mode='none')

    def _atom_names(self, atom_names_scaling):
            for i in range(self.molecular_system.nAtoms):
                self.mlab.text3d(self.molecular_system.atoms_R[i, 0],
                                 self.molecular_system.atoms_R[i, 1],
                                 self.molecular_system.atoms_R[i, 2],
                                 self.molecular_system.atoms_Name[i], 
                                 scale=(.9*atom_names_scaling, .9*atom_names_scaling, .9*atom_names_scaling)  )

    def plot_orbitals_MO(self,  orbital_numbers = [0],  
                                plot_atoms = True, 
                                atom_scaling = 1.0, 
                                atom_names =True, 
                                atom_names_scaling = 1.0, 
                                plot_bonds = True , 
                                bond_scaling = 1.0, 
                                contours = 6, 
                                background_color = None, 
                                sclalarbar = False, 
                                auto_show = True, 
                                figure = None ):  
        
        if background_color is None:
            _background_color = self.visualization_data.background_colors['White']
        else:
            _background_color = background_color

        if figure is None:
            _figure = self.mlab.figure( "Dispersion", 
                                        bgcolor=_background_color,
                                        size=(600, 600) )

            self.mlab.clf()

        else:

            _figure = figure

        if plot_atoms:
            self._plot_atoms(atom_scaling)

        if atom_names:
            self._atom_names(atom_names_scaling)

        if plot_bonds:
            self._plot_bonds(bond_scaling)

        for number in orbital_numbers:

            X, Y, Z = self.orbital_generator.grid.return_grid_arrays()
            self.mlab.contour3d( X, Y, Z, self.molecular_system.MOs[number], 
                                contours= self.contur_process( contours, self.molecular_system.MOs[number] ) , opacity=0.5)

        if auto_show:
            self.mlab.show()

        return _figure

    def plot_orbitals_MO_square(self,  orbital_numbers = [0],  
                                plot_atoms = True, 
                                atom_scaling = 1.0, 
                                atom_names =True, 
                                atom_names_scaling = 1.0, 
                                plot_bonds = True , 
                                bond_scaling = 1.0, 
                                contours = 6, 
                                background_color = None, 
                                sclalarbar = False, 
                                auto_show = True, 
                                figure = None ):  
        
        if background_color is None:
            _background_color = self.visualization_data.background_colors['White']
        else:
            _background_color = background_color

        if figure is None:
            _figure = self.mlab.figure( "Dispersion", 
                                        bgcolor=_background_color,
                                        size=(600, 600) )

            self.mlab.clf()

        else:

            _figure = figure

        if plot_atoms:
            self._plot_bonds(atom_scaling)

        if atom_names:
            self._atom_names(atom_names_scaling)

        if plot_bonds:
            self._plot_bonds(bond_scaling)

        for number in orbital_numbers:

            X, Y, Z = self.orbital_generator.grid.return_grid_arrays()
            self.mlab.contour3d( X, Y, Z, self.molecular_system.MOs[number]**2, 
                                contours= self.contur_process( contours, self.molecular_system.MOs[number]**2 ) , opacity=0.5)

        if auto_show:
            self.mlab.show()

        return _figure

    def get_geminals(self, R_max_multip = 3.0, x_n= 50, y_n= 50, z_n= 50, gpu=False):

        self.get_orbital_data()
        self.orbital_generator.grid.R_max_multip = R_max_multip
        self.orbital_generator.grid.x_n= x_n 
        self.orbital_generator.grid.y_n= y_n
        self.orbital_generator.grid.z_n= z_n
        self.orbital_generator.init_grid( )

        self.orbital_generator.init_AOs()

        if gpu:
            self.generate_AO_orbitals_gpu()            
            self.generate_MO_orbitals_gpu()

        else:
            self.generate_AO_orbitals()            
            self.generate_MO_orbitals()

        self.get_geminal_data()
        self.get_generate_geminals()

    def dummy_plot(self):
        return 0

    def plot_any_data(self,  any_data_to_plot=None,  
                                plot_atoms = True, 
                                atom_scaling = 1.0, 
                                atom_names =True, 
                                atom_names_scaling = 1.0, 
                                plot_bonds = True , 
                                bond_scaling = 1.0, 
                                contours = 6, 
                                background_color = None, 
                                sclalarbar = False, 
                                auto_show = True,
                                opacity = 0.5,
                                figure = None ):  
        
        if background_color is None:
            _background_color = self.visualization_data.background_colors['White']
        else:
            _background_color = background_color

        if figure is None:
            _figure = self.mlab.figure( "Dispersion", 
                                        bgcolor=_background_color,
                                        size=(600, 600) )

            self.mlab.clf()

        else:

            _figure = figure

        if plot_atoms:
            self._plot_bonds(atom_scaling)

        if atom_names:
            self._atom_names(atom_names_scaling)

        if plot_bonds:
            self._plot_bonds(bond_scaling)


        X, Y, Z = self.orbital_generator.grid.return_grid_arrays()
        self.mlab.contour3d( X, Y, Z, any_data_to_plot, 
                            contours= self.contur_process( contours, any_data_to_plot ) , opacity=opacity)

        if auto_show:
            self.mlab.show()

        return _figure
