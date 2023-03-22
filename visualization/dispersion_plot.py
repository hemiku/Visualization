
import visualization.visualization


class DispersionPlot(visualization.visualization.Visualization):

    GAMMCOR_output_filename = None
    GAMMCOR_output = None

    number_of_fragments = None

    nFragments = None
        
    Orb:np.ndarray = None
    Occ:np.ndarray = None
    Fragments:np.ndarray = None
    Monomers:np.ndarray = None

    Eps_A:np.ndarray = None
    Eps_AB:np.ndarray = None    
    
    D_AB:np.ndarray = None
        

    def set_GAMMCOR_filename(self, filename):
        
        self.GAMMCOR_output_filename = filename

    def get_GAMMCOR_Output(self):
        
        f=open(self.GAMMCOR_output_filename,'r')
        Out=f.read()
        f.close()
        
        return Out

    def get_NumberOfFragments(self):
        
        Out = self.get_GAMMCOR_Output( )
        
        self.NumberOfFragments = int(Out[ Out.find("NUMBER OF FRAGMENTS:"):Out.find("NUMBER OF FRAGMENTS:")+50 ].split()[3]) 
    
    def get_Fragments(self):
        
        
        nParts =  2 * self.MolecularSystem.n_geminals + self.MolecularSystem.inactive
        
        Orb = self.np.zeros([ nParts ], dtype=self.np.int32)
        Occ = self.np.zeros([ nParts ], dtype=self.np.float64)        
        Fragments = self.np.zeros([ nParts ], dtype=self.np.int32)
        Monomers = self.np.zeros([ nParts ], dtype=self.np.int32)

        
        Out = self.get_GAMMCOR_Output( )
        
        buf = Out[ Out.find("NUMBER OF FRAGMENTS:"):Out.find("GVB ONE-BODY CORRELATION ENERGY")].splitlines()[2:]
        for i in range( nParts ):
            
            Orb[i] = int( buf[i].split()[0] ) -1 #- self.MolecularSystem.inactive
            Occ[i]         = float( buf[i].split()[1] )
            Fragments[i] = int( buf[i].split()[2] ) -1  #- self.MolecularSystem.inactive
            Monomers[i] = int( buf[i].split()[3] ) -1 
            
        self.nParts = nParts 
        
        self.Orb  = Orb-1
        self.Occ  = Occ
        self.Fragments  = Fragments
        self.Monomers  = Monomers
    
    def get_Epsilons(self):
        
        Eps_A  = self.np.zeros([self.nParts ,self.nParts], dtype=self.np.float32)
        Eps_AB = self.np.zeros([self.nParts ,self.nParts], dtype=self.np.float32)
        
        Out = self.get_GAMMCOR_Output( )
        
        Out_Eps_A = Out[ Out.find("*** 2-body correlation for each pair of fragments"):Out.find("*** Total 2-body correlation")].splitlines()[1:-2]
        Out_Eps_AB = Out[ Out.find("*** 2-body correlation for fragments in different monomers"):Out.find("CPU  time in AC1")].splitlines()[1:-2]
        
        for i in Out_Eps_A:
            
            try:
                p =  int(i.split()[0]) - 1 
                q =  int(i.split()[1]) - 1 
                Eps_A[p , q ] = float(i.split()[2])
                Eps_A[q , p ] = float(i.split()[2])
            except:
                print("break at line:", i  )
                break
            
        for i in Out_Eps_AB:
            try:
                p =  int(i.split()[0]) - 1 
                q =  int(i.split()[1]) - 1 
            
                Eps_AB[p , q ] = float(i.split()[2])
                Eps_AB[q , p ] = float(i.split()[2])
            
            except:
                print("break at line:", i  )
                break

        self.Eps_A = Eps_A
        self.Eps_AB = Eps_AB    


    def Calc_D(self, Monomer_A, Monomer_B ):

        import logging

        D_A = self.np.zeros_like(self.MolecularSystem.AOs[0])
        D_B = self.np.zeros_like(self.MolecularSystem.AOs[0])

        grid = self.orbital_generator.grid

        dx = ( grid.x_max - grid.x_min )/ ( grid.x_n -1 )
        dy = ( grid.y_max - grid.y_min )/ ( grid.y_n -1 )
        dz = ( grid.z_max - grid.z_min )/ ( grid.z_n -1 )

        dv = dx * dy * dz

        fragments_values, fragments_counts = self.np.unique(self.Fragments, return_counts=True)


        fragments_multiplier = self.np.array( [ [1, 1,  1,  1,  1,  1,  1],
                                                [1, 1,  1,  1,  1,  1,  1],
                                                [1, 1,  0.5,1,  1,  1,  1],
                                                [1, 1,  1,  1,  1,  1,  1],
                                                [1, 1,  1,  1,  1,  1,  1],
                                                [1, 1,  1,  1,  1,  1,  1],
                                                [1, 1,  1,  1,  1,  1,  1],
                                                [1, 1,  1,  1,  1,  1,  1] ])

        counter_6_2 = 0
        counter_6_6 = 0
        counter = 0

        for i in range(self.nParts):
            for j in range(self.nParts):

                if  ( fragments_counts[self.Fragments[i]] == 2 and fragments_counts[self.Fragments[j]] == 2 ):
                    multipier = 1.0/2
                
                elif  ( ( fragments_counts[self.Fragments[i]] == 4 and fragments_counts[self.Fragments[j]] == 2 ) or 
                        ( fragments_counts[self.Fragments[j]] == 4 and fragments_counts[self.Fragments[i]] == 2 ) ):
                    multipier = 1.0/4
                
                elif  ( fragments_counts[self.Fragments[i]] == 4 and fragments_counts[self.Fragments[j]] == 4 ):
                    multipier = 1.0/8

                elif  ( fragments_counts[self.Fragments[i]] == 6 and fragments_counts[self.Fragments[j]] == 6 ):
                    counter_6_6 += 1
                    multipier = 1.0/18

                elif  ( ( fragments_counts[self.Fragments[i]] == 6 and fragments_counts[self.Fragments[j]] == 2 ) or 
                        ( fragments_counts[self.Fragments[j]] == 6 and fragments_counts[self.Fragments[i]] == 2 ) ):
                    counter_6_2 += 1
                    multipier = 1.0/6

                else:
                    print("Undefnied proportion", "Fragment A:", self.Fragments[i],fragments_counts[self.Fragments[i]] , "Fragment B:", self.Fragments[j],fragments_counts[self.Fragments[j]]  )
                    counter += 1
                    multipier = 1.0


                if self.Monomers[i] == Monomer_A and self.Monomers[j] == Monomer_B :

                    epsilon = self.Eps_AB[self.Fragments[i],self.Fragments[j]]

                    if epsilon != 0.0:
                        logging.info(f'A: {i},{j}, {self.Orb[i]}, {epsilon} ')
                        print(f'A: {i},{j}, {self.Orb[i]}, {epsilon} ')
                    D_A -= multipier * self.Occ[i] * self.MolecularSystem.MOs[i]**2 * epsilon
                
                if self.Monomers[i] == Monomer_B and self.Monomers[j] == Monomer_A:

                    epsilon = self.Eps_AB[self.Fragments[i],self.Fragments[j]]

                    if epsilon != 0.0:
                        logging.info(f'B: {i},{j}, {self.Orb[i]}, {epsilon} ')
                        print(f'B: {i},{j}, {self.Orb[i]}, {epsilon} ')
                    D_B -= multipier * self.Occ[i] * self.MolecularSystem.MOs[i]**2 * epsilon

        print("counter 6 6 :",counter_6_6)
        print("counter 6 2 :",counter_6_2)
        print("counter:",counter)

        D_AB  = 0.5 *( D_A + D_B )
        
        self.D_AB = D_AB
        self.D_A = D_A
        self.D_B = D_B


    def Plot_D_AB(self, plot_atoms = 1, atom_names = 1, plot_bonds = 1 , atom_scaling = 1.0, bond_scaling = 1.0, contours = 6, opacity = 0.5, background_color = None, sclalarbar = False, auto_show = True, figure = None ):  
        
        #from mayavi import mlab

        if background_color is None:
            _background_color = self.visualization_data.background_colors['White']
        else:
            _background_color = background_color

        if figure is None:
            _figure = self.mlab.figure(   "Dispersion", 
                            bgcolor=_background_color,
                            size=(600, 400) )
            self.mlab.clf()

        else:

            _figure = figure

        if plot_atoms:
            self._plot_atoms(atom_scaling)

        if atom_names:
            self._atom_names( atom_names_scaling )

        if plot_bonds:
            self._plot_bonds( plot_bonds, bond_scaling)

        X, Y, Z = self.orbital_generator.grid.return_grid_arrays()

        D_AB_contur = self.mlab.contour3d(  X, Y, Z, (self.D_AB) , 
                                            contours=self.contour_process( contours, self.D_AB ), 
                                            opacity=opacity) 
    
        if sclalarbar:

            _scalarbar = self.mlab.scalarbar(object=D_AB_contur, title='D^{AB}' ,orientation='vertical' )
            _scalarbar.scalar_bar.unconstrained_font_size = True
            _scalarbar.scalar_bar.label_text_property.color = (0.0, 0.0, 0.0) 
            _scalarbar.scalar_bar.label_text_property.italic = False
            _scalarbar.scalar_bar.label_text_property.bold = False 
            _scalarbar.scalar_bar.number_of_labels = contours
            _scalarbar.scalar_bar.label_text_property.font_size = 14

            _scalarbar.scalar_bar.position = self.np.array([0.01      , 0.01])
            _scalarbar.scalar_bar.position2 = self.np.array([0.2       , 0.95])

        if auto_show:
            self.mlab.show()

        return _figure

    def Plot_D_A(self, plot_atoms = 1, atom_names = 1, plot_bonds = 1 , atom_scaling = 1.0, bond_scaling = 1.0, contours = 6, background_color = None, sclalarbar = False, auto_show = True, figure = None ):  
        
        #from mayavi import mlab

        if background_color is None:
            _background_color = self.visualization_data.background_colors['White']
        else:
            _background_color = background_color

        if figure is None:
            _figure = self.mlab.figure(   "Dispersion", 
                            bgcolor=_background_color,
                            size=(600, 400) )
            self.mlab.clf()

        else:

            _figure = figure

        if plot_atoms:
            for i in range(self.MolecularSystem.nAtoms):
                self.mlab.points3d( self.MolecularSystem.atoms_R[i, 0],
                                    self.MolecularSystem.atoms_R[i, 1],
                                    self.MolecularSystem.atoms_R[i, 2],
                                    scale_factor= atom_scaling * self.visualization_data.Atoms_Scale[ self.u.letters(self.MolecularSystem.atoms_Name[i])],
                                    resolution=20,
                                    color=self.visualization_data.Atoms_Color[ self.u.letters(self.MolecularSystem.atoms_Name[i])],
                                    scale_mode='none')
            
        if atom_names:
            for i in range( self.MolecularSystem.nAtoms ):
                self.mlab.text3d(   self.MolecularSystem.atoms_R[i, 0],
                                    self.MolecularSystem.atoms_R[i, 1],
                                    self.MolecularSystem.atoms_R[i, 2],
                                    self.MolecularSystem.atoms_Name[i], scale=(.9*atom_scaling , .9*atom_scaling , .9*atom_scaling ))
        
        if plot_bonds:        
            for i, bond in enumerate( self.MolecularSystem.bonds ) :

                bond_begin = self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(bond[0]) ]
                bond_end   = self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(bond[1]) ]

                bond_half = 0.5 * ( self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(bond[0]) ] + self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(bond[1]) ] )

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

        X, Y, Z = self.orbital_generator.grid.return_grid_arrays()

        D_A_contur = self.mlab.contour3d(X, Y, Z,   (self.D_A) ,
                                                    contours=self.contour_process( contours, self.D_A ),
                                                    opacity=0.5) 

        if auto_show:
            self.mlab.show()

        return _figure

    def Plot_D_B(self, plot_atoms = 1, atom_names = 1, plot_bonds = 1 , atom_scaling = 1.0, bond_scaling = 1.0, contours = 6, background_color = None, sclalarbar = False, auto_show = True, figure = None ):  
        
        #from mayavi import mlab

        if background_color is None:
            _background_color = self.visualization_data.background_colors['White']
        else:
            _background_color = background_color

        if figure is None:
            _figure = self.mlab.figure(   "Dispersion", 
                            bgcolor=_background_color,
                            size=(600, 400) )
            self.mlab.clf()

        else:

            _figure = figure

        if plot_atoms:
            for i in range(self.MolecularSystem.nAtoms):
                self.mlab.points3d( self.MolecularSystem.atoms_R[i, 0],
                                    self.MolecularSystem.atoms_R[i, 1],
                                    self.MolecularSystem.atoms_R[i, 2],
                                    scale_factor= atom_scaling * self.visualization_data.Atoms_Scale[ self.u.letters(self.MolecularSystem.atoms_Name[i])],
                                    resolution=20,
                                    color=self.visualization_data.Atoms_Color[ self.u.letters(self.MolecularSystem.atoms_Name[i])],
                                    scale_mode='none')
            
        if atom_names:
            for i in range( self.MolecularSystem.nAtoms ):
                self.mlab.text3d(   self.MolecularSystem.atoms_R[i, 0],
                                    self.MolecularSystem.atoms_R[i, 1],
                                    self.MolecularSystem.atoms_R[i, 2],
                                    self.MolecularSystem.atoms_Name[i], scale=(.9*atom_scaling , .9*atom_scaling , .9*atom_scaling ))
        
        if plot_bonds:        
            for i, bond in enumerate( self.MolecularSystem.bonds ) :

                bond_begin = self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(bond[0]) ]
                bond_end   = self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(bond[1]) ]

                bond_half = 0.5 * ( self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(bond[0]) ] + self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(bond[1]) ] )

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

        X, Y, Z = self.orbital_generator.grid.return_grid_arrays()

        D_B_contur = self.mlab.contour3d(X, Y, Z, (self.D_B) ,
                                                    contours=self.contour_process( contours, self.D_B ),
                                                    opacity=0.5) 
        if auto_show:
            self.mlab.show()

        return _figure

    def Plot_D_AB_log(self, Plot_Atoms = 1, Atom_Names = 1, Plot_Bonds = 1 , Atom_Scaling = 1.0, Bond_Scaling = 1.0, contours = 6 ):  

        return _ 
        #from mayavi import mlab
        
        # self.mlab.figure("Dispersion ", bgcolor=(.5, .5, .75), size=(1000, 1000))
        # self.mlab.clf()
        
        # if Plot_Atoms:
        #     for i in xrange(self.MolecularSystem.nAtoms):
        #         self.mlab.points3d(self.MolecularSystem.atoms_R[i,0], self.MolecularSystem.atoms_R[i,1], self.MolecularSystem.atoms_R[i,2], 
        #             scale_factor=self.MolecularSystem.Atoms_Scale[letters(self.MolecularSystem.atoms_Name[i])] * Atom_Scaling  , resolution=20, 
        #             color=self.MolecularSystem.Atoms_Color[letters(self.MolecularSystem.atoms_Name[i])], scale_mode='none')
            
        # if Atom_Names:
        #     for i in xrange(self.MolecularSystem.nAtoms):
        #         self.mlab.text3d(self.MolecularSystem.atoms_R[i,0], self.MolecularSystem.atoms_R[i,1], self.MolecularSystem.atoms_R[i,2], 
        #                     self.MolecularSystem.atoms_Name[i],scale=(.2, .2, .2))
        
        # if Plot_Bonds:        
        #     for i in xrange(len(self.MolecularSystem.bonds)):
        #         self.mlab.plot3d(
        #             self.np.array([self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(self.MolecularSystem.bonds[i][0]), 0] , 
        #                            self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(self.MolecularSystem.bonds[i][1]), 0]]),
        #             self.np.array([self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(self.MolecularSystem.bonds[i][0]), 1] , 
        #                            self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(self.MolecularSystem.bonds[i][1]), 1]]), 
        #             self.np.array([self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(self.MolecularSystem.bonds[i][0]), 2] , 
        #                            self.MolecularSystem.atoms_R[ self.MolecularSystem.atoms_Name.index(self.MolecularSystem.bonds[i][1]), 2]]),
        #             tube_radius=0.2 * Bond_Scaling , color=(233.0/255, 165.0/255, 165.0/255))
        
        # X, Y, Z = self.orbital_generator.grid.return_grid_arrays()

        # self.mlab.contour3d(X, Y, Z, self.np.log(self.D_AB), contours = contours, opacity=0.5, vmax = 4.8312953302547844e-05) 
    
        # self.mlab.show()
    

    

    def get_dispersion_index(self, R_max_multip = 3.0, x_n= 50, y_n= 50, z_n= 50, monomer_A=None, monomer_B=None, get_bonds=True, gpu=False):


        self.get_geometry( get_bonds = get_bonds)
        self.get_geminals( R_max_multip = R_max_multip, x_n= x_n, y_n= y_n, z_n= z_n, gpu=gpu)

        self.get_NumberOfFragments()
        self.get_Fragments()
        self.get_Epsilons()
        Monomers_set = set( self.Monomers )

        if monomer_A is None and monomer_B is None: 
            self.Calc_D( Monomer_A= Monomers_set.pop(), Monomer_B=Monomers_set.pop() )
        else:
            self.Calc_D( Monomer_A= monomer_A, Monomer_B=monomer_B )

    def print_something(self):

        print("orint something")