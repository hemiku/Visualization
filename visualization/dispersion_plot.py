
import visualization.visualization


class DispersionPlot(visualization.visualization.Visualization):

    GAMMCOR_output_filename = None
    GAMMCOR_output = None

    number_of_fragments = None

    nFragments = None
        
    Orb = None
    Occ = None
    Fragments = None
    Monomers = None

    Eps_A = None
    Eps_AB = None    
    
    D_AB = None
        

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
        
        
        nFragments =  2 * self.molecular_system.n_geminals + self.molecular_system.inactive
        
        Orb = self.np.zeros([ nFragments ], dtype=self.np.int32)
        Occ = self.np.zeros([ nFragments ], dtype=self.np.float64)        
        Fragments = self.np.zeros([ nFragments ], dtype=self.np.int32)
        Monomers = self.np.zeros([ nFragments ], dtype=self.np.int32)

        
        Out = self.get_GAMMCOR_Output( )
        
        buf = Out[ Out.find("NUMBER OF FRAGMENTS:"):Out.find("GVB ONE-BODY CORRELATION ENERGY")].splitlines()[2:]
        for i in range( nFragments ):
            
            Orb[i] = int( buf[i].split()[0] ) #- self.molecular_system.inactive
            Occ[i]         = float( buf[i].split()[1] )
            Fragments[i] = int( buf[i].split()[2] ) #- self.molecular_system.inactive
            Monomers[i] = int( buf[i].split()[3] )
            
        self.nFragments = nFragments 
        
        self.Orb  = Orb
        self.Occ  = Occ
        self.Fragments  = Fragments
        self.Monomers  = Monomers
    
    def get_Epsilons(self):
        
        Eps_A  = self.np.zeros([self.nFragments ,self.nFragments], dtype=self.np.float32)
        Eps_AB = self.np.zeros([self.nFragments ,self.nFragments], dtype=self.np.float32)
        
        Out = self.get_GAMMCOR_Output( )
        
        Out_Eps_A = Out[ Out.find("*** 2-body correlation for each pair of fragments"):Out.find("*** Total 2-body correlation")].splitlines()[1:-2]
        Out_Eps_AB = Out[ Out.find("*** 2-body correlation for fragments in different monomers"):Out.find("CPU  time in AC1")].splitlines()[1:-4]
        
        for i in Out_Eps_A:
            p =  int(i.split()[0])
            q =  int(i.split()[1])
            Eps_A[p , q ] = float(i.split()[2])
            Eps_A[q , p ] = float(i.split()[2])
        
#        print Out_Eps_AB
        for i in Out_Eps_AB:
#            print p,q
            p =  int(i.split()[0])
            q =  int(i.split()[1])
            
            Eps_AB[p , q ] = float(i.split()[2])
            Eps_AB[q , p ] = float(i.split()[2])
            
        self.Eps_A = Eps_A
        self.Eps_AB = Eps_AB    


    def Calc_D(self, Monomer_A, Monomer_B ):

        D_A = self.np.zeros_like(self.molecular_system.AOs[0])
        D_B = self.np.zeros_like(self.molecular_system.AOs[0])

        for i in range(self.nFragments):
            for j in range(self.nFragments):
#            for j in xrange( i ):
                if self.Monomers[i] == Monomer_A and self.Monomers[j] == Monomer_B:

                    D_A -= self.Occ[i] * self.molecular_system.MOs[i]**2 * self.Eps_AB[self.Fragments[i],self.Fragments[j]]
                
                if self.Monomers[i] == Monomer_B and self.Monomers[j] == Monomer_A:
                    
                    D_B -= self.Occ[i] * self.molecular_system.MOs[i]**2 * self.Eps_AB[self.Fragments[i],self.Fragments[j]]

        D_AB  = 0.5 *( D_A + D_B )
        
        self.D_AB = D_AB
        

    def Plot_D_AB(self, plot_atoms = 1, atom_names = 1, plot_bonds = 1 , atom_scaling = 1.0, bond_scaling = 1.0, contours = 6, background_color = None, sclalarbar = False, auto_show = True, figure = None ):  
        
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

        #if Plot_Atoms:
        #    for i in range(self.molecular_system.nAtoms):
        #        self.mlab.points3d(self.molecular_system.atoms_R[i,0], self.molecular_system.atoms_R[i,1], self.molecular_system.atoms_R[i,2], 
        #            scale_factor=self.visualization_data.Atoms_Scale[letters(self.molecular_system.atoms_Name[i])] * Atom_Scaling  , resolution=20, 
        #            color=self.visualization_data.Atoms_Color[letters(self.molecular_system.atoms_Name[i])], scale_mode='none')
        if plot_atoms:
            for i in range(self.molecular_system.nAtoms):
                self.mlab.points3d(self.molecular_system.atoms_R[i, 0],
                                   self.molecular_system.atoms_R[i, 1],
                                   self.molecular_system.atoms_R[i, 2],
                                   scale_factor= atom_scaling * self.visualization_data.Atoms_Scale[ self.u.letters(self.molecular_system.atoms_Name[i])],
                                   resolution=20,
                                   color=self.visualization_data.Atoms_Color[ self.u.letters(self.molecular_system.atoms_Name[i])],
                                   scale_mode='none')
            
        #if Atom_Names:
        #    for i in range(self.molecular_system.nAtoms):
        #        self.mlab.text3d(self.molecular_system.atoms_R[i,0], self.molecular_system.atoms_R[i,1], self.molecular_system.atoms_R[i,2], 
        #                    self.molecular_system.atoms_Name[i],scale=(.2, .2, .2))
        
        if atom_names:
            for i in range( self.molecular_system.nAtoms ):
                self.mlab.text3d(self.molecular_system.atoms_R[i, 0],
                                 self.molecular_system.atoms_R[i, 1],
                                 self.molecular_system.atoms_R[i, 2],
                                 self.molecular_system.atoms_Name[i], scale=(.9*atom_scaling , .9*atom_scaling , .9*atom_scaling ))
        
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



        X, Y, Z = self.orbital_generator.grid.return_grid_arrays()

        D_AB_contur = self.mlab.contour3d(X, Y, Z, (self.D_AB) ,contours=contours,opacity=0.5) 
    
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

    def Plot_D_AB_log(self, Plot_Atoms = 1, Atom_Names = 1, Plot_Bonds = 1 , Atom_Scaling = 1.0, Bond_Scaling = 1.0, contours = 6 ):  
        
        #from mayavi import mlab
        
        self.mlab.figure("Dispersion ", bgcolor=(.5, .5, .75), size=(1000, 1000))
        self.mlab.clf()
        
        if Plot_Atoms:
            for i in xrange(self.molecular_system.nAtoms):
                self.mlab.points3d(self.molecular_system.atoms_R[i,0], self.molecular_system.atoms_R[i,1], self.molecular_system.atoms_R[i,2], 
                    scale_factor=self.molecular_system.Atoms_Scale[letters(self.molecular_system.atoms_Name[i])] * Atom_Scaling  , resolution=20, 
                    color=self.molecular_system.Atoms_Color[letters(self.molecular_system.atoms_Name[i])], scale_mode='none')
            
        if Atom_Names:
            for i in xrange(self.molecular_system.nAtoms):
                self.mlab.text3d(self.molecular_system.atoms_R[i,0], self.molecular_system.atoms_R[i,1], self.molecular_system.atoms_R[i,2], 
                            self.molecular_system.atoms_Name[i],scale=(.2, .2, .2))
        
        if Plot_Bonds:        
            for i in xrange(len(self.molecular_system.bonds)):
                self.mlab.plot3d(
                    self.np.array([self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][0]), 0] , 
                                   self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][1]), 0]]),
                    self.np.array([self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][0]), 1] , 
                                   self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][1]), 1]]), 
                    self.np.array([self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][0]), 2] , 
                                   self.molecular_system.atoms_R[ self.molecular_system.atoms_Name.index(self.molecular_system.bonds[i][1]), 2]]),
                    tube_radius=0.2 * Bond_Scaling , color=(233.0/255, 165.0/255, 165.0/255))
        
        X, Y, Z = self.orbital_generator.grid.return_grid_arrays()

        self.mlab.contour3d(X, Y, Z, self.np.log(self.D_AB), contours = contours, opacity=0.5, vmax = 4.8312953302547844e-05) 
    
        self.mlab.show()
    

    

    def get_dispersion_index(self, R_max_multip = 3.0, x_n= 50, y_n= 50, z_n= 50, monomer_A=None, monomer_B=None):

        self.get_geminals( R_max_multip = R_max_multip, x_n= x_n, y_n= y_n, z_n= z_n)

        self.get_NumberOfFragments()
        self.get_Fragments()
        self.get_Epsilons()
        Monomers_set = set( self.Monomers )

        if monomer_A is None and monomer_B is None: 
            self.Calc_D( Monomer_A= Monomers_set.pop(), Monomer_B=Monomers_set.pop() )
        else:
            self.Calc_D( Monomer_A= monomer_A, Monomer_B=monomer_B )