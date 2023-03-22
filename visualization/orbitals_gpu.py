
#from visualisation.grid import Grid


class AOParameters_gpu():
    
    import numpy as np
    
    universal_norm = ( 2 / np.pi)**(3/4)
    sqrt3 = np.sqrt(3.0)
    sqrt3p4 = np.sqrt(3.0 / 4.0)
    sqrt5 = np.sqrt(5.0)
    sqrt15 = np.sqrt(15.0)
    sqrt10p16 = np.sqrt(10.0 / 16.0)
    sqrt45p100 = np.sqrt(45.0 / 100)
    sqrt40_25 = np.sqrt(40.0 / 25)
    sqrt125p108_1440 = np.sqrt(125.0 / (108 - np.sqrt(1440)))
    sqrt30p25 = np.sqrt(30.0 / 25.0)
    sqrt50p16 = np.sqrt(50.0 / 16.0)


class OrbitalsGenerator_gpu( ):

    import numpy as np
    import visualization.grid

    grid:visualization.grid.Grid = None

    nAtoms: int = None
    atoms_R:np.ndarray = None

    spherical = None
    nb:int = None
    coeff:np.ndarray = None
    basis:list = None
    basis_norm:list = None

    AOs:np.ndarray = None
    MOs:np.ndarray = None


    def __init__(self, nAtoms = None, atoms_R = None, spherical=None, nb = None, coeff = None, basis = None, basis_norm = None, grid = None ):

        self.nAtoms = nAtoms
        self.atoms_R = atoms_R    

        self.spherical=spherical
        self.nb = nb
        self.coeff = coeff
        self.basis = basis
        self.basis_norm = basis_norm

        if grid is None:
    
            self.grid = self.visualization.grid.Grid()

        else:

            self.grid = grid


    def set_Atoms_R( self, Atoms_R):

        self.Atoms_R = Atoms_R

    def init_grid(self):

        self.grid.generate_grid_boundaries( atoms_R= self.atoms_R, basis= self.basis )

    def init_AOs(self):

        X, Y, Z = self.grid.return_grid_arrays()

        self.AOs = self.np.zeros([self.nb, self.np.shape(X)[0], self.np.shape(Y)[1], self.np.shape(Z)[2]], dtype=self.np.float64)

    def init_MOs(self):

        X, Y, Z = self.grid.return_grid_arrays()

        self.MOs = self.np.zeros([self.nb, self.np.shape(X)[0], self.np.shape(Y)[1], self.np.shape(Z)[2]], dtype=self.np.float64)

    def calc_MOs(self, start = None, stop = None):

        X, Y, Z = self.grid.return_grid_arrays()

        self.MOs = self.np.zeros( [self.nb, self.np.shape(X)[0], self.np.shape(Y)[1], self.np.shape(Z)[2]], dtype=self.np.float64)

#        for i in range(self.inactive + 2 * self.nGeminal):
        for i in range(self.nb):
            for j in range(self.nb):
                self.MOs[i, :, :, :] += self.coeff[i, j] * self.AOs[j, :, :, :]

    def calc_AOs(self, AO, grid=None):

        # from scipy.special import sph_harm  
        import cupy as cp      

        if grid is None:
            grid = self.grid

        X, Y, Z = self.grid.return_grid_arrays()

        dx = ( grid.x_max - grid.x_min )/ ( grid.x_n -1 )
        dy = ( grid.y_max - grid.y_min )/ ( grid.y_n -1 )
        dz = ( grid.z_max - grid.z_min )/ ( grid.z_n -1 )

        dv = dx * dy * dz

        nOrb = 0

        for n in range( self.nAtoms ):

            x = X - self.atoms_R[n, 0]
            y = Y - self.atoms_R[n, 1]
            z = Z - self.atoms_R[n, 2]

            RR = x * x + y * y + z * z

            r = self.np.sqrt(RR)
            theta = self.np.arctan(y / x)
            phi = self.np.arccos(z / r)

            
            for j in range(len(self.basis[n])):

                #if nOrb > 20:
                #    break

                if (j == 0):

                    for k in range(self.np.shape(self.basis[n][j])[1] - 1):

                        for l in range(self.np.shape(self.basis[n][j])[0]):
                    
                            alpha = (self.basis[n][j])[l, 0]
                            coefficient = (self.basis[n][j])[l, k + 1]

                            norm = 2**(3/4)*alpha**(3/4)/self.np.pi**(3/4)

                            AO[nOrb, :, :, :] +=   norm * coefficient * self.np.exp( -alpha * RR ) 
                        
                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] * (self.basis_norm[n][j])[k]
                        nOrb += 1

                if (j == 1):

                    for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                        for m in ['x', 'y', 'z']:

                            if (m == 'x'):
                                
                                for l in range(self.np.shape(self.basis[n][j])[0]):
                                    
                                    alpha = (self.basis[n][j])[l, 0]
                                    coefficient = (self.basis[n][j])[l, k + 1]

                                    norm =  2*2**(3/4)*alpha**(5/4)/self.np.pi**(3/4)

                                    AO[nOrb, :, :, :] += norm * coefficient * self.np.exp( -alpha * RR) * (x)

                                #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] * (self.basis_norm[n][j])[k]

                            if (m == 'y'):
                                
                                for l in range(self.np.shape(self.basis[n][j])[0]):

                                    alpha = (self.basis[n][j])[l, 0]
                                    coefficient = (self.basis[n][j])[l, k + 1]

                                    norm = 2*2**(3/4)*alpha**(5/4)/self.np.pi**(3/4)
                                    
                                    AO[nOrb, :, :, :] += norm * coefficient  * self.np.exp( -alpha * RR) * (y) 
                                #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / (self.basis_norm[n][j])[k]


                            if (m == 'z'):
                                
                                for l in range(self.np.shape(self.basis[n][j])[0]):
                                    alpha = (self.basis[n][j])[l, 0]
                                    coefficient = (self.basis[n][j])[l, k + 1]
                                    norm = 2*2**(3/4)*alpha**(5/4)/self.np.pi**(3/4)
                                    #norm = 1 
                                    AO[nOrb, :, :, :] += norm * coefficient * self.np.exp( -alpha * RR) * (z)
                                #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]

                            #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] * (self.basis_norm[n][j])[k]   
                            nOrb += 1

                if (j == 2):
                    if (self.spherical == 1):

                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in self.np.arange(-j,j+1):
                                    if (m == -2):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            
                                            alpha = (self.basis[n][j])[l, 0]
                                            coefficient = (self.basis[n][j])[l, k + 1]

                                            norm = 1

                                            AO[nOrb, :, :, :] += norm  * coefficient *  ( x * y ) *self.np.exp(-alpha*RR)   

                                    if (m == -1):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            
                                            alpha = (self.basis[n][j])[l, 0]
                                            coefficient = (self.basis[n][j])[l, k + 1]

                                            norm = 1

                                            AO[nOrb, :, :, :] += norm * coefficient * ( y * z ) *self.np.exp(-alpha * RR )   
                                            
                                    if (m == 0):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):

                                            alpha = (self.basis[n][j])[l, 0]
                                            coefficient = (self.basis[n][j])[l, k + 1]

                                            norm = 1 / self.np.sqrt(12)

                                            AO[nOrb, :, :, :] += norm * coefficient * (-x**2 - y**2 + 2*z**2)*self.np.exp( -alpha * RR )   
                                    
                                    if (m == 1):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            
                                            alpha = (self.basis[n][j])[l, 0]
                                            coefficient = (self.basis[n][j])[l, k + 1]

                                            norm = 1

                                            AO[nOrb, :, :, :] += norm  * coefficient * ( x * z ) *self.np.exp(-alpha * RR )   
                                            
                                    if (m == 2):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            
                                            alpha = (self.basis[n][j])[l, 0]
                                            coefficient = (self.basis[n][j])[l, k + 1]

                                            norm = 0.5
                                            
                                            AO[nOrb, :, :, :] += norm * coefficient * (x**2 - y**2) *self.np.exp( -alpha * RR )   
                                    
                                    AO[nOrb, :, :, :] = AO[nOrb, :, :, :] * (self.basis_norm[n][j])[k]
                                    nOrb += 1


                    else:
                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']:

                                    if (m == 'xx'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * ( (self.basis[n][j])[l, 0])**((1+2*(j+1))/4) * (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x ** 2)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv 
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'yy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * ( (self.basis[n][j])[l, 0])**((1+2*(j+1))/4) * (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y ** 2)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv 
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'zz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * ( (self.basis[n][j])[l, 0])**((1+2*(j+1))/4) * (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (z ** 2)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'xy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * ( (self.basis[n][j])[l, 0])**((1+2*(j+1))/4) * (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * y * AOParameters_gpu.sqrt3)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv 
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'xz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * ( (self.basis[n][j])[l, 0])**((1+2*(j+1))/4) * (self.basis[n][j])[l, k + 1] * self.np.exp( -( self.basis[n][j])[l, 0] * RR) * (x * z * AOParameters_gpu.sqrt3)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv 
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'yz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * ( (self.basis[n][j])[l, 0])**((1+2*(j+1))/4) * (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y * z * AOParameters_gpu.sqrt3)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv 
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1 

                if (j == 3):
                    if (self.spherical == 1):

                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in range(-j, j + 1, 1):

                                    if (m == -3):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (AOParameters_gpu.sqrt10p16 * ( 3 * x ** 2 * y - y ** 3))
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == -2):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( AOParameters_gpu.sqrt15 * (x * y * z))
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == -1):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( AOParameters_gpu.sqrt125p108_1440 * AOParameters_gpu.sqrt45p100 * ( -x ** 2 * y - y ** 3 + AOParameters_gpu.sqrt40_25 * y * z ** 2))
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 0):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( z ** 3 - AOParameters_gpu.sqrt45p100 * ( AOParameters_gpu.sqrt5 * x ** 2 * z + AOParameters_gpu.sqrt5 * y ** 2 * z))
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 1):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( AOParameters_gpu.sqrt125p108_1440 * AOParameters_gpu.sqrt45p100 * ( -x ** 3 - x * y ** 2 + AOParameters_gpu.sqrt40_25 * x * z ** 2))
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 2):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( AOParameters_gpu.sqrt30p25 * AOParameters_gpu.sqrt50p16 * (x ** 2 * z - y ** 2 * z))
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 3):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += AOParameters_gpu.universal_norm * (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( AOParameters_gpu.sqrt10p16 * ( x ** 3 - 3 * x * y ** 2))
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                    else:
                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in ['xxx', 'yyy', 'zzz', 'xyy', 'xxy', 'xxz', 'xzz', 'yzz', 'yyz', 'xyz']:

                                    if (m == 'xxx'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * x * x)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'yyy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y * y * y)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'zzz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (z * z * z)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'xyy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * y * y * AOParameters_gpu.sqrt5)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'xxy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * x * y * AOParameters_gpu.sqrt5)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'xxz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * x * z * AOParameters_gpu.sqrt5)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'xzz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * z * z * AOParameters_gpu.sqrt5)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'yzz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y * z * z * AOParameters_gpu.sqrt5)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'yyz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y * y * z * AOParameters_gpu.sqrt5)
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1

                                    if (m == 'xyz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * y * z * AOParameters_gpu.sqrt15 )
                                        #AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :] * sqrt_dv
                                        #AO[nOrb, :, :, :] = AO[nOrb, :, :, :] / self.np.sqrt( self.np.sum(AO[nOrb, :, :, :] ** 2 ) )
                                        nOrb += 1
        
        #for n in range( self.nb ):
        #    AO[n, :, :, :] = AO[n, :, :, :] / self.np.sqrt( self.np.sum(AO[n, :, :, :] ** 2 * dv ) )

