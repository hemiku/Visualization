
#from visualisation.grid import Grid


class AOParameters():
    
    import numpy as np
    
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


class OrbitalsGenerator( ):

    import numpy as np
    import visualization.grid

    grid = None

    nAtoms: int = None
    atoms_R = None

    spherical = None
    nb = None
    coeff = None
    basis = None
    basis_norm = None

    AOs = None
    MOs = None


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

        if grid is None:
            X, Y, Z = self.grid.return_grid_arrays()
        else:
            X, Y, Z = grid.return_grid_arrays()

        if self.spherical:
            dx = X[1, 0, 0] - X[0, 0, 0]
            dy = Y[0, 1, 0] - Y[0, 0, 0]
            dz = Z[0, 0, 1] - Z[0, 0, 0]

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

                if (j == 0):

                    for k in range(self.np.shape(self.basis[n][j])[1] - 1):

                        for l in range(self.np.shape(self.basis[n][j])[0]):
                            # if (self.basis[n][j])[l, k+1]  > 10**-16:
                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR)  # / ( (self.basis[n][j])[l,0]  **( 3.0/2.0))

                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                        nOrb += 1

                if (j == 1):

                    for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                        for m in ['x', 'y', 'z']:
                            if (m == 'x'):

                                for l in range(self.np.shape(self.basis[n][j])[0]):
                                    AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x)
                                AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                nOrb += 1

                            if (m == 'y'):

                                for l in range(self.np.shape(self.basis[n][j])[0]):
                                    AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y)
                                AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                nOrb += 1

                            if (m == 'z'):

                                for l in range(self.np.shape(self.basis[n][j])[0]):
                                    AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (z)
                                AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                nOrb += 1
                if (j == 2):
                    if (self.spherical == 1):

                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in range(-j, j + 1, 1):

                                    if (m == -2):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (AOParameters.sqrt3 * x * y)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == -1):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (AOParameters.sqrt3 * y * z)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 0):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (-0.5 * (x * x + y * y) + z * z)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 1):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (AOParameters.sqrt3 * x * z)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 2):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (
                                                                             AOParameters.sqrt3p4 * (x * x - y * y))
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1
                    else:
                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in ['xx', 'yy', 'zz', 'xy', 'xz', 'yz']:

                                    if (m == 'xx'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x ** 2)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'yy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y ** 2)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'zz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (z ** 2)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * y * AOParameters.sqrt3)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -( self.basis[n][j])[l, 0] * RR) * (x * z * AOParameters.sqrt3)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'yz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y * z * AOParameters.sqrt3)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                if (j == 3):
                    if (self.spherical == 1):

                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in range(-j, j + 1, 1):

                                    if (m == -3):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (AOParameters.sqrt10p16 * ( 3 * x ** 2 * y - y ** 3))
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == -2):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( AOParameters.sqrt15 * (x * y * z))
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == -1):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( AOParameters.sqrt125p108_1440 * AOParameters.sqrt45p100 * ( -x ** 2 * y - y ** 3 + AOParameters.sqrt40_25 * y * z ** 2))
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 0):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( z ** 3 - AOParameters.sqrt45p100 * ( AOParameters.sqrt5 * x ** 2 * z + AOParameters.sqrt5 * y ** 2 * z))
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 1):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( AOParameters.sqrt125p108_1440 * AOParameters.sqrt45p100 * ( -x ** 3 - x * y ** 2 + AOParameters.sqrt40_25 * x * z ** 2))
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 2):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( AOParameters.sqrt30p25 * AOParameters.sqrt50p16 * (x ** 2 * z - y ** 2 * z))
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 3):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * ( AOParameters.sqrt10p16 * ( x ** 3 - 3 * x * y ** 2))
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                    else:
                        if (len(self.np.shape(self.basis[n][j])) == 2):
                            for k in range(self.np.shape(self.basis[n][j])[1] - 1):
                                for m in ['xxx', 'yyy', 'zzz', 'xyy', 'xxy', 'xxz', 'xzz', 'yzz', 'yyz', 'xyz']:

                                    if (m == 'xxx'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * x * x)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'yyy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y * y * y)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'zzz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (z * z * z)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xyy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * y * y * AOParameters.sqrt5)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xxy'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * x * y * AOParameters.sqrt5)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xxz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * x * z * AOParameters.sqrt5)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xzz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * z * z * AOParameters.sqrt5)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'yzz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y * z * z * AOParameters.sqrt5)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'yyz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (y * y * z * AOParameters.sqrt5)
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1

                                    if (m == 'xyz'):
                                        for l in range(self.np.shape(self.basis[n][j])[0]):
                                            AO[nOrb, :, :, :] += (self.basis[n][j])[l, k + 1] * self.np.exp( -(self.basis[n][j])[l, 0] * RR) * (x * y * z * AOParameters.sqrt15 )
                                        AO[nOrb, :, :, :] = (self.basis_norm[n][j])[k] * AO[nOrb, :, :, :]
                                        nOrb += 1
