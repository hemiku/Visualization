
class Grid():

    import numpy as np

    initail = None

    R_max_multip = None

    x_min = None
    x_max = None
    x_n = None
    y_min = None
    y_max = None
    y_n = None
    z_min = None
    z_max = None
    z_n = None


    def __init__(self, R_max_multip=12.0, x_min=None, x_max=None, x_n=None, y_min=None, y_max=None, y_n=None, z_min=None, z_max=None, z_n=None):

        self.R_max_multip = R_max_multip

        self.x_min = x_min
        self.x_max = x_max
        self.x_n = x_n
        self.y_min = y_min
        self.y_max = y_max
        self.y_n = y_n
        self.z_min = z_min
        self.z_max = z_max
        self.z_n = z_n

    def return_grid_arrays(self):

        return self.np.mgrid[self.x_min:self.x_max:self.x_n*1j, self.y_min:self.y_max:self.y_n*1j, self.z_min:self.z_max:self.z_n*1j]

    def generate_grid_boundaries(self, basis, atoms_R):

        Basis_max = []
        for i in range(len(basis)):
            for j in range(len(basis[i])):
                Basis_max.append(1.64526336574595 / self.np.sqrt(( basis[i][j] ).max()))

        R_max = self.R_max_multip * self.np.max(Basis_max)


        self.x_max = self.np.max( atoms_R[:, 0]) + R_max
        self.x_min = self.np.min( atoms_R[:, 0]) - R_max

        self.y_max = self.np.max( atoms_R[:, 1]) + R_max
        self.y_min = self.np.min( atoms_R[:, 1]) - R_max

        self.z_max = self.np.max( atoms_R[:, 2]) + R_max
        self.z_min = self.np.min( atoms_R[:, 2]) - R_max

    def generate_grid_boundaries_obsolete(self, basis, Atoms_R):

        Basis_max = []
        for i in range(len(basis)):
            for j in range(len(basis[i])):
                Basis_max.append(1.64526336574595 / self.np.sqrt(( basis[i][j] ).max()))

        R_max = self.R_max_multip * self.np.max(Basis_max)

        if self.x_max is None:
            self.x_max = self.np.max( Atoms_R[:, 0]) + R_max
        if self.x_min is None:
            self.x_min = self.np.min( Atoms_R[:, 0]) - R_max

        if self.y_max is None:
            self.y_max = self.np.max( Atoms_R[:, 1]) + R_max
        if self.y_min is None:
            self.y_min = self.np.min( Atoms_R[:, 1]) - R_max

        if self.z_max is None:
            self.z_max = self.np.max( Atoms_R[:, 2]) + R_max
        if self.z_min is None:
            self.z_min = self.np.min( Atoms_R[:, 2]) - R_max
