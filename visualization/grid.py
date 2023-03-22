
class Grid():

    import numpy as np

    initail = None

    R_max_multip = None

    x_min:np.float64 
    x_max:np.float64 
    x_n:np.int32
    y_min:np.float64 
    y_max:np.float64 
    y_n:np.int32 
    z_min:np.float64 
    z_max:np.float64 
    z_n:np.int32


    def __init__(self,  R_max_multip:float =12.0, 
                        x_min:np.float32=None, 
                        x_max:np.float32=None, 
                        x_n:np.int32= 100, 
                        y_min:np.float32=None, 
                        y_max:np.float32=None, 
                        y_n:np.int32=100,
                        z_min:np.float32=None, 
                        z_max:np.float32=None, 
                        z_n:np.int32=100):

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

    def return_grid_arrays(self) -> np.ndarray:

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

    # def generate_grid_boundaries_obsolete(self, basis, Atoms_R):

    #     Basis_max = []
    #     for i in range(len(basis)):
    #         for j in range(len(basis[i])):
    #             Basis_max.append(1.64526336574595 / self.np.sqrt(( basis[i][j] ).max()))

    #     R_max = self.R_max_multip * self.np.max(Basis_max)

    #     if self.x_max is None:
    #         self.x_max = self.np.max( Atoms_R[:, 0]) + R_max
    #     if self.x_min is None:
    #         self.x_min = self.np.min( Atoms_R[:, 0]) - R_max

    #     if self.y_max is None:
    #         self.y_max = self.np.max( Atoms_R[:, 1]) + R_max
    #     if self.y_min is None:
    #         self.y_min = self.np.min( Atoms_R[:, 1]) - R_max

    #     if self.z_max is None:
    #         self.z_max = self.np.max( Atoms_R[:, 2]) + R_max
    #     if self.z_min is None:
    #         self.z_min = self.np.min( Atoms_R[:, 2]) - R_max
