class Geminal():
    
    grid = None
    values = None

    def __init__(self, grid, values):
        self.grid = grid
        self.values = values

import numpy as np


class GeminalGenerator():

    MOs:np.ndarray
    inactive: int 
    n_geminal: int
    Orb2Gem:np.ndarray
    G_coeff:np.ndarray
    geminals:np.ndarray

    def __init__(self,  MOs:np.ndarray, 
                        n_geminal:int, 
                        inactive:int, 
                        Orb2Gem:np.ndarray, 
                        G_coeff:np.ndarray, 
                        geminals:np.ndarray):

        self.inactive = inactive
        self.MOs = MOs
        self.nGeminal = n_geminal
        self.Orb2Gem = Orb2Gem
        self.G_coeff = G_coeff
        self.geminals = geminals


    def calc_geminals(self):

        geminals_shape = np.array( self.MOs.shape )
        geminals_shape[0] = self.nGeminal

        self.geminals = np.zeros(geminals_shape, dtype=np.float64)

        for i in range( len( self.Orb2Gem )):
            self.geminals[self.Orb2Gem[i], :, :, :] += ( self.G_coeff[i] * self.MOs[i + self.inactive, :, :, :] ) ** 2
