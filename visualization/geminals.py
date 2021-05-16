class Geminal():
    
    grid = None
    values = None

    def __init__(self, grid, values):
        self.grid = grid
        self.values = values

class GeminalGenerator():

    import numpy as np

    MOs = None
    inactive = None
    n_geminal = None
    Orb2Gem = None
    G_coeff = None
    geminals = None

    def __init__(self, MOs = MOs, n_geminal = n_geminal, inactive = inactive, Orb2Gem = Orb2Gem, G_coeff = G_coeff, geminals = geminals ):

        self.inactive = inactive
        self.MOs = MOs
        self.nGeminal = n_geminal
        self.Orb2Gem = Orb2Gem
        self.G_coeff = G_coeff
        self.geminals = geminals


    def Calc_Geminals(self):

        geminals_shape = self.np.array( self.MOs.shape )
        geminals_shape[0] = self.nGeminal

        self.geminals = self.np.zeros(geminals_shape, dtype=self.np.float64)

        for i in range( len( self.Orb2Gem )):
            self.geminals[self.Orb2Gem[i], :, :, :] += ( self.G_coeff[i] * self.MOs[i + self.inactive, :, :, :] ) ** 2
