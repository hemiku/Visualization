class Geminal():
    
    grid = None
    values = None

    def __init__(self, grid, values):
        self.grid = grid
        self.values = values

class GeminalGenerator():

    nGeminal = None
    Orb2Gem = None
    G_coeff = None
    geminals = None

    def __init__(self):

        pass

    def Calc_Geminals(self):

        self.geminals = self.np.zeros( [self.nGeminal, self.np.shape(self.X)[0], self.np.shape(self.Y)[1], self.np.shape(self.Z)[2]], dtype=self.np.float64)

        for i in range(len(self.Orb2Gem)):
            self.GEMINALS_init[self.Orb2Gem[i], :, :, :] += ( self.G_coeff[i] * self.MO[i + self.inactive, :, :, :] ) ** 2
