
class DaltonPlt:
    import numpy as np

    def __init__(self, path):
        self.path = path

        self.Atoms_R = None
        self.Atoms_Charge = None
        self.Atoms_Name = None

        self.basis = None
        self.basis_norm = None
        self.basis_norm2 = None

        self.X = None
        self.Y = None
        self.Z = None

        with open(path, 'r') as f:
            self.Out = f.read()


    def read_data(self):
        self.record1 = self.np.fromstring(self.Out[:4], dtype=self.np.int32)[0]
        self.record2 = self.np.fromstring(self.Out[4:8], dtype=self.np.int32)[0]

        self.x_n = self.np.fromstring(self.Out[8:12], dtype=self.np.int32)[0]
        self.y_n = self.np.fromstring(self.Out[12:16], dtype=self.np.int32)[0]
        self.z_n = self.np.fromstring(self.Out[16:120], dtype=self.np.int32)[0]

        self.z_min = self.np.fromstring(self.Out[20:24], dtype=self.np.float32)[0]
        self.z_max = self.np.fromstring(self.Out[24:28], dtype=self.np.float32)[0]

        self.y_min = self.np.fromstring(self.Out[28:32], dtype=self.np.float32)[0]
        self.y_max = self.np.fromstring(self.Out[32:36], dtype=self.np.float32)[0]

        self.x_min = self.np.fromstring(self.Out[36:40], dtype=self.np.float32)[0]
        self.x_max = self.np.fromstring(self.Out[40:44], dtype=self.np.float32)[0]

        self.values = self.np.reshape(self.np.fromstring(self.Out[44:], dtype=self.np.float32, count=self.x_n * self.y_n * self.z_n), [self.x_n, self.y_n, self.z_n])

        self.XYZ_str = "  self.np.mgrid[ " + str(self.x_min) + ":" + str(self.x_max) + ":" + str(self.x_n) + "j, " + \
                       str(self.y_min) + ":" + str(self.y_max) + ":" + str(self.y_n) + "j, " + \
                       str(self.z_min) + ":" + str(self.z_max) + ":" + str(self.z_n) + "j ]"

    def get_grid(self):
        exec("X, Y ,Z = " + self.XYZ_str)

        return X, Y, Z

    def gen_grid(self):
        exec("self.X, self.Y ,self.Z = " + self.XYZ_str)

    def get_self_grid(self):
        if not ('X' in dir(self) and 'Y' in dir(self) and 'Z' in dir(self)):
            self.gen_grid()

        return self.X, self.Y, self.Z
