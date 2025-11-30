from typing import Dict, List, Optional

import numpy as np

from visualization.backends import get_backend, VisualizationBackend
from visualization.geminals import GeminalGenerator
from visualization.inputs import Input
from visualization.orbitals import OrbitalsGenerator
from visualization.molecular_system import MolecularSystem
from visualization.constants import (
    BACKGROUND_COLORS,
    ATOM_COLORS,
    ATOM_SCALES,
    ATOM_MASSES,
    ELEMENT_TO_ATOMIC_NUMBER,
)
import visualization.utils as u


class VisualizationData:
    """Data class holding visualization configuration."""

    background_colors = BACKGROUND_COLORS
    Atoms_Color = ATOM_COLORS
    Atoms_Scale = ATOM_SCALES
    Atoms_Mass = ATOM_MASSES
    Name_to_Atom_Number = ELEMENT_TO_ATOMIC_NUMBER


class Visualization:
    """Main class for molecular orbital visualization.

    This class orchestrates the complete visualization pipeline:
    1. Parse input files (Dalton, Molpro, Molden)
    2. Generate 3D grids
    3. Calculate atomic and molecular orbitals
    4. Render visualizations

    Attributes
    ----------
    input_type : str
            Type of input file ('Dalton', 'Molpro', 'MolproSapt')
    input_name : str
            Base name of input file (without extension)
    data_input : Input
            Parsed input file data
    molecular_system : MolecularSystem
            Container for molecular structure and orbital data
    orbital_generator : OrbitalsGenerator
            Calculator for AOs and MOs

    Examples
    --------
    >>> vis = Visualization(input_type='Dalton', input_name='water')
    >>> vis.get_orbitals(x_n=50, y_n=50, z_n=50)
    >>> vis.plot_Orbitals(orbital_numbers=[0, 1, 2])
    """

    def __init__(
        self,
        input_type: Optional[str] = None,
        input_sub_type: str = "Output",
        input_name: Optional[str] = None,
        file_string: Optional[str] = None,
        BAS_filename: Optional[str] = None,
        data_source: Optional[str] = None,
        backend: Optional[str] = None,
    ) -> None:
        """Initialize the visualization system.

        Parameters
        ----------
        input_type : str, optional
                Input file type ('Dalton', 'Molpro', 'MolproSapt')
        input_sub_type : str, default 'Output'
                Sub-type for input parsing
        input_name : str, optional
                Base name of input file
        file_string : str, optional
                Raw file content (alternative to file)
        BAS_filename : str, optional
                Basis set file name
        data_source : str, optional
                Data source identifier
        backend : str, optional
                Visualization backend ('pyvista' or 'mayavi')
        """
        self.input_type = input_type
        self.input_sub_type = input_sub_type
        self.input_name = input_name
        self.file_string = file_string
        self.BAS_filename = BAS_filename
        self.data_source = data_source

        self.data_input: Optional[Input] = None
        self.molecular_system: Optional[MolecularSystem] = None
        self.orbital_generator: Optional[OrbitalsGenerator] = None
        self.visualization_data = VisualizationData()

        # Initialize visualization backend (default: pyvista)
        self._backend: VisualizationBackend = get_backend(backend)

        # Initialize data input if both type and name are provided
        if self.input_type is not None and self.input_name is not None:
            self.initialize_data_input()

        self.initialize_molecular_system()

    def set_input(
        self,
        input_type=None,
        input_sub_type=None,
        input_name=None,
        file_string=None,
        BAS_filename=None,
        data_source=None,
    ):

        import visualization.inputs
        import visualization.input_molpro

        INPUT_TYPES: dict[str, visualization.inputs.Input]

        INPUT_TYPES = {
            "Dalton": visualization.inputs.DaltonInput,
            "Molpro": visualization.input_molpro.MolproInput,
            "MolproSapt": visualization.input_molpro.MolproSaptInput,
        }

        try:
            return INPUT_TYPES[input_type](
                input_type=input_sub_type,
                input_name=input_name,
                file_string=file_string,
                BAS_filename=BAS_filename,
                data_source=data_source,
            )
        except KeyError:
            raise Exception(f"Input_type: {input_type} not found")

    def set_input_type(self, input_type):

        self.input_type = input_type

    def set_input_sub_type(self, input_sub_type):

        self.input_sub_type = input_sub_type

    def set_input_name(self, input_name):

        self.input_name = input_name

    def set_source(self, source):

        self.source = source

    def set_data_source(self, data_source):

        self.data_source = data_source

    def initialize_data_input(self):

        # import visualization.inputs

        self.data_input = self.set_input(
            input_type=self.input_type,
            input_sub_type=self.input_sub_type,
            input_name=self.input_name,
        )

    def initialize_molecular_system(self):

        import visualization.molecular_system

        self.molecular_system = visualization.molecular_system.MolecularSystem()

    def get_geometry(self, get_bonds=True):

        if self.data_input is None or self.molecular_system is None:

            print("Either self.data_input or self.MolecularSystem")

        else:

            nAtoms = self.data_input.get_nAtoms()
            atoms_R, atoms_Charge, atoms_Name = self.data_input.get_atoms()

            self.molecular_system.set_nAtoms(nAtoms=nAtoms)
            self.molecular_system.set_atoms_R(atoms_R=atoms_R)
            self.molecular_system.set_atoms_Charge(atoms_Charge=atoms_Charge)
            self.molecular_system.set_atoms_Name(atoms_Name=atoms_Name)

            if get_bonds:
                bonds = self.data_input.get_bonds()
                self.molecular_system.set_bonds(bonds=bonds)

    def get_orbital_data(self, get_bonds=True):

        if self.data_input is None or self.molecular_system is None:

            print("Either self.data_input or self.MolecularSystem")

        else:

            from visualization.orbitals import OrbitalsGenerator

            self.get_geometry(get_bonds=get_bonds)

            self.molecular_system.set_spherical(spherical=self.data_input.get_spherical())

            self.molecular_system.set_nb(nb=self.data_input.get_nb())

            self.molecular_system.set_coeff(coeff=self.data_input.get_coeff())

            self.molecular_system.set_basis_and_norms(input=self.data_input.get_basis())

            self.orbital_generator = OrbitalsGenerator(
                nAtoms=self.molecular_system.get_nAtoms(),
                atoms_R=self.molecular_system.get_atoms_R(),
                spherical=self.molecular_system.get_spherical(),
                nb=self.molecular_system.get_nb(),
                coeff=self.molecular_system.get_coeff(),
                basis=self.molecular_system.get_basis(),
                basis_norm=self.molecular_system.get_basis_norm(),
            )

    def generate_ao_orbitals(self):

        self.orbital_generator.calc_aos(AO=self.orbital_generator.AOs)
        self.molecular_system.AOs = self.orbital_generator.AOs

    def generate_ao_orbitals_gpu(self):

        self.orbital_generator.calc_aos_gpu(AO=self.orbital_generator.AOs)
        self.molecular_system.AOs = self.orbital_generator.AOs

    def generate_ao_orbitals_numba(self):

        self.orbital_generator.calc_aos(AO=self.orbital_generator.AOs)
        self.molecular_system.AOs = self.orbital_generator.AOs

    def generate_mo_orbitals(self, gpu=False, mos_limit=None):

        if mos_limit is None:
            self.orbital_generator.calc_mos()
        else:
            self.orbital_generator.calc_mos(mos_limit=mos_limit)

        self.molecular_system.MOs = self.orbital_generator.MOs

    def generate_mo_orbitals_gpu(self, mos_limit=None):

        if mos_limit is None:
            self.orbital_generator.calc_mos_gpu()
        else:
            self.orbital_generator.calc_mos_gpu(mos_limit=mos_limit)
        self.molecular_system.MOs = self.orbital_generator.MOs

    def generate_mo_orbitals_gpu_fast(self, mos_limit=None):

        if mos_limit is None:
            self.orbital_generator.calc_mos_gpu_fast()
        else:
            self.orbital_generator.calc_mos_gpu_fast(mos_limit=mos_limit)

        self.molecular_system.MOs = self.orbital_generator.MOs

    def generate_mo_orbitals_gpu_low_memory(self, mos_limit=None):

        self.orbital_generator.calc_mos_gpu_low_memory()
        self.molecular_system.MOs = self.orbital_generator.MOs

    def get_geminal_data(self):

        if self.data_input is None or self.molecular_system is None:

            print("Either self.data_input or self.MolecularSystem")

        else:

            self.data_input.get_inactive()
            self.data_input.get_electrons()
            self.data_input.get_occ()

            self.data_input.get_g_coeff()
            self.data_input.get_orb2gem()

            self.molecular_system.inactive = self.data_input.inactive
            self.molecular_system.G_coeff = self.data_input.G_coeff
            self.molecular_system.Orb2Gem = self.data_input.Orb2Gem
            self.molecular_system.n_geminals = self.data_input.nGeminal

            self.geminal_generator = GeminalGenerator(
                MOs=self.molecular_system.MOs,
                n_geminal=self.molecular_system.n_geminals,
                inactive=self.molecular_system.inactive,
                Orb2Gem=self.molecular_system.Orb2Gem,
                G_coeff=self.molecular_system.G_coeff,
                geminals=self.molecular_system.geminals,
            )

    def contour_process(self, contour, cube):

        if type(contour) == str:
            if contour[-1] == "%":

                contour_input_value = float(contour[:-1])

                cube_processed = np.sort(cube.flatten())[::-1]
                cube_csum = np.cumsum(cube_processed)
                cube_csum = cube_csum / cube_csum[-1]

                limit_value = cube_processed[cube_csum > ((contour_input_value / 100))][0]

                return [limit_value]

        if type(contour) == list:

            return_list = []

            for contour_element in contour:
                return_list.append(self.contour_process(contour_element, cube))

            return list(np.array(return_list).flat)

        else:
            return contour

    def get_generate_geminals(self):

        self.geminal_generator.calc_geminals()
        self.molecular_system.geminals = self.geminal_generator.geminals

    def plot_geometry(
        self,
        plot_atoms=True,
        atom_scaling=1.0,
        atom_names=True,
        atom_names_scaling=1.0,
        plot_bonds=False,
        bond_scaling=1.0,
        contours=6,
        background_color=None,
        scalarbar=False,
        auto_show=True,
        opacity=0.5,
        figure=None,
    ):

        _plot_type = "Geometry"
        _figure_title = _plot_type
        if background_color is None:
            _background_color = self.visualization_data.background_colors["White"]
        else:
            _background_color = background_color

        if figure is None:
            self._backend.create_figure(_figure_title, _background_color, (600, 600))

        if plot_atoms:
            self._plot_atoms(atom_scaling)

        if atom_names:
            self._atom_names(atom_names_scaling)

        if plot_bonds:
            self._plot_bonds(plot_bonds, bond_scaling)

        if auto_show:
            self._backend.show()

    def plot_orbitals(
        self,
        orbital_numbers,
        plot_atoms=True,
        atom_scaling=1.0,
        atom_names=True,
        atom_names_scaling=1.0,
        plot_bonds=True,
        bond_scaling=1.0,
        contours=6,
        background_color=None,
        scalarbar=False,
        auto_show=True,
        opacity=0.5,
        figure=None,
    ):

        _plot_type = "Orbitals_"
        _figure_title = _plot_type + str(orbital_numbers)

        if background_color is None:
            _background_color = self.visualization_data.background_colors["White"]
        else:
            _background_color = background_color

        if figure is None:
            self._backend.create_figure(_figure_title, _background_color, (600, 600))

        if plot_atoms:
            self._plot_atoms(atom_scaling)

        if atom_names:
            self._atom_names(atom_names_scaling)

        if plot_bonds:
            self._plot_bonds(plot_bonds, bond_scaling)

        for number in orbital_numbers:
            X, Y, Z = self.orbital_generator.grid.return_grid_arrays()
            self._backend.add_contour3d(
                X, Y, Z, self.molecular_system.MOs[number], contours=12, opacity=0.5
            )

        if auto_show:
            self._backend.show()

    def plot_orbitals_ao(
        self,
        orbital_numbers,
        plot_atoms=True,
        atom_scaling=1.0,
        atom_names=True,
        atom_names_scaling=1.0,
        plot_bonds=True,
        bond_scaling=1.0,
        contours=6,
        background_color=None,
        scalarbar=False,
        auto_show=True,
        opacity=0.5,
        figure=None,
    ):

        _plot_type = "Orbitals_AO_"
        _figure_title = _plot_type + str(orbital_numbers)

        if background_color is None:
            _background_color = self.visualization_data.background_colors["White"]
        else:
            _background_color = background_color

        if figure is None:
            self._backend.create_figure(_figure_title, _background_color, (600, 600))

        if plot_atoms:
            self._plot_atoms(atom_scaling)

        if atom_names:
            self._atom_names(atom_names_scaling)

        if plot_bonds:
            self._plot_bonds(plot_bonds, bond_scaling)

        for number in orbital_numbers:
            X, Y, Z = self.orbital_generator.grid.return_grid_arrays()
            self._backend.add_contour3d(
                X, Y, Z, self.molecular_system.AOs[number], contours=12, opacity=0.5
            )

        if auto_show:
            self._backend.show()

    def plot_geminals(
        self,
        geminal_numbers=[0],
        plot_atoms=True,
        atom_scaling=1.0,
        atom_names=True,
        atom_names_scaling=1.0,
        plot_bonds=False,
        bond_scaling=1.0,
        contours=6,
        background_color=None,
        scalarbar=False,
        auto_show=True,
        opacity=0.5,
        figure=None,
    ):

        _plot_type = "Geminals_"
        _figure_title = _plot_type + str(geminal_numbers)

        if background_color is None:
            _background_color = self.visualization_data.background_colors["White"]
        else:
            _background_color = background_color

        if figure is None:
            self._backend.create_figure(_figure_title, _background_color, (600, 600))

        if plot_atoms:
            self._plot_atoms(atom_scaling)

        if atom_names:
            self._atom_names(atom_names_scaling)

        if plot_bonds:
            self._plot_bonds(plot_bonds, bond_scaling)

        for number in geminal_numbers:
            X, Y, Z = self.orbital_generator.grid.return_grid_arrays()
            contour_values = self.contour_process(contours, self.molecular_system.geminals[number])
            self._backend.add_contour3d(
                X,
                Y,
                Z,
                self.molecular_system.geminals[number],
                contours=contour_values,
                opacity=opacity,
            )

        if auto_show:
            self._backend.show()

    def _plot_bonds(self, plot_bonds=True, bond_scaling=1.0):

        if plot_bonds and self.molecular_system.bonds is not None:
            for i, bond in enumerate(self.molecular_system.bonds):

                bond_begin = self.molecular_system.atoms_R[bond[0]]
                bond_end = self.molecular_system.atoms_R[bond[1]]

                bond_half = 0.5 * (
                    self.molecular_system.atoms_R[bond[0]] + self.molecular_system.atoms_R[bond[1]]
                )

                # First half of bond (from begin to midpoint)
                points_first = np.array([bond_begin, bond_half])
                self._backend.add_tube(
                    points_first,
                    radius=0.2 * bond_scaling,
                    color=self.visualization_data.Atoms_Color[
                        u.letters(self.molecular_system.atoms_Name[bond[0]])
                    ],
                )

                # Second half of bond (from end to midpoint)
                points_second = np.array([bond_end, bond_half])
                self._backend.add_tube(
                    points_second,
                    radius=0.2 * bond_scaling,
                    color=self.visualization_data.Atoms_Color[
                        u.letters(self.molecular_system.atoms_Name[bond[1]])
                    ],
                )

    def _plot_atoms(self, atom_scaling):

        for i in range(self.molecular_system.nAtoms):
            element = u.letters(self.molecular_system.atoms_Name[i])
            # Note: Mayavi's scale_factor is diameter, backend uses radius
            radius = 0.5 * atom_scaling * self.visualization_data.Atoms_Scale[element]
            self._backend.add_sphere(
                self.molecular_system.atoms_R[i, 0],
                self.molecular_system.atoms_R[i, 1],
                self.molecular_system.atoms_R[i, 2],
                radius=radius,
                color=self.visualization_data.Atoms_Color[element],
                resolution=20,
            )

    def _atom_names(self, atom_names_scaling):
        for i in range(self.molecular_system.nAtoms):
            element = u.letters(self.molecular_system.atoms_Name[i])
            self._backend.add_text3d(
                self.molecular_system.atoms_R[i, 0],
                self.molecular_system.atoms_R[i, 1],
                self.molecular_system.atoms_R[i, 2],
                self.molecular_system.atoms_Name[i],
                color=self.visualization_data.Atoms_Color[element],
                scale=0.9 * atom_names_scaling,
            )

    def plot_orbitals_MO(
        self,
        orbital_numbers=[0],
        plot_atoms=True,
        atom_scaling=1.0,
        atom_names=True,
        atom_names_scaling=1.0,
        plot_bonds=False,
        bond_scaling=1.0,
        contours=6,
        background_color=None,
        scalarbar=False,
        auto_show=True,
        figure=None,
    ):

        _plot_type = "Orbitals_MO_"
        _figure_title = _plot_type + str(orbital_numbers)

        if background_color is None:
            _background_color = self.visualization_data.background_colors["White"]
        else:
            _background_color = background_color

        if figure is None:
            self._backend.create_figure(_figure_title, _background_color, (600, 600))

        if plot_atoms:
            self._plot_atoms(atom_scaling)

        if atom_names:
            self._atom_names(atom_names_scaling)

        if plot_bonds:
            self._plot_bonds(plot_bonds, bond_scaling)

        for number in orbital_numbers:
            X, Y, Z = self.orbital_generator.grid.return_grid_arrays()
            contour_values = self.contour_process(contours, self.molecular_system.MOs[number])
            self._backend.add_contour3d(
                X, Y, Z, self.molecular_system.MOs[number], contours=contour_values, opacity=0.5
            )

        if auto_show:
            self._backend.show()

    def plot_orbitals_MO_square(
        self,
        orbital_numbers=[0],
        plot_atoms=True,
        atom_scaling=1.0,
        atom_names=True,
        atom_names_scaling=1.0,
        plot_bonds=True,
        bond_scaling=1.0,
        contours=6,
        background_color=None,
        scalarbar=False,
        auto_show=True,
        figure=None,
    ):

        if background_color is None:
            _background_color = self.visualization_data.background_colors["White"]
        else:
            _background_color = background_color

        if figure is None:
            self._backend.create_figure("MO_Square", _background_color, (600, 600))

        if plot_atoms:
            self._plot_atoms(atom_scaling)

        if atom_names:
            self._atom_names(atom_names_scaling)

        if plot_bonds:
            self._plot_bonds(plot_bonds, bond_scaling)

        for number in orbital_numbers:
            X, Y, Z = self.orbital_generator.grid.return_grid_arrays()
            mo_squared = self.molecular_system.MOs[number] ** 2
            contour_values = self.contour_process(contours, mo_squared)
            self._backend.add_contour3d(X, Y, Z, mo_squared, contours=contour_values, opacity=0.5)

        if auto_show:
            self._backend.show()

    def get_orbitals(self, R_max_multip=3.0, x_n=50, y_n=50, z_n=50, gpu=False):

        self.get_orbital_data()
        self.orbital_generator.grid.R_max_multip = R_max_multip
        self.orbital_generator.grid.x_n = x_n
        self.orbital_generator.grid.y_n = y_n
        self.orbital_generator.grid.z_n = z_n
        self.orbital_generator.init_grid()

        self.orbital_generator.init_aos()

        if gpu:
            self.generate_ao_orbitals_gpu()
            self.generate_mo_orbitals_gpu()

        else:
            self.generate_ao_orbitals()
            self.generate_mo_orbitals()

    def get_geminals(self, R_max_multip=3.0, x_n=50, y_n=50, z_n=50, gpu=False):

        self.get_orbitals(R_max_multip=R_max_multip, x_n=x_n, y_n=y_n, z_n=z_n, gpu=gpu)

        self.get_geminal_data()
        self.get_generate_geminals()

    def plot_any_data(
        self,
        any_data_to_plot=None,
        plot_atoms=True,
        atom_scaling=1.0,
        atom_names=True,
        atom_names_scaling=1.0,
        plot_bonds=True,
        bond_scaling=1.0,
        contours=6,
        background_color=None,
        scalarbar=False,
        auto_show=True,
        opacity=0.5,
        figure=None,
    ):

        if background_color is None:
            _background_color = self.visualization_data.background_colors["White"]
        else:
            _background_color = background_color

        if figure is None:
            self._backend.create_figure("Data", _background_color, (600, 600))

        if plot_atoms:
            self._plot_atoms(atom_scaling)

        if atom_names:
            self._atom_names(atom_names_scaling)

        if plot_bonds:
            self._plot_bonds(plot_bonds, bond_scaling)

        X, Y, Z = self.orbital_generator.grid.return_grid_arrays()
        contour_values = self.contour_process(contours, any_data_to_plot)
        self._backend.add_contour3d(
            X, Y, Z, any_data_to_plot, contours=contour_values, opacity=opacity
        )

        if auto_show:
            self._backend.show()

    def dummy_plot(self):
        return 0

    def plot_grid(
        self,
        plot_boundaries=True,
        plot_points=True,
        point_scaling=1.0,
        points_color=(1.0, 0, 0),
        plot_atoms=True,
        atom_scaling=1.0,
        atom_names=True,
        atom_names_scaling=1.0,
        plot_bonds=True,
        bond_scaling=1.0,
        background_color=None,
        auto_show=True,
        opacity=0.5,
        figure=None,
    ):

        if background_color is None:
            _background_color = self.visualization_data.background_colors["White"]
        else:
            _background_color = background_color

        if figure is None:
            self._backend.create_figure("Grid", _background_color, (600, 600))

        if plot_atoms:
            self._plot_atoms(atom_scaling)

        if atom_names:
            self._atom_names(atom_names_scaling)

        if plot_bonds:
            self._plot_bonds(plot_bonds, bond_scaling)

        X, Y, Z = self.orbital_generator.grid.return_grid_arrays()

        point_size = 0.05 * np.sqrt(
            (X[0, 0, 0] - X[1, 0, 0]) ** 2
            + (Y[0, 0, 0] - Y[0, 1, 0]) ** 2
            + (Z[0, 0, 0] - Z[0, 0, 1]) ** 2
        )

        if plot_points:
            # Plot grid points as small spheres
            # Note: This iterates through all grid points which can be slow for large grids
            radius = 0.5 * point_scaling * point_size
            for i in range(X.shape[0]):
                for j in range(X.shape[1]):
                    for k in range(X.shape[2]):
                        self._backend.add_sphere(
                            X[i, j, k],
                            Y[i, j, k],
                            Z[i, j, k],
                            radius=radius,
                            color=points_color,
                            resolution=10,
                        )

        if plot_boundaries:
            self._backend.add_outline()

        if auto_show:
            self._backend.show()
