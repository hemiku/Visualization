import tarfile
from abc import ABC, abstractmethod
from typing import List, Optional, Tuple

import numpy as np
from numpy import float32

from visualization.molecular_system import MolecularSystem
from visualization.basis_normalization import (
    norm_s,
    norm_p,
    norm_d,
    norm_f,
    norm_g,
    norm_h,
    calc_norm_from_basis as _calc_norm,
)


class Input(ABC):
    """Abstract base class for quantum chemistry file parsers.

    This class defines the interface that all input parsers must implement.
    Subclasses handle specific file formats (Dalton, Molpro, Molden, etc.)
    but all provide the same core data: atoms, basis sets, and coefficients.
    """

    input_type: Optional[str] = None
    input_name: Optional[str] = None
    file_string: Optional[str] = None
    BAS_filename: Optional[str] = None
    data_source: Optional[str] = None

    output: Optional[str] = None

    nb: Optional[int] = None
    spherical: bool = False
    nAtoms: Optional[int] = None
    inactive: int = 0
    electrons: Optional[int] = None
    Occ: Optional[np.ndarray] = None

    basis: Optional[List] = None
    basis_norm: Optional[List] = None
    basis_norm2: Optional[List] = None

    Coeff: Optional[np.ndarray] = None

    Atoms_R: Optional[np.ndarray] = None
    Atoms_Charge: Optional[np.ndarray] = None
    Atoms_Name: Optional[List[str]] = None

    Bonds: Optional[List] = None

    def __init__(
        self,
        input_type=None,
        input_sub_type=None,
        input_name=None,
        file_string=None,
        BAS_filename=None,
        data_source=None,
        **kwargs,
    ):
        if input_type is not None:
            self.input_type = input_type
        if input_sub_type is not None:
            self.input_sub_type = input_sub_type
        if input_name is not None:
            self.input_name = input_name
        if file_string is not None:
            self.file_string = file_string
        if BAS_filename is not None:
            self.BAS_filename = BAS_filename
        if data_source is not None:
            self.data_source = data_source

    # Abstract methods - must be implemented by subclasses
    @abstractmethod
    def get_nb(self) -> int:
        """Return number of basis functions."""
        pass

    @abstractmethod
    def get_nAtoms(self) -> int:
        """Return number of atoms."""
        pass

    @abstractmethod
    def get_atoms(self) -> Tuple[np.ndarray, np.ndarray, List[str]]:
        """Return atomic data (positions, charges, names)."""
        pass

    @abstractmethod
    def get_basis(self) -> Tuple[List, List, List]:
        """Return basis set data (basis, basis_norm, basis_norm2)."""
        pass

    @abstractmethod
    def get_coeff(self) -> np.ndarray:
        """Return MO coefficients."""
        pass

    # Optional methods - subclasses may override
    def get_spherical(self) -> bool:
        """Return True if using spherical harmonics."""
        return self.spherical

    def get_inactive(self) -> int:
        """Return number of inactive orbitals."""
        return self.inactive

    def get_electrons(self) -> Optional[int]:
        """Return number of electrons."""
        return self.electrons

    def get_occ(self) -> Optional[np.ndarray]:
        """Return orbital occupation numbers."""
        return self.Occ

    def get_bonds(self) -> Optional[List]:
        """Return bond information."""
        return self.Bonds

    def get_output(self) -> str:
        """Read and return output file contents."""
        if self.output is not None:
            return self.output
        with open(self.input_name + ".out", "r", encoding="utf-8") as f:
            self.output = f.read()
        return self.output

    def calc_norm_from_basis(self, basis: List) -> List:
        """Calculate normalization factors for basis set."""
        return _calc_norm(basis)


class DaltonInput(Input):

    input_name = "Dalton"

    dalton_output = None
    F_BAS = None

    def set_source(self, source):

        self.input_type = source

    def get_data_source(self, data_source):

        self.data_source = data_source

    def get_dalton_output(self):

        if self.dalton_output is not None:

            return self.dalton_output

        else:
            with open(self.input_name + ".out", "r", encoding="utf-8") as f:
                self.dalton_output = f.read()

            return self.dalton_output

    def get_spherical(self):

        Out = self.get_dalton_output()

        if Out.find("Spherical harmonic basis used.") > 0:
            self.spherical = 1
        else:
            self.spherical = 0

        return self.spherical

    def get_nb(self):

        if self.nb is not None:
            return self.nb

        Out = self.get_dalton_output()

        self.nb = int(
            Out[
                Out.find("Number of basis functions") : Out.find("Number of basis functions") + 38
            ].split()[-1]
        )

        return self.nb

    def get_nAtoms(self):

        if self.nAtoms is not None:
            return self.nAtoms

        Out = self.get_dalton_output()

        self.nAtoms = int(
            Out[
                Out.find("Total number of atoms:") : Out.find("Total number of atoms:") + 30
            ].split()[-1]
        )

        return self.nAtoms

    def get_inactive(self):

        if self.inactive is not None:
            return self.inactive

        Out = self.get_dalton_output()

        self.inactive = int(
            Out[
                Out.find("@    Inactive orbitals") : Out.find("@    Inactive orbitals") + 100
            ].split()[3]
        )

        return self.inactive

    def get_electrons(self):

        if self.electrons is not None:
            return self.electrons

        Out = self.get_dalton_output()

        self.electrons = int(
            Out[
                Out.find("@    Number of electrons in active shells") : Out.find(
                    "@    Total charge of the molecule"
                )
                + 100
            ].split()[7]
        )

        return self.electrons

    def get_occ(self):

        if self.Occ is not None:
            return self.Occ

        self.Occ = np.zeros([self.nb], dtype=np.float64)

        Out = self.get_dalton_output()

        Geminal_buf = (
            Out[
                Out.find("APSG geminal coefficients and natural orbital occupations:") : Out.rfind(
                    'NEWORB " orbitals punched.'
                )
            ]
            .split("====================================================================")[1]
            .split("\n")[1:-1]
        )

        for i in range(self.electrons):
            self.Occ[self.inactive + i] = float(Geminal_buf[i].split()[6])

        self.Occ[: self.inactive] = 1

        return self.Occ

    def get_coeff(self):

        F_MOPUN: str

        if self.Coeff is not None:
            return self.Coeff

        if self.input_type == "MOPUN":
            with open("DALTON.MOPUN", "r") as f:
                F_MOPUN = f.read()

        if self.input_type == "tar":
            tar = tarfile.open(self.input_name + ".tar.gz")
            f: str = tar.extractfile(tar.getmember("DALTON.MOPUN"))
            F_MOPUN = f.read().decode(encoding="utf-8")
            tar.close()

        F_MOPUN = F_MOPUN.replace("-", " -")

        a = " ".join(F_MOPUN[F_MOPUN.find("\n") :].split()).split("**NATOCC")[0]
        b = np.fromstring(a, dtype=np.float64, sep=" ")

        self.Coeff = np.reshape(
            np.fromstring(
                " ".join(F_MOPUN[F_MOPUN.find("\n") :].split()).split("**NATOCC")[0],
                dtype=np.float64,
                sep=" ",
            ),
            [self.nb, self.nb],
        )

        return self.Coeff

    def get_basis_file(self) -> str:

        if self.F_BAS is not None:
            return self.F_BAS

        if self.input_type == "MOPUN":
            with open("DALTON.BAS", "r") as f:
                self.F_BAS = f.read()

            return self.F_BAS

        if self.input_type == "tar":
            tar = tarfile.open(self.input_name + ".tar.gz")
            f = tar.extractfile(tar.getmember("DALTON.BAS"))
            self.F_BAS = f.read().decode(encoding="utf-8")
            tar.close()

            return self.F_BAS

    def find_header_end_line(self, lines, Atom_header_pattern=r" {1,9}\d{1,9}\. "):

        import re

        for i, line in enumerate(lines):
            if re.match(Atom_header_pattern, line):
                return i

    def get_atoms(self):

        import re

        self.Atomsa_R = np.zeros([self.nAtoms, 3], dtype=np.float64)
        self.Atoms_Charge = np.zeros(self.nAtoms, dtype=np.int64)
        self.Atoms_Name = []

        F_BAS = self.get_basis_file()
        F_BAS_split_lines = F_BAS.splitlines()

        Atom_header_pattern = r" {1,9}\d{1,9}\. "
        Atom_Geometries_pattern = r"^\w{1,6} \D \D \D"
        Basis_header_pattern = r"^H {1,3}\d{1,3} {1,4}\d{1,3}$"

        Atoms_Groups = []

        header_end_line = self.find_header_end_line(
            lines=F_BAS_split_lines, Atom_header_pattern=Atom_header_pattern
        )

        for line in F_BAS_split_lines[header_end_line:]:

            if re.match(Atom_header_pattern, line):

                Atoms_Group = {"headerLine": line, "Geometries": [], "Basis": []}
                Atoms_Groups.append(Atoms_Group)

            elif re.match(Atom_Geometries_pattern, line):

                Atoms_Group["Geometries"].append(line)

            elif re.match(Basis_header_pattern, line):

                basis_part = {"header": line, "data": []}
                Atoms_Group["Basis"].append(basis_part)

            else:
                pass
                # basis_part['data'].append(line)

        self.map_atoms_geometry_from_BAS_file(Atoms_Groups)

        return self.Atoms_R, self.Atoms_Charge, self.Atoms_Name

    def get_bonds(self):

        self.Bonds = []

        Out = self.get_dalton_output()

        start_bond_section = "Bond distances "
        start_next_section = "| Starting in Integral Section (HERMIT) |"

        Bonds_str = Out[Out.find(start_bond_section) : Out.find(start_next_section)]
        Bonds_str = Bonds_str[Bonds_str.find("bond") :]
        Bonds_str_lines = Bonds_str.splitlines()

        for bond_str in Bonds_str_lines:
            bond_str_split = bond_str.split()

            try:
                self.Bonds.append(
                    [
                        self.Atoms_Name.index(bond_str_split[2]),
                        self.Atoms_Name.index(bond_str_split[3]),
                        float(bond_str_split[4]),
                    ]
                )
                continue
            except:
                pass

            try:
                self.Bonds.append(
                    [
                        self.Atoms_Name.index(bond_str_split[2]),
                        self.Atoms_Name.index(bond_str_split[4]),
                        float(bond_str_split[6]),
                    ]
                )
                continue
            except:
                pass

            if bond_str == "":
                break

        return self.Bonds

    def get_basis(self):

        import re

        self.basis = []
        self.basis_norm = []
        self.basis_norm2 = []

        F_BAS = self.get_basis_file()
        F_BAS_split_lines = F_BAS.splitlines()

        Atom_header_pattern = r" {1,9}\d{1,9}\. "
        Atom_Geometries_pattern = r"^\w{1,6} \D \D \D"
        Basis_header_pattern = r"^H {1,3}\d{1,3} {1,4}\d{1,3}$"

        Atoms_Groups: list = []
        Atoms_Group: dict

        header_end_line = self.find_header_end_line(
            lines=F_BAS_split_lines, Atom_header_pattern=Atom_header_pattern
        )

        for line in F_BAS_split_lines[header_end_line:]:

            if re.match(Atom_header_pattern, line):

                Atoms_Group: dict = {"headerLine": line, "Geometries": [], "Basis": []}
                Atoms_Groups.append(Atoms_Group)

            elif re.match(Atom_Geometries_pattern, line):

                Atoms_Group["Geometries"].append(line)

            elif re.match(Basis_header_pattern, line):

                basis_part = {"header": line, "data": []}
                Atoms_Group["Basis"].append(basis_part)

            else:
                basis_part["data"].append(line)

        self.map_atoms_from_BAS_file(Atoms_Groups)

        return self.basis, self.basis_norm, self.basis_norm2

    def map_atoms_geometry_from_BAS_file(self, Atoms_Groups):
        """
        docstring
        """

        self.Atoms_R = np.zeros([self.nAtoms, 3], dtype=np.float64)
        self.Atoms_Charge = np.zeros(self.nAtoms, dtype=np.int64)
        self.Atoms_Name = []

        i = 0

        for Atoms_Group in Atoms_Groups:

            Atom_Charge = int(Atoms_Group["headerLine"][: Atoms_Group["headerLine"].find(".")])

            for Atom in Atoms_Group["Geometries"]:
                Atom_split = Atom.split()

                Atom_name = Atom_split[0]
                Atom_R = np.array(
                    [float(Atom_split[1]), float(Atom_split[2]), float(Atom_split[3])]
                )

                self.Atoms_Name.append(Atom_name)
                self.Atoms_Charge[i] = Atom_Charge
                self.Atoms_R[i] = Atom_R

                i += 1

    def map_atoms_from_BAS_file(self, Atoms_Groups):
        """
        docstring
        """

        self.Atoms_R = np.zeros([self.nAtoms, 3], dtype=np.float64)
        self.Atoms_Charge = np.zeros(self.nAtoms, dtype=np.int64)
        self.Atoms_Name = []

        self.basis = []
        self.basis_norm = []
        self.basis_norm2 = []

        i = 0

        for Atoms_Group in Atoms_Groups:

            # TODO basis
            Orbitals = []

            for Orbital in Atoms_Group["Basis"]:
                dim = np.fromstring(Orbital["header"][1:], dtype=np.int64, sep=" ")
                dim[1] += 1

                if len(Orbital["data"]):
                    Orbital = np.reshape(
                        np.fromstring("".join(Orbital["data"]), dtype=np.float64, sep=" "), dim
                    )
                    Orbitals.append(Orbital)

            Atom_Charge = int(Atoms_Group["headerLine"][: Atoms_Group["headerLine"].find(".")])

            for Atom in Atoms_Group["Geometries"]:
                Atom_split = Atom.split()

                Atom_name = Atom_split[0]
                Atom_R = np.array(
                    [float(Atom_split[1]), float(Atom_split[2]), float(Atom_split[3])]
                )

                self.Atoms_Name.append(Atom_name)
                self.Atoms_Charge[i] = Atom_Charge
                self.Atoms_R[i] = Atom_R

                if len(Orbitals):
                    self.basis.append(Orbitals)

                i += 1

        if len(self.basis):
            norm_funcs = [norm_s, norm_p, norm_d, norm_f, norm_g, norm_h]
            for n in range(self.nAtoms):
                self.basis_norm.append([])
                for j in range(len(self.basis[n])):
                    if j < len(norm_funcs):
                        self.basis_norm[n].append(norm_funcs[j](self.basis[n][j]))

    def get_printout_of_final_geminals(self, dalton_output):

        GEMINAL_PART_START = "Printout of final geminals"
        WAWE_FUNCTION_SECTION_END = "| End of Wave Function Section (SIRIUS) |"

        _geminal_part = dalton_output[
            dalton_output.find(GEMINAL_PART_START) : dalton_output.find(WAWE_FUNCTION_SECTION_END)
        ]

        return _geminal_part

    def get_g_coeff(self):

        Coeff_separator = "===================================================================="
        Out = self.get_dalton_output()
        Geminal_buf = self.get_printout_of_final_geminals(dalton_output=Out)
        coeff_buf = Geminal_buf.split(Coeff_separator)[1].split("\n")[1:-1]

        G_coeff = []

        for coeff in coeff_buf:
            G_coeff.append(float(coeff.split()[5]))

        self.G_coeff = np.array(G_coeff, dtype=np.float32)

        return self.G_coeff

    def get_orb2gem(self):

        Coeff_separator = "===================================================================="
        Out = self.get_dalton_output()
        Geminal_buf = self.get_printout_of_final_geminals(dalton_output=Out)
        coeff_buf = Geminal_buf.split(Coeff_separator)[1].split("\n")[1:-1]

        Orb2Gem = []

        for coeff in coeff_buf:
            Orb2Gem.append(int(coeff.split()[3]) - 1)

        self.Orb2Gem = np.array(Orb2Gem, dtype=np.int32)
        self.nGeminal = int(np.max(self.Orb2Gem)) + 1

        # self.Orb2Gem = np.zeros([self.electrons], dtype=np.int64)

        # Out = self.get_dalton_output()
        # Geminal_buf = Out[Out.find("APSG geminal coefficients and natural orbital occupations:"):
        #                   Out.rfind('NEWORB " orbitals punched.')].split(
        #     "====================================================================")[1].split("\n")[1:-1]

        # for i in range(self.electrons):
        #     self.Orb2Gem[i] = int(Geminal_buf[i].split()[3]) - 1

        # self.nGeminal = int(len(self.Orb2Gem) / 2)

        return self.Orb2Gem, self.nGeminal


class MoldenInput(Input):

    input_name = "molden"

    output = None

    def get_output(self, filename):
        file_type = ".molden"
        if filename[-len(file_type) :] == file_type:
            self.file_name = filename[: -len(file_type)]
        else:
            self.file_name = filename

    def get_file(self):

        if self.output is not None:

            return self.output

        else:
            with open(self.input_name, "r", encoding="utf-8") as f:
                self.output = f.read()

            return self.output

    def get_data(self):

        Atoms_begin = "Atoms]"
        GTO_begin = "GTO]"
        MO_begin = "MO]"

        Out = self.get_file()

        if Out.find("[Molden Format]") >= 0:
            Out_sections = Out.split("[")

            sections = {}

            for section in Out_sections:

                section_key = section[: section.find("]")]
                section_data = section[section.find("]") + 1 :]

                sections[section_key] = section_data

            self.get_atoms(sections["Atoms"])
            self.get_GTOs(sections["GTO"])
            self.get_MOs(sections["MO"])

        self.Spherical = True
        self.calc_nb()

    def calc_nb(self):
        self.nb = len(self.MOs)

    def get_mos(self, section):

        MO_sym = "Sym"
        MO_Ene = "Ene"
        MO_Spin = "Spin"
        MO_Occup = "Occup"

        self.MOs_list = []
        section_lines = section.splitlines()
        for j, line in enumerate(section_lines[1:]):
            if line.count(MO_Occup):
                # MO = []
                # self.MOs_list.append(MO)
                MO_dict = {}

            elif (
                line.count(MO_sym) + line.count(MO_Ene) + line.count(MO_Spin) + line.count(MO_Occup)
            ) == 0:
                # MO.append( float( line.split()[1] ) )

                MO_dict[int(line.split()[0])] = float(line.split()[1])

        self.MOs = np.zeros([len(self.MOs_list), len(self.MOs_list)], dtype=np.float64)

        for i, MO_dict in enumerate(self.MOs_list):
            for j in MO_dict:
                self.MOs[i, j] = MO_dict[j]

        # self.MOs = np.asarray( self.MOs_list )

    def get_gtos(self, section):
        self.GTOs = []
        section_lines = section.splitlines()
        for j, line in enumerate(section_lines[1:]):
            line_split = line.split()
            if len(line_split) == 2 and len(line) < 18:
                atom_GTOs = [[], [], [], [], [], []]
                self.GTOs.append(atom_GTOs)
            elif len(line_split) == 3 and line_split[0] == "s":
                GTO = []
                harmonic_GTOs = atom_GTOs[0]
                harmonic_GTOs.append(GTO)
            #            elif line[:2] == ' p':
            elif len(line_split) == 3 and line_split[0] == "p":
                GTO = []
                harmonic_GTOs = atom_GTOs[1]
                harmonic_GTOs.append(GTO)
            # elif line[:2] == ' d':
            elif len(line_split) == 3 and line_split[0] == "d":
                GTO = []
                harmonic_GTOs = atom_GTOs[2]
                harmonic_GTOs.append(GTO)
            # elif line[:2] == ' f':
            elif len(line_split) == 3 and line_split[0] == "f":
                GTO = []
                harmonic_GTOs = atom_GTOs[3]
                harmonic_GTOs.append(GTO)
            #            elif line[:2] == '  ' and len(line) > 2:
            elif len(line_split) == 2 and line[18] == " ":

                GTO.append(
                    [float(line_split[0].replace("D", "E")), float(line_split[1].replace("D", "E"))]
                )

        self.basis = []

        for atom_GTOs in self.GTOs:
            atom_basis = []
            self.basis.append(atom_basis)
            for harmonic_GTOs in atom_GTOs:
                if len(harmonic_GTOs) > 0:
                    atom_basis.append(self.transform_molden_basis(harmonic_GTOs))

    def transform_molden_basis(self, harmonic_GTOs):

        Sigmas_list = []
        number_of_orbital = len(harmonic_GTOs)

        for GTO in harmonic_GTOs:
            for GTO_line in GTO:
                Sigmas_list.append(GTO_line[0])

        Sigmas = np.array(sorted(set(Sigmas_list), reverse=True))

        basis = np.zeros([Sigmas.shape[0], number_of_orbital + 1], dtype=np.float64)
        basis[:, 0] = Sigmas
        for i, GTO in enumerate(harmonic_GTOs):
            for GTO_line in GTO:
                mask = Sigmas == GTO_line[0]
                basis[mask, i + 1] = GTO_line[1]

        return basis

    def get_atoms(self, section):
        self.Atoms = []
        section_lines = section.splitlines()
        for j, line in enumerate(section_lines[1:]):
            self.Atoms.append(line)

        self.nAtoms = len(self.Atoms)

        self.Atoms_R = np.zeros([self.nAtoms, 3], dtype=np.float64)
        self.Atoms_Charge = np.zeros(self.nAtoms, dtype=np.int64)
        self.Atoms_Name = []

        for i, atom in enumerate(self.Atoms):

            atom_split = atom.split()
            self.Atoms_Name.append(atom_split[0])
            self.Atoms_Charge[i] = int(atom_split[2])
            self.Atoms_R[i] = np.array(
                [float(atom_split[-3]), float(atom_split[-2]), float(atom_split[-1])]
            )

    def get_nb(self):

        if self.nb is None:
            self.get_data()
        return self.nb

    def get_ao_number(self):
        return self.MOs.shape[1]

    def get_atoms_r(self):
        return self.Atoms_R

    def get_atoms_name(self):
        return self.Atoms_Name

    def get_atoms_charge(self):
        return self.Atoms_Charge

    def get_basis(self):
        # return self.Basis
        return self.basis, self.basis_norm, self.basis_norm2

    def get_spherical(self):
        return self.spherical

    def get_nAtoms(self):
        return self.nAtoms

    def get_coeff(self):
        return self.Coeff

    def get_atoms(self):

        self.get_data()
        return self.Atoms_R, self.Atoms_Charge, self.Atoms_Name


MolproInput = None
INPUT_TYPES = {"Dalton": DaltonInput, "molden": MoldenInput, "Molpro": MolproInput}


def get_input(
    input_type: str,
    input_sub_type: str,
    input_name=None,
    file_string=None,
    BAS_filename=None,
    data_source=None,
):

    try:
        data_input = INPUT_TYPES[input_type](
            input_type=input_sub_type,
            input_name=input_name,
            file_string=file_string,
            BAS_filename=BAS_filename,
            data_source=data_source,
        )
        return data_input

    except KeyError:
        raise Exception(f"Input_type: {input_type} not found")
