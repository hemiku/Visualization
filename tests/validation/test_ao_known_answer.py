"""Known-answer validation tests for atomic orbital calculations.

These tests verify that atomic orbitals computed by OrbitalsGenerator
match analytical reference values for Gaussian-type orbitals.

Analytical Reference Values (alpha=1.0, coefficient=1.0):
=========================================================

S-orbital (l=0):
- Normalization: N_s = (2*alpha/pi)^(3/4)
- At origin: phi_s(0,0,0) = (2/pi)^(3/4) = 0.7127054703549901
- At (1,0,0): phi_s = N_s * exp(-1) = 0.2621589005446064

P-orbital (l=1):
- Normalization: N_p = 2 * (2*alpha/pi)^(3/4) * alpha^(1/2) = 2 * (2/pi)^(3/4)
- N_p(alpha=1) = 1.4254109407099802
- At origin: phi_p = 0 (node)
- Px at (1,0,0): N_p * x * exp(-alpha*r^2) = N_p * exp(-1) = 0.5243178010892128

D-orbital (l=2, Cartesian):
- For d_xy: phi = N_d * x * y * exp(-alpha*r^2)
- At origin: phi_d = 0 (node)
- Symmetry: d_xy(x,y,z) = d_xy(-x,-y,z) (even in xy)

F-orbital (l=3):
- At origin: phi_f = 0 (node)
- Various angular symmetries
"""

import numpy as np
import pytest
from scipy import integrate

from visualization.grid import Grid
from visualization.orbitals import OrbitalsGenerator


class TestSOrbitalAtOrigin:
    """Test S-orbital values at the origin against analytical reference."""

    @pytest.fixture
    def simple_s_orbital_system(self):
        """Create minimal system with one S-orbital at origin.

        Single atom at origin with one S-type primitive:
        - alpha (exponent) = 1.0
        - coefficient = 1.0
        """
        # Single atom at origin
        nAtoms = 1
        atoms_R = np.array([[0.0, 0.0, 0.0]])
        spherical = 1
        nb = 1  # One basis function

        # Basis: [[s_data], [p_data], [d_data], ...]
        # s_data has shape (n_primitives, 1 + n_contractions)
        # Format: [exponent, coeff1, coeff2, ...]
        s_data = np.array([[1.0, 1.0]])  # alpha=1.0, coeff=1.0
        basis = [[s_data]]  # One atom with only S orbital

        # Normalization factors (not used for S orbitals in calc_aos)
        basis_norm = [[[1.0]]]

        # Create grid centered at origin with origin as a grid point
        # Using odd number ensures origin is exactly on a grid point
        grid = Grid(x_n=21, y_n=21, z_n=21)
        grid.x_min, grid.x_max = -5.0, 5.0
        grid.y_min, grid.y_max = -5.0, 5.0
        grid.z_min, grid.z_max = -5.0, 5.0

        return {
            "nAtoms": nAtoms,
            "atoms_R": atoms_R,
            "spherical": spherical,
            "nb": nb,
            "basis": basis,
            "basis_norm": basis_norm,
            "grid": grid,
        }

    def test_s_orbital_value_at_origin(self, simple_s_orbital_system):
        """Test S-orbital (alpha=1.0) gives (2/pi)^(3/4) at origin."""
        config = simple_s_orbital_system

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        # Find origin indices (center of grid)
        center_idx = config["grid"].x_n // 2

        # Expected value: N = (2*alpha/pi)^(3/4) for alpha=1.0
        expected = (2.0 / np.pi) ** 0.75  # = 0.7127054703549901
        actual = gen.AOs[0, center_idx, center_idx, center_idx]

        assert np.isclose(
            actual, expected, rtol=1e-10
        ), f"S-orbital at origin: expected {expected}, got {actual}"

    def test_s_orbital_value_at_origin_alpha_2(self):
        """Test S-orbital (alpha=2.0) gives (4/pi)^(3/4) at origin."""
        alpha = 2.0
        nAtoms = 1
        atoms_R = np.array([[0.0, 0.0, 0.0]])

        s_data = np.array([[alpha, 1.0]])
        basis = [[s_data]]
        basis_norm = [[[1.0]]]

        grid = Grid(x_n=21, y_n=21, z_n=21)
        grid.x_min, grid.x_max = -5.0, 5.0
        grid.y_min, grid.y_max = -5.0, 5.0
        grid.z_min, grid.z_max = -5.0, 5.0

        gen = OrbitalsGenerator(
            nAtoms=nAtoms,
            atoms_R=atoms_R,
            spherical=1,
            nb=1,
            basis=basis,
            basis_norm=basis_norm,
            grid=grid,
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        center_idx = grid.x_n // 2
        expected = (2.0 * alpha / np.pi) ** 0.75
        actual = gen.AOs[0, center_idx, center_idx, center_idx]

        assert np.isclose(actual, expected, rtol=1e-10)

    def test_s_orbital_decay_with_distance(self, simple_s_orbital_system):
        """Test S-orbital decays as exp(-alpha*r^2) from origin."""
        config = simple_s_orbital_system

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        # Get values along x-axis (y=0, z=0)
        center_idx = config["grid"].x_n // 2

        # Value at origin
        val_origin = gen.AOs[0, center_idx, center_idx, center_idx]

        # Value at x=1 (one grid step from center depends on grid spacing)
        x_vals = np.linspace(config["grid"].x_min, config["grid"].x_max, config["grid"].x_n)
        dx = x_vals[1] - x_vals[0]

        # Value at some distance from origin
        offset = 2  # 2 grid points from center
        val_at_offset = gen.AOs[0, center_idx + offset, center_idx, center_idx]

        x_at_offset = offset * dx
        expected_ratio = np.exp(-1.0 * x_at_offset**2)  # alpha=1.0
        actual_ratio = val_at_offset / val_origin

        assert np.isclose(actual_ratio, expected_ratio, rtol=1e-8)


class TestSOrbitalNormalization:
    """Test numerical normalization of S-orbitals."""

    def test_s_orbital_integral_approximately_one(self):
        """Test that integral of |phi_s|^2 dV ≈ 1 for normalized S-orbital."""
        nAtoms = 1
        atoms_R = np.array([[0.0, 0.0, 0.0]])

        # Use alpha=1.0 for simplicity
        s_data = np.array([[1.0, 1.0]])
        basis = [[s_data]]
        basis_norm = [[[1.0]]]

        # Use a fine grid for accurate integration
        # Grid must extend far enough to capture most of the density
        grid = Grid(x_n=51, y_n=51, z_n=51)
        grid.x_min, grid.x_max = -6.0, 6.0
        grid.y_min, grid.y_max = -6.0, 6.0
        grid.z_min, grid.z_max = -6.0, 6.0

        gen = OrbitalsGenerator(
            nAtoms=nAtoms,
            atoms_R=atoms_R,
            spherical=1,
            nb=1,
            basis=basis,
            basis_norm=basis_norm,
            grid=grid,
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        # Compute numerical integral using trapezoidal rule
        dx = (grid.x_max - grid.x_min) / (grid.x_n - 1)
        dy = (grid.y_max - grid.y_min) / (grid.y_n - 1)
        dz = (grid.z_max - grid.z_min) / (grid.z_n - 1)
        dv = dx * dy * dz

        integral = np.sum(gen.AOs[0] ** 2) * dv

        # For a normalized Gaussian, integral should be 1
        # Allow some tolerance due to grid truncation
        assert np.isclose(
            integral, 1.0, rtol=0.05
        ), f"S-orbital normalization integral: expected ~1.0, got {integral}"


class TestPOrbitalAnalyticalValues:
    """Test P-orbital values at specific points against analytical reference.

    P-orbital formula: phi_p = N_p * x * exp(-alpha * r^2)
    where N_p = 2 * (2*alpha/pi)^(3/4) * alpha^(1/2)

    For alpha=1.0: N_p = 2 * (2/pi)^(3/4) = 1.4254109407099802
    """

    @pytest.fixture
    def p_orbital_system(self):
        """Create system with S and P orbitals."""
        nAtoms = 1
        atoms_R = np.array([[0.0, 0.0, 0.0]])
        spherical = 1
        nb = 4  # 1 S + 3 P

        s_data = np.array([[1.0, 1.0]])
        p_data = np.array([[1.0, 1.0]])
        basis = [[s_data, p_data]]
        basis_norm = [[[1.0], [1.0]]]

        # Fine grid with x=1 as a grid point
        grid = Grid(x_n=41, y_n=41, z_n=41)
        grid.x_min, grid.x_max = -5.0, 5.0
        grid.y_min, grid.y_max = -5.0, 5.0
        grid.z_min, grid.z_max = -5.0, 5.0

        return {
            "nAtoms": nAtoms,
            "atoms_R": atoms_R,
            "spherical": spherical,
            "nb": nb,
            "basis": basis,
            "basis_norm": basis_norm,
            "grid": grid,
        }

    def test_p_orbital_normalization_constant(self):
        """Verify P-orbital normalization constant: N_p = 2*(2/pi)^(3/4) for alpha=1."""
        alpha = 1.0
        N_p_expected = 2.0 * (2.0 * alpha / np.pi) ** 0.75 * (alpha**0.5)
        assert np.isclose(N_p_expected, 1.4254109407099802, rtol=1e-10)

    def test_px_orbital_at_x_equals_1(self, p_orbital_system):
        """Test Px at (1,0,0): N_p * 1 * exp(-1) = 0.5243178010892128."""
        config = p_orbital_system

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        # Find grid point closest to x=1
        x_vals = np.linspace(config["grid"].x_min, config["grid"].x_max, config["grid"].x_n)
        x_idx = np.argmin(np.abs(x_vals - 1.0))
        center_y = config["grid"].y_n // 2
        center_z = config["grid"].z_n // 2

        # Analytical value: N_p * x * exp(-alpha * x^2) for x=1, alpha=1
        N_p = 2.0 * (2.0 / np.pi) ** 0.75
        x_actual = x_vals[x_idx]
        expected = N_p * x_actual * np.exp(-1.0 * x_actual**2)

        # Px is orbital index 1 (after S)
        actual = gen.AOs[1, x_idx, center_y, center_z]

        assert np.isclose(
            actual, expected, rtol=1e-6
        ), f"Px at ({x_actual},0,0): expected {expected}, got {actual}"

    def test_py_orbital_at_y_equals_1(self, p_orbital_system):
        """Test Py at (0,1,0): N_p * 1 * exp(-1)."""
        config = p_orbital_system

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        y_vals = np.linspace(config["grid"].y_min, config["grid"].y_max, config["grid"].y_n)
        y_idx = np.argmin(np.abs(y_vals - 1.0))
        center_x = config["grid"].x_n // 2
        center_z = config["grid"].z_n // 2

        N_p = 2.0 * (2.0 / np.pi) ** 0.75
        y_actual = y_vals[y_idx]
        expected = N_p * y_actual * np.exp(-1.0 * y_actual**2)

        # Py is orbital index 2
        actual = gen.AOs[2, center_x, y_idx, center_z]

        assert np.isclose(actual, expected, rtol=1e-6)

    def test_pz_orbital_at_z_equals_1(self, p_orbital_system):
        """Test Pz at (0,0,1): N_p * 1 * exp(-1)."""
        config = p_orbital_system

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        z_vals = np.linspace(config["grid"].z_min, config["grid"].z_max, config["grid"].z_n)
        z_idx = np.argmin(np.abs(z_vals - 1.0))
        center_x = config["grid"].x_n // 2
        center_y = config["grid"].y_n // 2

        N_p = 2.0 * (2.0 / np.pi) ** 0.75
        z_actual = z_vals[z_idx]
        expected = N_p * z_actual * np.exp(-1.0 * z_actual**2)

        # Pz is orbital index 3
        actual = gen.AOs[3, center_x, center_y, z_idx]

        assert np.isclose(actual, expected, rtol=1e-6)

    def test_p_orbital_maximum_location(self, p_orbital_system):
        """Test Px maximum is at x = 1/sqrt(2*alpha) for alpha=1."""
        config = p_orbital_system

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        # For Px = N_p * x * exp(-alpha*x^2), max is at x = 1/sqrt(2*alpha)
        # For alpha=1: x_max = 1/sqrt(2) = 0.7071...
        alpha = 1.0
        x_max_expected = 1.0 / np.sqrt(2.0 * alpha)

        x_vals = np.linspace(config["grid"].x_min, config["grid"].x_max, config["grid"].x_n)
        center_y = config["grid"].y_n // 2
        center_z = config["grid"].z_n // 2

        # Find max along positive x-axis
        px_along_x = gen.AOs[1, :, center_y, center_z]
        x_max_idx = np.argmax(px_along_x)
        x_max_actual = x_vals[x_max_idx]

        assert np.isclose(
            x_max_actual, x_max_expected, atol=0.3
        ), f"Px maximum at x={x_max_actual}, expected {x_max_expected}"


class TestDOrbitalAnalyticalValues:
    """Test D-orbital values at specific points (Cartesian basis).

    For Cartesian D-orbitals (spherical=0):
    - d_xx: N_d * x^2 * exp(-alpha*r^2)
    - d_xy: N_d * sqrt(3) * x * y * exp(-alpha*r^2)
    etc.
    """

    @pytest.fixture
    def d_orbital_system_cartesian(self):
        """Create system with S, P, D orbitals in Cartesian basis."""
        nAtoms = 1
        atoms_R = np.array([[0.0, 0.0, 0.0]])
        spherical = 0  # Cartesian basis
        nb = 10  # 1 S + 3 P + 6 D (Cartesian)

        s_data = np.array([[1.0, 1.0]])
        p_data = np.array([[1.0, 1.0]])
        d_data = np.array([[1.0, 1.0]])
        basis = [[s_data, p_data, d_data]]
        basis_norm = [[[1.0], [1.0], [1.0]]]

        grid = Grid(x_n=41, y_n=41, z_n=41)
        grid.x_min, grid.x_max = -4.0, 4.0
        grid.y_min, grid.y_max = -4.0, 4.0
        grid.z_min, grid.z_max = -4.0, 4.0

        return {
            "nAtoms": nAtoms,
            "atoms_R": atoms_R,
            "spherical": spherical,
            "nb": nb,
            "basis": basis,
            "basis_norm": basis_norm,
            "grid": grid,
        }

    def test_d_orbital_zero_at_origin(self, d_orbital_system_cartesian):
        """All D-orbitals should be zero at origin."""
        config = d_orbital_system_cartesian

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        center = config["grid"].x_n // 2

        # D orbitals are indices 4-9 (after 1 S + 3 P)
        for d_idx in range(4, 10):
            value = gen.AOs[d_idx, center, center, center]
            assert np.isclose(
                value, 0.0, atol=1e-12
            ), f"D-orbital {d_idx} at origin: expected 0, got {value}"

    def test_d_xy_symmetry(self, d_orbital_system_cartesian):
        """Test d_xy symmetry: d_xy(x,y) = d_xy(-x,-y) and d_xy(x,-y) = -d_xy(x,y)."""
        config = d_orbital_system_cartesian

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        center = config["grid"].x_n // 2
        offset = 3

        # d_xy is orbital index 7 (after xx=4, yy=5, zz=6, xy=7)
        # Actually the order is: xx, yy, zz, xy, xz, yz
        d_xy_idx = 7

        # d_xy(x,y,0) should equal d_xy(-x,-y,0)
        val_pp = gen.AOs[d_xy_idx, center + offset, center + offset, center]
        val_mm = gen.AOs[d_xy_idx, center - offset, center - offset, center]
        assert np.isclose(
            val_pp, val_mm, rtol=1e-10
        ), f"d_xy symmetry: d_xy(+,+)={val_pp}, d_xy(-,-)={val_mm}"

        # d_xy(x,-y,0) should equal -d_xy(x,y,0)
        val_pm = gen.AOs[d_xy_idx, center + offset, center - offset, center]
        assert np.isclose(
            val_pm, -val_pp, rtol=1e-10
        ), f"d_xy antisymmetry: d_xy(+,-)={val_pm}, -d_xy(+,+)={-val_pp}"


class TestDOrbitalSpherical:
    """Test D-orbital values in spherical harmonic basis."""

    @pytest.fixture
    def d_orbital_system_spherical(self):
        """Create system with D orbitals in spherical basis."""
        nAtoms = 1
        atoms_R = np.array([[0.0, 0.0, 0.0]])
        spherical = 1  # Spherical harmonic basis
        nb = 9  # 1 S + 3 P + 5 D (spherical)

        s_data = np.array([[1.0, 1.0]])
        p_data = np.array([[1.0, 1.0]])
        d_data = np.array([[1.0, 1.0]])
        basis = [[s_data, p_data, d_data]]
        basis_norm = [[[1.0], [1.0], [1.0]]]

        grid = Grid(x_n=41, y_n=41, z_n=41)
        grid.x_min, grid.x_max = -4.0, 4.0
        grid.y_min, grid.y_max = -4.0, 4.0
        grid.z_min, grid.z_max = -4.0, 4.0

        return {
            "nAtoms": nAtoms,
            "atoms_R": atoms_R,
            "spherical": spherical,
            "nb": nb,
            "basis": basis,
            "basis_norm": basis_norm,
            "grid": grid,
        }

    def test_d_spherical_zero_at_origin(self, d_orbital_system_spherical):
        """All spherical D-orbitals should be zero at origin."""
        config = d_orbital_system_spherical

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        center = config["grid"].x_n // 2

        # D orbitals are indices 4-8 (after 1 S + 3 P, 5 spherical D)
        for d_idx in range(4, 9):
            value = gen.AOs[d_idx, center, center, center]
            assert np.isclose(
                value, 0.0, atol=1e-12
            ), f"Spherical D-orbital {d_idx} at origin: expected 0, got {value}"

    def test_d_z2_axial_symmetry(self, d_orbital_system_spherical):
        """Test d_z2 (m=0) has axial symmetry: same value at (x,0,z) and (0,y,z)."""
        config = d_orbital_system_spherical

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        center = config["grid"].x_n // 2
        offset = 3

        # d_z2 (m=0) is orbital index 6 (after S, Px, Py, Pz, d_{-2}, d_{-1})
        # Order: m = -2, -1, 0, 1, 2 -> indices 4, 5, 6, 7, 8
        d_z2_idx = 6

        # Same distance from z-axis should give same value
        # (x,0,z) and (0,x,z) at same r should be equal
        val_xz = gen.AOs[d_z2_idx, center + offset, center, center + offset]
        val_yz = gen.AOs[d_z2_idx, center, center + offset, center + offset]

        assert np.isclose(
            val_xz, val_yz, rtol=1e-8
        ), f"d_z2 axial symmetry: val(x,0,z)={val_xz}, val(0,y,z)={val_yz}"


class TestPOrbitalSymmetry:
    """Test P-orbital symmetry properties."""

    @pytest.fixture
    def simple_p_orbital_system(self):
        """Create minimal system with P-orbitals at origin."""
        nAtoms = 1
        atoms_R = np.array([[0.0, 0.0, 0.0]])
        spherical = 1
        nb = 4  # 1 S + 3 P

        # S and P orbitals
        s_data = np.array([[1.0, 1.0]])
        p_data = np.array([[1.0, 1.0]])
        basis = [[s_data, p_data]]
        basis_norm = [[[1.0], [1.0]]]

        # Grid centered at origin
        grid = Grid(x_n=21, y_n=21, z_n=21)
        grid.x_min, grid.x_max = -5.0, 5.0
        grid.y_min, grid.y_max = -5.0, 5.0
        grid.z_min, grid.z_max = -5.0, 5.0

        return {
            "nAtoms": nAtoms,
            "atoms_R": atoms_R,
            "spherical": spherical,
            "nb": nb,
            "basis": basis,
            "basis_norm": basis_norm,
            "grid": grid,
        }

    def test_px_orbital_antisymmetry(self, simple_p_orbital_system):
        """Test Px orbital: phi(-x,y,z) = -phi(x,y,z)."""
        config = simple_p_orbital_system

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        # Px is orbital index 1 (after S)
        px = gen.AOs[1]

        # Check antisymmetry: px(-x) = -px(x)
        center = config["grid"].x_n // 2
        offset = 3

        # px at (x, 0, 0) and (-x, 0, 0)
        px_positive = px[center + offset, center, center]
        px_negative = px[center - offset, center, center]

        assert np.isclose(
            px_positive, -px_negative, rtol=1e-10
        ), f"Px antisymmetry failed: px(+x)={px_positive}, px(-x)={px_negative}"

    def test_py_orbital_antisymmetry(self, simple_p_orbital_system):
        """Test Py orbital: phi(x,-y,z) = -phi(x,y,z)."""
        config = simple_p_orbital_system

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        # Py is orbital index 2
        py = gen.AOs[2]

        center = config["grid"].y_n // 2
        offset = 3

        py_positive = py[center, center + offset, center]
        py_negative = py[center, center - offset, center]

        assert np.isclose(py_positive, -py_negative, rtol=1e-10)

    def test_pz_orbital_antisymmetry(self, simple_p_orbital_system):
        """Test Pz orbital: phi(x,y,-z) = -phi(x,y,z)."""
        config = simple_p_orbital_system

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        # Pz is orbital index 3
        pz = gen.AOs[3]

        center = config["grid"].z_n // 2
        offset = 3

        pz_positive = pz[center, center, center + offset]
        pz_negative = pz[center, center, center - offset]

        assert np.isclose(pz_positive, -pz_negative, rtol=1e-10)

    def test_p_orbital_zero_at_origin(self, simple_p_orbital_system):
        """Test all P-orbitals are zero at the origin."""
        config = simple_p_orbital_system

        gen = OrbitalsGenerator(
            nAtoms=config["nAtoms"],
            atoms_R=config["atoms_R"],
            spherical=config["spherical"],
            nb=config["nb"],
            basis=config["basis"],
            basis_norm=config["basis_norm"],
            grid=config["grid"],
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        center = config["grid"].x_n // 2

        # All P orbitals (indices 1, 2, 3) should be zero at origin
        for p_idx in [1, 2, 3]:
            value_at_origin = gen.AOs[p_idx, center, center, center]
            assert np.isclose(
                value_at_origin, 0.0, atol=1e-14
            ), f"P-orbital {p_idx} at origin: expected 0, got {value_at_origin}"


class TestGridConvergence:
    """Test that orbital values converge with increasing grid resolution."""

    def test_s_orbital_origin_converges(self):
        """Test S-orbital at origin converges to analytical value."""
        expected = (2.0 / np.pi) ** 0.75

        nAtoms = 1
        atoms_R = np.array([[0.0, 0.0, 0.0]])
        s_data = np.array([[1.0, 1.0]])
        basis = [[s_data]]
        basis_norm = [[[1.0]]]

        # Test with increasing grid resolution
        resolutions = [11, 21, 41]
        values = []

        for res in resolutions:
            grid = Grid(x_n=res, y_n=res, z_n=res)
            grid.x_min, grid.x_max = -5.0, 5.0
            grid.y_min, grid.y_max = -5.0, 5.0
            grid.z_min, grid.z_max = -5.0, 5.0

            gen = OrbitalsGenerator(
                nAtoms=nAtoms,
                atoms_R=atoms_R,
                spherical=1,
                nb=1,
                basis=basis,
                basis_norm=basis_norm,
                grid=grid,
            )

            gen.init_aos()
            gen.calc_aos(gen.AOs)

            center = res // 2
            values.append(gen.AOs[0, center, center, center])

        # All values should be close to expected (grid origin is exact)
        for i, val in enumerate(values):
            assert np.isclose(
                val, expected, rtol=1e-10
            ), f"Resolution {resolutions[i]}: expected {expected}, got {val}"


class TestAtomPositioning:
    """Test orbitals centered at non-origin positions."""

    def test_s_orbital_off_origin(self):
        """Test S-orbital centered away from coordinate origin."""
        # Place atom at (2, 0, 0)
        atom_pos = np.array([[2.0, 0.0, 0.0]])

        nAtoms = 1
        s_data = np.array([[1.0, 1.0]])
        basis = [[s_data]]
        basis_norm = [[[1.0]]]

        # Grid that includes both origin and atom position
        grid = Grid(x_n=41, y_n=21, z_n=21)
        grid.x_min, grid.x_max = -3.0, 7.0
        grid.y_min, grid.y_max = -5.0, 5.0
        grid.z_min, grid.z_max = -5.0, 5.0

        gen = OrbitalsGenerator(
            nAtoms=nAtoms,
            atoms_R=atom_pos,
            spherical=1,
            nb=1,
            basis=basis,
            basis_norm=basis_norm,
            grid=grid,
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        # Find grid index closest to atom position
        x_vals = np.linspace(grid.x_min, grid.x_max, grid.x_n)
        atom_x_idx = np.argmin(np.abs(x_vals - 2.0))
        center_y = grid.y_n // 2
        center_z = grid.z_n // 2

        # Value at atom position should be maximum
        value_at_atom = gen.AOs[0, atom_x_idx, center_y, center_z]
        expected = (2.0 / np.pi) ** 0.75

        # Allow some tolerance since atom might not be exactly on grid point
        assert np.isclose(value_at_atom, expected, rtol=0.01)

    def test_two_atoms(self):
        """Test system with two atoms at different positions."""
        # Two hydrogen-like atoms
        atoms_R = np.array(
            [
                [-1.5, 0.0, 0.0],
                [1.5, 0.0, 0.0],
            ]
        )

        nAtoms = 2
        nb = 2  # One S per atom

        s_data = np.array([[1.0, 1.0]])
        basis = [[s_data], [s_data]]  # Same basis on both atoms
        basis_norm = [[[1.0]], [[1.0]]]

        grid = Grid(x_n=41, y_n=21, z_n=21)
        grid.x_min, grid.x_max = -5.0, 5.0
        grid.y_min, grid.y_max = -4.0, 4.0
        grid.z_min, grid.z_max = -4.0, 4.0

        gen = OrbitalsGenerator(
            nAtoms=nAtoms,
            atoms_R=atoms_R,
            spherical=1,
            nb=nb,
            basis=basis,
            basis_norm=basis_norm,
            grid=grid,
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        # Should have 2 AOs
        assert gen.AOs.shape[0] == 2

        # Each AO should have its maximum near its atom
        x_vals = np.linspace(grid.x_min, grid.x_max, grid.x_n)
        center_y = grid.y_n // 2
        center_z = grid.z_n // 2

        # First AO centered at x=-1.5
        ao0_max_idx = np.argmax(gen.AOs[0, :, center_y, center_z])
        assert np.isclose(x_vals[ao0_max_idx], -1.5, atol=0.5)

        # Second AO centered at x=1.5
        ao1_max_idx = np.argmax(gen.AOs[1, :, center_y, center_z])
        assert np.isclose(x_vals[ao1_max_idx], 1.5, atol=0.5)


class TestEdgeCases:
    """Test edge cases and numerical stability."""

    def test_very_tight_gaussian(self):
        """Test Gaussian with large exponent (tight function)."""
        alpha = 100.0
        nAtoms = 1
        atoms_R = np.array([[0.0, 0.0, 0.0]])

        s_data = np.array([[alpha, 1.0]])
        basis = [[s_data]]
        basis_norm = [[[1.0]]]

        # Fine grid near origin (tight Gaussians are very localized)
        grid = Grid(x_n=41, y_n=41, z_n=41)
        grid.x_min, grid.x_max = -1.0, 1.0
        grid.y_min, grid.y_max = -1.0, 1.0
        grid.z_min, grid.z_max = -1.0, 1.0

        gen = OrbitalsGenerator(
            nAtoms=nAtoms,
            atoms_R=atoms_R,
            spherical=1,
            nb=1,
            basis=basis,
            basis_norm=basis_norm,
            grid=grid,
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        center = grid.x_n // 2
        expected = (2.0 * alpha / np.pi) ** 0.75
        actual = gen.AOs[0, center, center, center]

        assert np.isfinite(actual)
        assert np.isclose(actual, expected, rtol=1e-10)

    def test_diffuse_gaussian(self):
        """Test Gaussian with small exponent (diffuse function)."""
        alpha = 0.1
        nAtoms = 1
        atoms_R = np.array([[0.0, 0.0, 0.0]])

        s_data = np.array([[alpha, 1.0]])
        basis = [[s_data]]
        basis_norm = [[[1.0]]]

        # Large grid (diffuse Gaussians extend far)
        grid = Grid(x_n=21, y_n=21, z_n=21)
        grid.x_min, grid.x_max = -15.0, 15.0
        grid.y_min, grid.y_max = -15.0, 15.0
        grid.z_min, grid.z_max = -15.0, 15.0

        gen = OrbitalsGenerator(
            nAtoms=nAtoms,
            atoms_R=atoms_R,
            spherical=1,
            nb=1,
            basis=basis,
            basis_norm=basis_norm,
            grid=grid,
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        center = grid.x_n // 2
        expected = (2.0 * alpha / np.pi) ** 0.75
        actual = gen.AOs[0, center, center, center]

        assert np.isfinite(actual)
        assert np.isclose(actual, expected, rtol=1e-10)

    def test_contracted_gaussian(self):
        """Test contracted Gaussian (multiple primitives)."""
        # STO-3G like contraction
        nAtoms = 1
        atoms_R = np.array([[0.0, 0.0, 0.0]])

        # Multiple primitives for one contracted function
        s_data = np.array(
            [
                [3.42525091, 0.15432897],
                [0.62391373, 0.53532814],
                [0.16885540, 0.44463454],
            ]
        )
        basis = [[s_data]]
        basis_norm = [[[1.0]]]

        grid = Grid(x_n=21, y_n=21, z_n=21)
        grid.x_min, grid.x_max = -5.0, 5.0
        grid.y_min, grid.y_max = -5.0, 5.0
        grid.z_min, grid.z_max = -5.0, 5.0

        gen = OrbitalsGenerator(
            nAtoms=nAtoms,
            atoms_R=atoms_R,
            spherical=1,
            nb=1,
            basis=basis,
            basis_norm=basis_norm,
            grid=grid,
        )

        gen.init_aos()
        gen.calc_aos(gen.AOs)

        center = grid.x_n // 2
        value_at_origin = gen.AOs[0, center, center, center]

        # Should be sum of primitives
        expected = sum((2.0 * s_data[i, 0] / np.pi) ** 0.75 * s_data[i, 1] for i in range(3))

        assert np.isclose(value_at_origin, expected, rtol=1e-10)
        assert np.isfinite(value_at_origin)
        assert value_at_origin > 0
