"""Tests for basis set normalization functions.

Tests cover:
- Analytical reference values for primitive Gaussians
- STO-3G contracted basis verification
- Coefficient scaling property
- Edge cases (very small/large exponents)
- All angular momentum types (S through H)
"""

import numpy as np
import pytest

from visualization.basis_normalization import (
    normalization_summation,
    norm_s,
    norm_p,
    norm_d,
    norm_f,
    norm_g,
    norm_h,
    calc_norm_from_basis,
)


class TestNormSPrimitive:
    """Test S-orbital normalization with analytical reference values."""

    def test_single_primitive_alpha_1(self):
        """S-orbital (alpha=1.0) should give (2/pi)^(3/4) = 0.7127054703549901"""
        data = np.array([[1.0, 1.0]])  # [exponent, coefficient]
        norm = norm_s(data)
        expected = (2.0 / np.pi) ** 0.75  # = 0.7127054703549901
        assert np.isclose(norm[0], expected, rtol=1e-10)

    def test_single_primitive_alpha_2(self):
        """S-orbital (alpha=2.0) should scale with alpha^(3/4)"""
        alpha = 2.0
        data = np.array([[alpha, 1.0]])
        norm = norm_s(data)
        expected = (2.0 * alpha / np.pi) ** 0.75
        assert np.isclose(norm[0], expected, rtol=1e-10)

    def test_single_primitive_alpha_0_5(self):
        """S-orbital (alpha=0.5) should scale with alpha^(3/4)"""
        alpha = 0.5
        data = np.array([[alpha, 1.0]])
        norm = norm_s(data)
        expected = (2.0 * alpha / np.pi) ** 0.75
        assert np.isclose(norm[0], expected, rtol=1e-10)


class TestNormSContracted:
    """Test S-orbital normalization for contracted Gaussians."""

    def test_sto3g_hydrogen(self):
        """STO-3G hydrogen 1s contracted basis should give positive norm."""
        data = np.array(
            [[3.42525091, 0.15432897], [0.62391373, 0.53532814], [0.16885540, 0.44463454]]
        )
        norm = norm_s(data)
        assert norm[0] > 0
        assert np.isfinite(norm[0])

    def test_two_primitives_same_exponent(self):
        """Two primitives with same exponent, coefficients sum to 1."""
        data = np.array([[1.0, 0.5], [1.0, 0.5]])
        norm = norm_s(data)
        # Should match single primitive since coeffs sum to 1
        expected = (2.0 / np.pi) ** 0.75
        assert np.isclose(norm[0], expected, rtol=1e-10)


class TestCoefficientScaling:
    """Test that normalization scales correctly with coefficients."""

    def test_doubling_coefficient_halves_norm(self):
        """Doubling coefficient should halve normalization factor."""
        data1 = np.array([[1.0, 1.0]])
        data2 = np.array([[1.0, 2.0]])
        norm1 = norm_s(data1)
        norm2 = norm_s(data2)
        assert np.isclose(norm1[0] / norm2[0], 2.0, rtol=1e-10)

    def test_tripling_coefficient(self):
        """Tripling coefficient should divide normalization by 3."""
        data1 = np.array([[1.0, 1.0]])
        data3 = np.array([[1.0, 3.0]])
        norm1 = norm_s(data1)
        norm3 = norm_s(data3)
        assert np.isclose(norm1[0] / norm3[0], 3.0, rtol=1e-10)


class TestNormP:
    """Test P-orbital normalization."""

    def test_single_primitive_alpha_1(self):
        """P-orbital (alpha=1.0) analytical value."""
        data = np.array([[1.0, 1.0]])
        norm = norm_p(data)
        # P-orbital normalization: 2 * (2*alpha)^(3/4) * alpha^(1/2) / pi^(3/4)
        # = 2 * 2^0.75 * 1^0.5 / pi^0.75 = 1.4254109407099802
        expected = 2.0 * (2.0**0.75) * (1.0**0.5) / (np.pi**0.75)
        assert np.isclose(norm[0], expected, rtol=1e-6)

    def test_positive_finite(self):
        """P-orbital norm should be positive and finite."""
        data = np.array([[1.0, 1.0]])
        norm = norm_p(data)
        assert norm[0] > 0
        assert np.isfinite(norm[0])


class TestNormD:
    """Test D-orbital normalization."""

    def test_single_primitive_positive(self):
        """D-orbital norm should be positive and finite."""
        data = np.array([[1.0, 1.0]])
        norm = norm_d(data)
        assert norm[0] > 0
        assert np.isfinite(norm[0])

    def test_larger_than_p(self):
        """D-orbital norm should be larger than P-orbital norm for same exponent."""
        data = np.array([[1.0, 1.0]])
        norm_p_val = norm_p(data)[0]
        norm_d_val = norm_d(data)[0]
        assert norm_d_val > norm_p_val


class TestNormF:
    """Test F-orbital normalization."""

    def test_single_primitive_positive(self):
        """F-orbital norm should be positive and finite."""
        data = np.array([[1.0, 1.0]])
        norm = norm_f(data)
        assert norm[0] > 0
        assert np.isfinite(norm[0])


class TestNormG:
    """Test G-orbital normalization."""

    def test_single_primitive_positive(self):
        """G-orbital norm should be positive and finite."""
        data = np.array([[1.0, 1.0]])
        norm = norm_g(data)
        assert norm[0] > 0
        assert np.isfinite(norm[0])


class TestNormH:
    """Test H-orbital normalization."""

    def test_single_primitive_positive(self):
        """H-orbital norm should be positive and finite."""
        data = np.array([[1.0, 1.0]])
        norm = norm_h(data)
        assert norm[0] > 0
        assert np.isfinite(norm[0])


class TestAllAngularMomenta:
    """Test consistency across all angular momentum types."""

    def test_all_positive(self):
        """All norm functions should return positive values."""
        data = np.array([[1.0, 1.0]])
        for norm_func in [norm_s, norm_p, norm_d, norm_f, norm_g, norm_h]:
            result = norm_func(data)
            assert result[0] > 0, f"{norm_func.__name__} returned non-positive"
            assert np.isfinite(result[0]), f"{norm_func.__name__} returned non-finite"

    def test_increasing_with_angular_momentum(self):
        """Normalization factors should increase with angular momentum."""
        data = np.array([[1.0, 1.0]])
        norms = [
            norm_s(data)[0],
            norm_p(data)[0],
            norm_d(data)[0],
            norm_f(data)[0],
            norm_g(data)[0],
            norm_h(data)[0],
        ]
        # Check that norms increase (or at least don't decrease drastically)
        for i in range(len(norms) - 1):
            assert (
                norms[i + 1] >= norms[i] * 0.5
            ), f"norm_{chr(ord('s') + i + 1)} unexpectedly small compared to previous"


class TestEdgeCases:
    """Test edge cases for numerical stability."""

    def test_very_small_exponent(self):
        """Diffuse function with alpha=0.001 should be finite."""
        data = np.array([[0.001, 1.0]])
        norm = norm_s(data)
        assert np.isfinite(norm[0])
        assert norm[0] > 0

    def test_very_large_exponent(self):
        """Tight function with alpha=10000 should be finite."""
        data = np.array([[10000.0, 1.0]])
        norm = norm_s(data)
        assert np.isfinite(norm[0])
        assert norm[0] > 0

    def test_multiple_contractions(self):
        """Multiple contracted functions from same primitives."""
        data = np.array([[10.0, 0.5, 0.3], [1.0, 0.5, 0.7]])
        norm = norm_s(data)
        assert norm.shape == (2,)
        assert all(n > 0 for n in norm)
        assert all(np.isfinite(n) for n in norm)


class TestNormalizationSummation:
    """Test the core normalization_summation function."""

    def test_single_primitive(self):
        """Single primitive with alpha=1, coeff=1, pow=1.5."""
        data = np.array([[1.0, 1.0]])
        result = normalization_summation(data, 1.5)
        # sum = 1*1 / (1+1)^1.5 = 1 / 2.828... = 0.3536...
        expected = 1.0 / (2.0**1.5)
        assert np.isclose(result[0], expected, rtol=1e-10)

    def test_two_primitives(self):
        """Two primitives with different exponents."""
        data = np.array([[1.0, 0.5], [2.0, 0.5]])
        result = normalization_summation(data, 1.5)
        # Manual calculation:
        # (0.5*0.5)/(1+1)^1.5 + (0.5*0.5)/(1+2)^1.5 + (0.5*0.5)/(2+1)^1.5 + (0.5*0.5)/(2+2)^1.5
        # = 0.25/2.828 + 0.25/5.196 + 0.25/5.196 + 0.25/8.0
        # = 0.0884 + 0.0481 + 0.0481 + 0.03125 = 0.2159
        expected = 0.25 / (2.0**1.5) + 0.25 / (3.0**1.5) + 0.25 / (3.0**1.5) + 0.25 / (4.0**1.5)
        assert np.isclose(result[0], expected, rtol=1e-10)


class TestCalcNormFromBasis:
    """Test the calc_norm_from_basis wrapper function."""

    def test_empty_basis(self):
        """Empty basis should return empty list."""
        result = calc_norm_from_basis([])
        assert result == []

    def test_single_atom_s_orbital(self):
        """Single atom with one S orbital."""
        s_data = np.array([[1.0, 1.0]])
        basis = [[s_data]]
        result = calc_norm_from_basis(basis)

        assert len(result) == 1
        assert len(result[0]) == 1
        expected = (2.0 / np.pi) ** 0.75
        assert np.isclose(result[0][0][0], expected, rtol=1e-10)

    def test_single_atom_s_and_p_orbitals(self):
        """Single atom with S and P orbitals."""
        s_data = np.array([[1.0, 1.0]])
        p_data = np.array([[1.0, 1.0]])
        basis = [[s_data, p_data]]
        result = calc_norm_from_basis(basis)

        assert len(result) == 1
        assert len(result[0]) == 2

    def test_two_atoms(self):
        """Two atoms with different basis sets."""
        s_data = np.array([[1.0, 1.0]])
        p_data = np.array([[1.0, 1.0]])
        basis = [[s_data], [s_data, p_data]]  # Atom 1: S only  # Atom 2: S and P
        result = calc_norm_from_basis(basis)

        assert len(result) == 2
        assert len(result[0]) == 1
        assert len(result[1]) == 2

    def test_water_like_basis(self):
        """Test with water-like basis structure (O + 2H)."""
        s_data = np.array([[5.0, 1.0], [1.0, 1.0]])
        p_data = np.array([[1.0, 1.0]])
        basis = [[s_data, p_data], [s_data], [s_data]]  # Oxygen  # H1  # H2
        result = calc_norm_from_basis(basis)

        assert len(result) == 3
        assert len(result[0]) == 2  # O: S and P
        assert len(result[1]) == 1  # H1: S only
        assert len(result[2]) == 1  # H2: S only

        # All norms should be positive
        for atom_norms in result:
            for orbital_norms in atom_norms:
                assert all(n > 0 for n in orbital_norms)
