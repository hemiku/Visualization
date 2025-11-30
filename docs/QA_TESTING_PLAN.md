# Comprehensive QA Testing Strategy

## Overview

Add thorough test coverage for critical untested components to catch silent bugs in scientific calculations. Scientific software with untested parsers can produce subtly wrong results that propagate through all downstream calculations.

**Approach:** QA (thorough) - Prioritizing correctness in scientific software

---

# Phase 1: Quick Wins and Foundation

## 1.1 Fix Regex Warnings in inputs.py

**File:** `visualization/inputs.py`

Add `r''` prefix to regex patterns to fix Python 3.12+ deprecation warnings:

| Line | Current Pattern | Fixed Pattern |
|------|-----------------|---------------|
| 273 | `' {1,9}\d{1,9}\. '` | `r' {1,9}\d{1,9}\. '` |
| 292 | `' {1,9}\d{1,9}\. '` | `r' {1,9}\d{1,9}\. '` |
| 293 | `'^\w{1,6} \D \D \D'` | `r'^\w{1,6} \D \D \D'` |
| 294 | `'^H {1,3}\d{1,3} {1,4}\d{1,3}$'` | `r'^H {1,3}\d{1,3} {1,4}\d{1,3}$'` |
| 368 | `' {1,9}\d{1,9}\. '` | `r' {1,9}\d{1,9}\. '` |
| 369 | `'^\w{1,6} \D \D \D'` | `r'^\w{1,6} \D \D \D'` |
| 370 | `'^H {1,3}\d{1,3} {1,4}\d{1,3}$'` | `r'^H {1,3}\d{1,3} {1,4}\d{1,3}$'` |

---

## 1.2 Basis Normalization Unit Tests

**File to create:** `tests/unit/test_basis_normalization.py`
**Target:** `visualization/basis_normalization.py`

### Mathematical Background

For Gaussian-type orbitals:
```
phi(r) = N * x^i * y^j * z^k * exp(-alpha * r^2)
```

The normalization condition requires `integral |phi(r)|^2 d^3r = 1`

### Analytical Reference Values

**Primitive Gaussian Normalization Constants:**

| Orbital Type | l | Formula | alpha=1.0 Value |
|--------------|---|---------|-----------------|
| S-type | 0 | `(2*alpha/pi)^(3/4)` | **0.7127054703549901** |
| P-type | 1 | `2*(2*alpha)^(3/4)*alpha^(1/2)/pi^(3/4)` | **1.4254109407099802** |
| D-type | 2 | `(2048*alpha^7/pi^3)^(1/4)` | **2.0170449244449186** |
| F-type | 3 | `(32768*alpha^9/pi^3)^(1/4)` | **2.8508218814199605** |

### Test Cases

```python
# Test 1: Single primitive S-orbital (alpha=1.0)
def test_norm_s_single_primitive():
    data = np.array([[1.0, 1.0]])  # [exponent, coefficient]
    norm = norm_s(data)
    expected = (2.0 / np.pi) ** 0.75  # = 0.7127054703549901
    assert np.isclose(norm[0], expected, rtol=1e-10)

# Test 2: STO-3G hydrogen 1s contracted basis
def test_norm_s_sto3g_hydrogen():
    data = np.array([
        [3.42525091, 0.15432897],
        [0.62391373, 0.53532814],
        [0.16885540, 0.44463454]
    ])
    norm = norm_s(data)
    assert norm[0] > 0 and np.isfinite(norm[0])

# Test 3: P-orbital normalization
def test_norm_p_single_primitive():
    data = np.array([[1.0, 1.0]])
    norm = norm_p(data)
    expected = 1.4254109407099802
    assert np.isclose(norm[0], expected, rtol=1e-6)

# Test 4: Coefficient scaling property
def test_coefficient_scaling():
    """Doubling coefficient should halve normalization factor"""
    data1 = np.array([[1.0, 1.0]])
    data2 = np.array([[1.0, 2.0]])
    norm1 = norm_s(data1)
    norm2 = norm_s(data2)
    assert np.isclose(norm1[0] / norm2[0], 2.0, rtol=1e-10)

# Test 5: Multiple contracted functions
def test_norm_s_multiple_contractions():
    data = np.array([
        [10.0, 0.5, 0.3],
        [1.0, 0.5, 0.7]
    ])
    norm = norm_s(data)
    assert norm.shape == (2,)
    assert all(n > 0 for n in norm)
```

### Edge Cases

```python
def test_very_small_exponent():
    """Diffuse function with alpha=0.001"""
    data = np.array([[0.001, 1.0]])
    norm = norm_s(data)
    assert np.isfinite(norm[0]) and norm[0] > 0

def test_very_large_exponent():
    """Tight function with alpha=10000"""
    data = np.array([[10000.0, 1.0]])
    norm = norm_s(data)
    assert np.isfinite(norm[0]) and norm[0] > 0

def test_all_angular_momenta_positive():
    """All norm functions return positive values"""
    data = np.array([[1.0, 1.0]])
    for norm_func in [norm_s, norm_p, norm_d, norm_f, norm_g, norm_h]:
        result = norm_func(data)
        assert result[0] > 0 and np.isfinite(result[0])
```

---

# Phase 2: Input Parser Tests

## 2.1 Available Test Fixtures

| File Path | Format | Description | Parser |
|-----------|--------|-------------|--------|
| `Examples/H2O_H2O.out` | Dalton output | Water dimer HF | `DaltonInput` |
| `Examples/H2O_H2O.tar.gz` | Dalton tarball | Water dimer complete | `DaltonInput` |
| `Examples/C2H4_C2H4/ethylene.out` | Dalton output | Ethylene dimer GVB | `DaltonInput` |
| `Examples/C2H4_C2H4/ethylene.tar.gz` | Dalton tarball | Ethylene complete | `DaltonInput` |
| `Examples/benzcyclo/benzcyclo.out` | Molpro output | Benzene-cyclohexane SAPT | `MolproInput` |

---

## 2.2 DaltonInput Parser Tests

**File to create:** `tests/unit/test_input_dalton.py`

### Expected Values for Water Dimer (H2O_H2O.out)

| Property | Expected Value |
|----------|---------------|
| `nAtoms` | 6 |
| `nb` (basis functions) | 82 |
| `spherical` | True (1) |
| `Atoms_Name` | `['O1', 'H1', 'H2', 'O2', 'H3', 'H4']` |
| `Atoms_Charge` | `[8, 1, 1, 8, 1, 1]` |
| Bond count | 4 |

### Atomic Coordinates (Bohr)

```python
expected_atoms_R = np.array([
    [0.0614560743, 0.0057782875, -0.1088907039],   # O1
    [-1.7447100105, -0.0153318415, 0.1667182727],  # H1
    [0.7691565673, -0.0763926677, 1.5618122632],   # H2
    [-5.4716162668, -0.0148424025, 0.3212978585],  # O2
    [-6.2145544803, -1.4081241458, -0.5815659756], # H3
    [-6.1912447103, 1.4641100685, -0.4554362145],  # H4
])
```

### Test Cases

```python
class TestDaltonInputParser:

    def test_get_nb(self, dalton_input_water):
        nb = dalton_input_water.get_nb()
        assert nb == 82

    def test_get_nAtoms(self, dalton_input_water):
        nAtoms = dalton_input_water.get_nAtoms()
        assert nAtoms == 6

    def test_get_spherical(self, dalton_input_water):
        spherical = dalton_input_water.get_spherical()
        assert spherical == 1

    def test_get_atoms_coordinates(self, dalton_input_water):
        atoms_R, _, _ = dalton_input_water.get_atoms()
        np.testing.assert_allclose(
            atoms_R[0],
            [0.0614560743, 0.0057782875, -0.1088907039],
            rtol=1e-6
        )

    def test_get_atoms_charges(self, dalton_input_water):
        _, atoms_Charge, _ = dalton_input_water.get_atoms()
        expected = np.array([8, 1, 1, 8, 1, 1])
        np.testing.assert_array_equal(atoms_Charge, expected)

    def test_get_atoms_names(self, dalton_input_water):
        _, _, atoms_Name = dalton_input_water.get_atoms()
        expected = ['O1', 'H1', 'H2', 'O2', 'H3', 'H4']
        assert atoms_Name == expected

    def test_get_bonds(self, dalton_input_water):
        bonds = dalton_input_water.get_bonds()
        assert len(bonds) == 4

    def test_get_coeff_shape(self, dalton_input_tar):
        nb = dalton_input_tar.get_nb()
        coeff = dalton_input_tar.get_coeff()
        assert coeff.shape == (nb, nb)
```

### Ethylene Dimer Expected Values (GVB/APSG)

| Property | Expected Value |
|----------|---------------|
| `nAtoms` | 12 |
| `nb` | 96 |
| `inactive` | 4 |
| `electrons` | 24 |

```python
def test_get_inactive_ethylene(self, dalton_input_ethylene):
    inactive = dalton_input_ethylene.get_inactive()
    assert inactive == 4

def test_get_electrons_ethylene(self, dalton_input_ethylene):
    electrons = dalton_input_ethylene.get_electrons()
    assert electrons == 24

def test_get_g_coeff_ethylene(self, dalton_input_ethylene):
    g_coeff = dalton_input_ethylene.get_g_coeff()
    np.testing.assert_allclose(g_coeff[0], 0.996475299688, rtol=1e-6)
```

---

## 2.3 MolproInput Parser Tests

**File to create:** `tests/unit/test_input_molpro.py`

### Expected Values for Benzene-Cyclohexane

| Property | Expected Value |
|----------|---------------|
| `nAtoms` | 27 |
| `nb` | 874 |
| `spherical` | True |
| `electrons` (monomer A) | 42 |
| `electrons` (monomer B) | 40 |

### Test Cases

```python
class TestMolproInputParser:

    def test_get_nb(self, molpro_input_benzcyclo):
        nb = molpro_input_benzcyclo.get_nb()
        assert nb == 874

    def test_get_nAtoms(self, molpro_input_benzcyclo):
        nAtoms = molpro_input_benzcyclo.get_nAtoms()
        assert nAtoms == 27

    def test_get_spherical(self, molpro_input_benzcyclo):
        spherical = molpro_input_benzcyclo.get_spherical()
        assert spherical == 1

    def test_get_atoms_first_coordinate(self, molpro_input_benzcyclo):
        atoms_R, _, _ = molpro_input_benzcyclo.get_atoms()
        np.testing.assert_allclose(
            atoms_R[0],
            [1.446671259, 1.640743993, 1.551448056],
            rtol=1e-6
        )


class TestMolproSaptInput:

    def test_different_monomer_charges(self, molpro_sapt_input):
        molpro_sapt_input.monomer = 0
        electrons_A = molpro_sapt_input.get_electrons()

        molpro_sapt_input.monomer = 1
        molpro_sapt_input.output = None
        electrons_B = molpro_sapt_input.get_electrons()

        assert electrons_A == 42
        assert electrons_B == 40
```

---

## 2.4 Error Handling Tests

```python
class TestParserErrorHandling:

    def test_dalton_missing_file(self):
        dalton = DaltonInput(input_name='/nonexistent/path')
        with pytest.raises(FileNotFoundError):
            dalton.get_dalton_output()

    def test_molpro_missing_file(self):
        molpro = MolproInput(input_name='/nonexistent/path')
        with pytest.raises(FileNotFoundError):
            molpro.get_output()

    def test_invalid_input_type(self):
        with pytest.raises(Exception, match="Input_type: .* not found"):
            get_input(input_type='InvalidType', input_sub_type='test')
```

---

# Phase 3: Known-Answer Validation Tests

**File to create:** `tests/validation/test_ao_known_answer.py`

## 3.1 Atomic Orbital Validation

### S-orbital at Origin

For a normalized 1s Gaussian with alpha=1.0 at r=0:
```
phi_1s(0) = (2*alpha/pi)^(3/4) = 0.7127054703549901
```

```python
def test_s_orbital_value_at_origin():
    """S-orbital at nucleus should equal (2*alpha/pi)^(3/4)"""
    grid = Grid(x_n=11, y_n=11, z_n=11)
    grid.x_min, grid.x_max = -5.0, 5.0
    grid.y_min, grid.y_max = -5.0, 5.0
    grid.z_min, grid.z_max = -5.0, 5.0

    setup = {
        'nAtoms': 1,
        'atoms_R': np.array([[0.0, 0.0, 0.0]]),
        'spherical': False,
        'nb': 1,
        'basis': [[np.array([[1.0, 1.0]])]],
        'basis_norm': [[np.array([1.0])]],
    }

    og = OrbitalsGenerator(grid=grid, **setup)
    og.init_aos()
    og.calc_aos(og.AOs, grid)

    ao_at_origin = og.AOs[0, 5, 5, 5]
    expected = (2.0 / np.pi) ** 0.75

    assert np.isclose(ao_at_origin, expected, rtol=0.01)
```

### S-orbital Normalization Integral

```python
def test_s_orbital_normalization_integral():
    """Integral of |phi|^2 dV should be approximately 1"""
    grid = Grid(x_n=51, y_n=51, z_n=51)
    grid.x_min, grid.x_max = -8.0, 8.0
    grid.y_min, grid.y_max = -8.0, 8.0
    grid.z_min, grid.z_max = -8.0, 8.0

    # ... setup OrbitalsGenerator ...

    dx = 16.0 / 50
    dv = dx ** 3
    integral = np.sum(og.AOs[0]**2) * dv

    assert 0.9 < integral < 1.1
```

### P-orbital Antisymmetry

```python
def test_p_orbital_antisymmetry():
    """P_x orbital should be antisymmetric: phi(-x) = -phi(x)"""
    # ... setup with P orbitals ...

    center = 5
    val_plus = og.AOs[1, center+2, center, center]
    val_minus = og.AOs[1, center-2, center, center]

    assert np.isclose(val_plus, -val_minus, rtol=0.05)
```

### Grid Convergence

```python
def test_normalization_converges():
    """Normalization integral should converge with grid density."""
    integrals = []
    for n in [21, 41, 61]:
        # ... compute integral at each resolution ...
        integrals.append(integral)

    # Should be getting closer to 1.0
    assert abs(integrals[-1] - 1.0) < abs(integrals[0] - 1.0)
```

---

# Phase 4: CI/CD Improvements

**File:** `.github/workflows/tests.yml`

## 4.1 Add Python 3.12 to Test Matrix

```yaml
strategy:
  matrix:
    python-version: ["3.8", "3.9", "3.10", "3.11", "3.12"]
```

## 4.2 Make Code Quality Checks Hard-Fail

Remove `continue-on-error: true` from black and flake8 steps.

## 4.3 Add SyntaxWarning Check

```yaml
- name: Check for SyntaxWarnings
  run: |
    python -W error::SyntaxWarning -c "import visualization"
```

---

# Pytest Fixtures to Add

**File:** `tests/conftest.py`

```python
@pytest.fixture
def dalton_input_water():
    """DaltonInput for water dimer output file."""
    from visualization.inputs import DaltonInput
    return DaltonInput(
        input_type='Output',
        input_name='Examples/H2O_H2O'
    )

@pytest.fixture
def dalton_input_tar():
    """DaltonInput for water dimer tar archive."""
    from visualization.inputs import DaltonInput
    return DaltonInput(
        input_type='tar',
        input_name='Examples/H2O_H2O'
    )

@pytest.fixture
def dalton_input_ethylene():
    """DaltonInput for ethylene dimer."""
    from visualization.inputs import DaltonInput
    return DaltonInput(
        input_type='Output',
        input_name='Examples/C2H4_C2H4/ethylene'
    )

@pytest.fixture
def molpro_input_benzcyclo():
    """MolproInput for benzene-cyclohexane."""
    from visualization.input_molpro import MolproInput
    return MolproInput(input_name='Examples/benzcyclo/benzcyclo')

@pytest.fixture
def primitive_s_orbital():
    """Single primitive S-orbital with alpha=1.0"""
    return np.array([[1.0, 1.0]])

@pytest.fixture
def sto3g_hydrogen():
    """STO-3G hydrogen basis"""
    return np.array([
        [3.42525091, 0.15432897],
        [0.62391373, 0.53532814],
        [0.16885540, 0.44463454]
    ])

@pytest.fixture
def expected_water_dimer_atoms():
    """Expected atomic data for water dimer."""
    return {
        'nAtoms': 6,
        'names': ['O1', 'H1', 'H2', 'O2', 'H3', 'H4'],
        'charges': np.array([8, 1, 1, 8, 1, 1]),
        'positions': np.array([
            [0.0614560743, 0.0057782875, -0.1088907039],
            [-1.7447100105, -0.0153318415, 0.1667182727],
            [0.7691565673, -0.0763926677, 1.5618122632],
            [-5.4716162668, -0.0148424025, 0.3212978585],
            [-6.2145544803, -1.4081241458, -0.5815659756],
            [-6.1912447103, 1.4641100685, -0.4554362145],
        ])
    }
```

---

# Files Summary

## Files to Create

| File | Description |
|------|-------------|
| `tests/unit/test_basis_normalization.py` | Normalization function tests with analytical values |
| `tests/unit/test_input_dalton.py` | DaltonInput parser tests |
| `tests/unit/test_input_molpro.py` | MolproInput/MolproSaptInput tests |
| `tests/validation/test_ao_known_answer.py` | Atomic orbital validation |

## Files to Modify

| File | Changes |
|------|---------|
| `visualization/inputs.py` | Fix 7 regex patterns (add `r''` prefix) |
| `tests/conftest.py` | Add parser and normalization fixtures |
| `.github/workflows/tests.yml` | Add Python 3.12, remove soft-fail |

---

# Success Criteria

- [ ] All 7 regex warnings resolved in inputs.py
- [ ] basis_normalization.py has >90% test coverage
- [ ] Normalization tests use analytical reference values
- [ ] DaltonInput parser tests with known expected values
- [ ] MolproInput parser tests with known expected values
- [ ] At least one known-answer validation test per orbital type (S, P)
- [ ] CI passes on Python 3.8-3.12
- [ ] Code quality checks are hard-fail in CI

---

# Estimated Effort

| Phase | Time |
|-------|------|
| Phase 1: Regex fixes + Normalization tests | 2-3 hours |
| Phase 2: Parser tests | 4-6 hours |
| Phase 3: Known-answer validation | 3-4 hours |
| Phase 4: CI/CD improvements | 1 hour |
| **Total** | **2-3 days** |
