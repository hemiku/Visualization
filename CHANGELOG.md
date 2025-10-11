# Changelog

All notable changes to the GVT (General Visualization Tool) project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

#### Package Management & Environment
- UV package manager support with uv.lock for reproducible builds (10-100x faster than pip)
- Hybrid Conda+UV setup (Conda for Mayavi, UV for Python packages)
- Automated environment setup script (`setup_environment.sh`)
- Modern `pyproject.toml` with PEP 621 compliance and hatchling backend
- `.python-version` file for Python 3.12
- Comprehensive setup documentation (`UV_SETUP.md`)

#### Testing Infrastructure
- Comprehensive test suite with 34 tests across 5 categories:
  - **Unit tests** (14 tests): Orbital calculation, grid initialization
  - **Integration tests** (11 tests): Dispersion, geminal, and SAPT workflows
  - **Validation tests** (14 tests): Scientific correctness (normalization, orthogonality)
  - **Visual tests** (4 tests): Rendering validation with reference images
  - **Smoke tests** (11 tests): Basic import and instantiation checks
- `pytest.ini` configuration with custom test markers
- Shared test fixtures in `tests/conftest.py`
- Reference image for visual regression testing
- Testing documentation (`TESTING.md`)

#### CI/CD
- GitHub Actions workflow (`.github/workflows/tests.yml`)
- Automated testing on push and pull requests
- Fast dependency installation with UV in CI

#### Scripts & Tools
- `show_water.py`: Manual verification script with interactive and headless modes
  - Interactive mode: Opens 3D Mayavi window
  - Headless mode (`--save`): Renders to PNG file for automated testing

#### Documentation
- Comprehensive README with installation, examples, and testing instructions
- Testing guide (`TESTING.md`) with detailed test descriptions
- UV setup guide (`UV_SETUP.md`) with multiple installation approaches
- CHANGELOG for tracking project changes

### Changed

#### Packaging
- Migrated from `setup.py` to modern `pyproject.toml` (PEP 621)
- Changed build backend from setuptools to hatchling (faster builds)
- Reorganized dependency specifications:
  - Core runtime: numpy, scipy
  - Optional GPU: cupy
  - Dev tools: pytest, black, flake8, mypy

#### Testing
- Refactored all tests to use `pathlib` instead of string concatenation
- Improved test fixtures to return Path objects
- Fixed path handling to prevent `.out.out` file extension bugs
- Updated test data paths to use consistent patterns

#### Configuration
- Enhanced `.gitignore` with UV cache, testing artifacts, and temporary files
- Added pytest markers for better test categorization

### Fixed
- Path handling in SAPT tests (removed `.out` extension doubling)
- Test fixtures now properly use tarball data for complete Dalton inputs
- Grid array initialization with molecular systems

### Development

#### Test Organization
```
tests/
├── unit/               # Fast, isolated component tests
├── integration/        # End-to-end workflow tests
├── validation/         # Scientific correctness checks
├── visual/            # Rendering and output validation
└── conftest.py        # Shared fixtures and configuration
```

#### Performance
- Small test grids (5x5x5, 10x10x10) for fast feedback
- Slow tests marked with `@pytest.mark.slow`
- Test suite completes in ~25 seconds (excluding slow tests)

## [0.0.1] - Previous Release

### Features
- Molecular geometry visualization with Mayavi
- Atomic orbital (AO) and molecular orbital (MO) visualization
- Dispersion interaction calculation from EERPA-GVB
- Geminal wavefunction visualization
- SAPT(CAS) dispersion component visualization
- GPU acceleration support via CuPy
- Support for Dalton, Molpro, and GAMMCOR input formats

---

## Version History Notes

This CHANGELOG was introduced as part of the test suite development and UV migration.
Previous changes were not formally tracked.
