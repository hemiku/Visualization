"""
Smoke tests - Basic sanity checks that imports work and basic functionality exists.

These tests ensure the package can be imported and basic objects can be created
without errors. They don't validate correctness, just that nothing crashes.
"""

import pytest
import os


class TestImports:
    """Test that all major modules can be imported."""

    def test_import_visualization(self):
        """Test importing main visualization module."""
        import visualization.visualization
        assert hasattr(visualization.visualization, 'Visualization')

    def test_import_molecular_system(self):
        """Test importing molecular system module."""
        import visualization.molecular_system
        assert hasattr(visualization.molecular_system, 'MolecularSystem')

    def test_import_orbitals(self):
        """Test importing orbitals module."""
        import visualization.orbitals
        assert hasattr(visualization.orbitals, 'OrbitalsGenerator')

    def test_import_grid(self):
        """Test importing grid module."""
        import visualization.grid
        assert hasattr(visualization.grid, 'Grid')

    def test_import_geminals(self):
        """Test importing geminals module."""
        import visualization.geminals
        assert hasattr(visualization.geminals, 'GeminalGenerator')

    def test_import_inputs(self):
        """Test importing inputs module."""
        import visualization.inputs
        assert hasattr(visualization.inputs, 'Input')

    def test_import_dalton_input(self):
        """Test importing Dalton parser."""
        import visualization.inputs
        assert hasattr(visualization.inputs, 'DaltonInput')

    def test_import_backends(self):
        """Test importing backends module."""
        import visualization.backends
        assert hasattr(visualization.backends, 'get_backend')
        assert hasattr(visualization.backends, 'set_backend')
        assert hasattr(visualization.backends, 'VisualizationBackend')

    def test_import_pyvista_backend(self):
        """Test importing PyVista backend (default)."""
        from visualization.backends.pyvista_backend import PyVistaBackend
        assert PyVistaBackend is not None


class TestBasicInstantiation:
    """Test that basic objects can be instantiated without crashing."""

    def test_create_molecular_system(self):
        """Test creating MolecularSystem object."""
        from visualization.molecular_system import MolecularSystem
        system = MolecularSystem()
        assert system is not None

    def test_create_grid(self):
        """Test creating Grid object with default parameters."""
        from visualization.grid import Grid
        grid = Grid(x_n=10, y_n=10, z_n=10)
        assert grid is not None
        assert grid.x_n == 10
        assert grid.y_n == 10
        assert grid.z_n == 10

    def test_create_visualization(self):
        """Test creating Visualization object."""
        from visualization.visualization import Visualization
        # Don't pass input_name to avoid file access
        viz = Visualization()
        assert viz is not None

    def test_create_pyvista_backend(self):
        """Test creating PyVista backend (default)."""
        from visualization.backends import get_backend
        backend = get_backend('pyvista')
        assert backend is not None
        assert backend.name == 'pyvista'

    def test_visualization_has_backend(self):
        """Test Visualization object has a backend."""
        from visualization.visualization import Visualization
        viz = Visualization()
        assert hasattr(viz, '_backend')
        assert viz._backend is not None
        assert viz._backend.name == 'pyvista'  # Default backend


class TestBackendFactory:
    """Test backend factory functionality."""

    def test_get_default_backend(self):
        """Test getting default backend (pyvista)."""
        from visualization.backends import get_backend, set_backend
        # Reset to ensure clean state
        set_backend('pyvista')
        backend = get_backend()
        assert backend.name == 'pyvista'

    def test_get_pyvista_backend_explicitly(self):
        """Test getting pyvista backend by name."""
        from visualization.backends import get_backend
        backend = get_backend('pyvista')
        assert backend.name == 'pyvista'

    def test_backend_env_variable(self):
        """Test backend selection via environment variable."""
        from visualization.backends import set_backend, get_backend
        # First reset to clean state
        set_backend('pyvista')
        # Set env var (this only affects new get_backend calls with name=None)
        old_env = os.environ.get('VIZ_BACKEND')
        try:
            os.environ['VIZ_BACKEND'] = 'pyvista'
            # Need to clear cache to pick up env var
            set_backend('pyvista')
            backend = get_backend()
            assert backend.name == 'pyvista'
        finally:
            if old_env is None:
                os.environ.pop('VIZ_BACKEND', None)
            else:
                os.environ['VIZ_BACKEND'] = old_env

    def test_invalid_backend_raises_error(self):
        """Test that invalid backend name raises ValueError."""
        from visualization.backends import get_backend
        with pytest.raises(ValueError, match="Unknown backend"):
            get_backend('nonexistent_backend')


class TestNumpyAvailability:
    """Test that NumPy is available and working."""

    def test_numpy_import(self):
        """Test NumPy can be imported."""
        import numpy as np
        assert np is not None

    def test_numpy_array_creation(self):
        """Test basic NumPy array operations."""
        import numpy as np
        arr = np.array([1, 2, 3])
        assert arr.shape == (3,)
        assert arr.sum() == 6


@pytest.mark.gpu
class TestGPUAvailability:
    """Test GPU/CuPy availability (skipped if not available)."""

    def test_cupy_import(self):
        """Test CuPy can be imported (requires GPU)."""
        import cupy as cp
        assert cp is not None

    def test_cupy_array_creation(self):
        """Test basic CuPy operations (requires GPU)."""
        import cupy as cp
        arr = cp.array([1, 2, 3])
        assert arr.shape == (3,)
