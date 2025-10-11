"""
Visual validation tests for water dimer visualization.

Tests verify that:
1. Images render without errors
2. Output has expected properties (size, content, non-blank)
3. Molecular data is loaded correctly
4. (Optional) Visual similarity to reference image
"""

import pytest
import numpy as np
from pathlib import Path


@pytest.mark.visual
@pytest.mark.requires_examples
class TestWaterVisualization:
    """Test water dimer visualization output."""

    def test_render_water_geometry(self, tmp_path, water_dimer_tarball):
        """Test that water dimer renders without errors."""
        import visualization.visualization as V

        # Load water dimer
        input_path = Path(water_dimer_tarball).with_suffix('').with_suffix('')

        vis = V.Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(input_path)
        )

        # Get geometry
        vis.get_geometry()

        # Render to file (headless)
        output_file = tmp_path / "test_water.png"

        figure = vis.plot_Geometry(
            plot_atoms=True,
            atom_scaling=0.5,
            atom_names=True,
            plot_bonds=True,
            auto_show=False
        )

        vis.mlab.savefig(str(output_file))
        vis.mlab.close()

        # Verify file was created
        assert output_file.exists(), "Output image was not created"

    def test_image_properties(self, tmp_path, water_dimer_tarball):
        """Test that rendered image has expected properties."""
        import visualization.visualization as V
        from PIL import Image

        # Load and render
        input_path = Path(water_dimer_tarball).with_suffix('').with_suffix('')

        vis = V.Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(input_path)
        )

        vis.get_geometry()

        output_file = tmp_path / "test_water.png"

        figure = vis.plot_Geometry(
            plot_atoms=True,
            atom_scaling=0.5,
            atom_names=True,
            plot_bonds=True,
            auto_show=False
        )

        vis.mlab.savefig(str(output_file))
        vis.mlab.close()

        # Load image
        img = Image.open(output_file)
        img_array = np.array(img)

        # Test 1: Image has reasonable dimensions (around 600 pixels)
        width, height = img.size
        assert 500 <= width <= 700, f"Unexpected width: {width} (expected 500-700)"
        assert 500 <= height <= 700, f"Unexpected height: {height} (expected 500-700)"

        # Test 2: Image is not blank (has variation)
        std_dev = img_array.std()
        assert std_dev > 10, f"Image appears blank (std={std_dev:.2f}, expected >10)"

        # Test 3: Image is not all black (allow mostly white backgrounds)
        mean_color = img_array.mean()
        assert mean_color > 30, \
            f"Image appears all black: mean={mean_color:.1f} (expected >30)"

        # Test 4: Image has multiple colors (not monochrome)
        # Check RGB channels have different values
        if len(img_array.shape) == 3 and img_array.shape[2] >= 3:
            r_mean = img_array[:, :, 0].mean()
            g_mean = img_array[:, :, 1].mean()
            b_mean = img_array[:, :, 2].mean()

            # At least some color variation between channels
            # (relaxed threshold for white backgrounds with colored atoms)
            color_variation = max(r_mean, g_mean, b_mean) - min(r_mean, g_mean, b_mean)
            assert color_variation > 1, \
                f"Image appears monochrome (variation={color_variation:.1f})"

    def test_molecular_data_loaded(self, water_dimer_tarball):
        """Test that molecular system data is loaded correctly."""
        import visualization.visualization as V

        input_path = Path(water_dimer_tarball).with_suffix('').with_suffix('')

        vis = V.Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(input_path)
        )

        vis.get_geometry()

        # Test molecular system properties
        assert vis.molecular_system.nAtoms == 6, \
            f"Expected 6 atoms, got {vis.molecular_system.nAtoms}"

        assert vis.molecular_system.atoms_R is not None
        assert vis.molecular_system.atoms_R.shape == (6, 3), \
            "Atom positions should be 6x3 array"

        assert vis.molecular_system.atoms_Charge is not None
        assert len(vis.molecular_system.atoms_Charge) == 6

        assert vis.molecular_system.atoms_Name is not None
        assert len(vis.molecular_system.atoms_Name) == 6

        # Should have bonds
        assert vis.molecular_system.bonds is not None
        assert len(vis.molecular_system.bonds) > 0, "Should have at least one bond"

        # Verify atom types (2 oxygen, 4 hydrogen)
        charges = vis.molecular_system.atoms_Charge
        oxygen_count = sum(1 for c in charges if c == 8)
        hydrogen_count = sum(1 for c in charges if c == 1)

        assert oxygen_count == 2, f"Expected 2 oxygen atoms, got {oxygen_count}"
        assert hydrogen_count == 4, f"Expected 4 hydrogen atoms, got {hydrogen_count}"


@pytest.mark.visual
@pytest.mark.requires_examples
@pytest.mark.requires_reference
class TestVisualSimilarity:
    """Test visual similarity to reference images (optional)."""

    def test_compare_to_reference(self, tmp_path, water_dimer_tarball):
        """Compare rendered image to reference image."""
        pytest.importorskip("skimage", reason="scikit-image required for similarity testing")

        import visualization.visualization as V
        from PIL import Image
        from skimage.metrics import structural_similarity as ssim

        # Check if reference exists
        reference_path = Path("Examples/reference_outputs/water_dimer_geometry.png")
        if not reference_path.exists():
            pytest.skip(f"Reference image not found: {reference_path}")

        # Load and render
        input_path = Path(water_dimer_tarball).with_suffix('').with_suffix('')

        vis = V.Visualization(
            input_type='Dalton',
            input_sub_type='tar',
            input_name=str(input_path)
        )

        vis.get_geometry()

        output_file = tmp_path / "test_water.png"

        figure = vis.plot_Geometry(
            plot_atoms=True,
            atom_scaling=0.5,
            atom_names=True,
            plot_bonds=True,
            auto_show=False
        )

        vis.mlab.savefig(str(output_file))
        vis.mlab.close()

        # Load images
        reference = np.array(Image.open(reference_path))
        test_output = np.array(Image.open(output_file))

        # Ensure same shape
        assert reference.shape == test_output.shape, \
            f"Shape mismatch: reference={reference.shape}, output={test_output.shape}"

        # Calculate structural similarity
        # Use channel_axis=2 for RGB images, channel_axis=-1 for RGBA
        if len(reference.shape) == 3:
            similarity = ssim(reference, test_output, channel_axis=2)
        else:
            similarity = ssim(reference, test_output)

        # Allow some variation due to rendering differences
        # SSIM of 1.0 = identical, 0.0 = completely different
        min_similarity = 0.90

        assert similarity >= min_similarity, \
            f"Image too different from reference: SSIM={similarity:.3f} (expected >={min_similarity})"

        print(f"Visual similarity: {similarity:.3f}")
