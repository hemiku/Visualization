# Testing Guide for GVT Visualization

## Overview

This project includes comprehensive testing for the visualization package:
- **Unit tests**: Individual component testing
- **Integration tests**: Workflow testing (dispersion, geminals, SAPT)
- **Validation tests**: Scientific correctness (orbital properties)
- **Visual tests**: Rendering and output validation

## Visual Validation Tests

### Purpose
Visual tests ensure that:
1. Rendering completes without errors
2. Output images have correct properties
3. Molecular data loads correctly
4. (Optional) Visual similarity to reference images

### Running Visual Tests

**Quick property-based tests (fast):**
```bash
pytest tests/visual/test_water_visualization.py::TestWaterVisualization -v
```

**All visual tests:**
```bash
pytest -m visual -v
```

**With reference image comparison (requires scikit-image):**
```bash
pytest -m requires_reference -v
```

### Manual Verification Script

`show_water.py` - Simple script to visualize water dimer

**Interactive mode (opens 3D window):**
```bash
python show_water.py
```

**Headless mode (save to file):**
```bash
python show_water.py --save
```

## Test Categories

### Property-Based Tests (Always Run)
✅ `test_render_water_geometry` - Renders without errors
✅ `test_image_properties` - Image size, content, color validation
✅ `test_molecular_data_loaded` - Atom count, bonds, charges

These tests run fast and don't require reference images.

### Reference Comparison Tests (Optional)
⚠️ `test_compare_to_reference` - SSIM similarity check

Requires:
- `scikit-image` package
- Reference image at `Examples/reference_outputs/water_dimer_geometry.png`

## Test Markers

Use pytest markers to run specific test categories:

```bash
# Visual tests only
pytest -m visual

# Excluding slow tests
pytest -m "not slow"

# Integration + validation tests
pytest -m "integration or validation"

# Visual tests requiring examples
pytest -m "visual and requires_examples"
```

## Files Created

### Test Files
- `tests/visual/test_water_visualization.py` - Visual validation tests
- `tests/visual/__init__.py` - Package init

### Reference Images
- `Examples/reference_outputs/water_dimer_geometry.png` - Known-good reference

### Scripts
- `show_water.py` - Manual verification script with `--save` option

### Configuration
- `pytest.ini` - Added `visual` and `requires_reference` markers

## What the Tests Validate

### 1. Rendering Tests
- Mayavi can create figures without errors
- `mlab.savefig()` works correctly
- Files are created with expected format

### 2. Image Property Tests
- **Dimensions**: Image size is reasonable (500-700px)
- **Content**: Image is not blank (std dev > 10)
- **Not all black**: Mean color > 30
- **Has colors**: RGB channel variation > 1

### 3. Molecular Data Tests
- **Atom count**: 6 atoms (2 O + 4 H)
- **Positions**: 6x3 array of coordinates
- **Charges**: Correct atomic numbers (8 for O, 1 for H)
- **Bonds**: At least one bond detected
- **Names**: Atom labels present

### 4. Reference Comparison (Optional)
- **SSIM similarity**: >= 0.90 (90% similar to reference)
- Allows small rendering variations
- Detects major visual regressions

## CI/CD Integration

Visual tests can run in CI/CD with headless rendering:

```yaml
# Example GitHub Actions
- name: Run visual tests
  run: |
    pytest tests/visual/ -v --tb=short
  env:
    QT_QPA_PLATFORM: offscreen
```

## Troubleshooting

### "Image appears blank" error
- Check Mayavi installation
- Verify `mlab.savefig()` permissions
- Try interactive mode first: `python show_water.py`

### "Reference image not found" warning
- Run: `python show_water.py --save`
- Copy to: `Examples/reference_outputs/water_dimer_geometry.png`

### SSIM test failures
- Small differences are normal (rendering variations)
- Adjust `min_similarity` threshold if needed
- Regenerate reference image if visualization improved

## Adding New Visual Tests

1. **Create test function** in `tests/visual/test_*.py`
2. **Add markers**: `@pytest.mark.visual`, `@pytest.mark.requires_examples`
3. **Use tmp_path** for output files
4. **Close figures**: Always call `vis.mlab.close()`
5. **Add assertions** for expected properties

Example:
```python
@pytest.mark.visual
@pytest.mark.requires_examples
def test_my_visualization(tmp_path):
    # Setup
    vis = create_visualization()

    # Render
    output = tmp_path / "output.png"
    fig = vis.plot_something(auto_show=False)
    vis.mlab.savefig(str(output))
    vis.mlab.close()

    # Validate
    assert output.exists()
    img = Image.open(output)
    assert img.size[0] > 100  # Reasonable size
```

## Future Enhancements

- [ ] Add tests for dispersion visualization
- [ ] Add tests for geminal visualization
- [ ] Add tests for SAPT visualization
- [ ] Generate reference images for all workflows
- [ ] Add performance benchmarks
- [ ] Test different rendering backends
- [ ] Add 3D coordinate validation tests
