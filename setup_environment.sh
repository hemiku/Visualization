#!/bin/bash
# Setup script for Visualization project with Conda + UV

set -e  # Exit on error

echo "🚀 Setting up Visualization project environment..."
echo ""

# Check if conda is available
if ! command -v conda &> /dev/null; then
    echo "❌ Error: conda not found. Please install Anaconda or Miniconda first."
    exit 1
fi

# Environment name
ENV_NAME="viz-env"

echo "📦 Creating conda environment: $ENV_NAME"
echo "   This will install Python 3.12 and Mayavi (may take 5-10 minutes)..."
conda create -n $ENV_NAME python=3.12 mayavi -c conda-forge -y

echo ""
echo "✅ Conda environment created!"
echo ""
echo "🔧 Activating environment and installing UV..."

# Activate environment and install UV
eval "$(conda shell.bash hook)"
conda activate $ENV_NAME

# Install UV
pip install uv

echo ""
echo "📥 Installing project with UV..."
uv pip install -e .

echo ""
echo "🧪 Installing development dependencies..."
uv pip install -e ".[dev]"

echo ""
echo "✨ Testing installation..."
python -c "from mayavi import mlab; print('✅ Mayavi OK')"
python -c "import numpy; print('✅ NumPy OK')"
python -c "import scipy; print('✅ SciPy OK')"

echo ""
echo "🎉 Setup complete!"
echo ""
echo "To activate the environment, run:"
echo "   conda activate $ENV_NAME"
echo ""
echo "To run tests:"
echo "   pytest tests/unit/test_smoke.py -v"
echo ""
echo "To deactivate when done:"
echo "   conda deactivate"
