"""
Simple script to visualize the water dimer (H2O_H2O) using GVT tools.

This uses the built-in plot_Geometry() method to show molecular structure.
Run this to verify that the visualization package is working correctly.

Usage:
    python show_water.py           # Interactive mode (opens window)
    python show_water.py --save    # Headless mode (save to file)
"""

# %% Import GVT modules
import visualization.visualization as V
from pathlib import Path
import sys

# %% Configuration
SAVE_MODE = '--save' in sys.argv
OUTPUT_FILE = 'water_dimer_output.png'

# %% Load water dimer
print("Loading water dimer data...")
water_path = Path("Examples/H2O_H2O")

vis = V.Visualization(
    input_type='Dalton',
    input_sub_type='tar',
    input_name=str(water_path)
)

# %% Get geometry
print("Loading geometry...")
vis.get_geometry()
print(f"✓ Loaded {vis.molecular_system.nAtoms} atoms")
print(f"  Atom types: {', '.join(set(vis.molecular_system.atoms_Name))}")

# %% Visualize with GVT plot_Geometry method
if SAVE_MODE:
    print(f"\nRendering to {OUTPUT_FILE}...")

    # Headless rendering
    figure = vis.plot_Geometry(
        plot_atoms=True,
        atom_scaling=0.5,
        atom_names=True,
        atom_names_scaling=1.0,
        plot_bonds=True,
        bond_scaling=1.0,
        auto_show=False  # Don't open window
    )

    # Save to file
    vis.mlab.savefig(OUTPUT_FILE)
    vis.mlab.close()

    print(f"✅ Saved to {OUTPUT_FILE}")

else:
    print("\nOpening 3D visualization...")
    print("(Close the Mayavi window when done)\n")

    vis.plot_Geometry(
        plot_atoms=True,
        atom_scaling=0.5,
        atom_names=True,
        atom_names_scaling=1.0,
        plot_bonds=True,
        bond_scaling=1.0,
        auto_show=True  # Open window
    )

    print("✅ Visualization complete!")
    print("💡 You can rotate the molecule with your mouse")
