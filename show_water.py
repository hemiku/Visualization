"""
Simple script to visualize the water dimer (H2O_H2O) using GVT tools.

This script demonstrates:
1. Molecular geometry visualization
2. Orbital calculation (AO and MO)
3. Molecular orbital visualization

Run this to verify that the visualization package is working correctly.

Usage:
    python show_water.py                # Interactive mode (opens window)
    python show_water.py --save         # Headless mode (save to file)
    python show_water.py --orbitals     # Include orbital calculations
"""

# %% Import GVT modules
import visualization.visualization as V
from pathlib import Path
import sys

# %% Configuration
SAVE_MODE = '--save' in sys.argv
CALC_ORBITALS = '--orbitals' in sys.argv
OUTPUT_FILE = 'water_dimer_output.png'
ORBITAL_OUTPUT_FILE = 'water_dimer_orbitals.png'

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

# %% Orbital calculations (optional)
if CALC_ORBITALS:
    print("\n" + "="*60)
    print("ORBITAL CALCULATIONS")
    print("="*60)

    # Calculate orbitals with small grid for demonstration
    print("\nCalculating orbitals (this may take a minute)...")
    print("Grid size: 20x20x20 (small for demonstration)")

    try:
        vis.get_orbitals(
            x_n=20,
            y_n=20,
            z_n=20,
            R_max_multip=3.0,
            gpu=False
        )

        print(f"✓ Calculated {len(vis.molecular_system.AOs)} atomic orbitals (AOs)")
        print(f"✓ Calculated {len(vis.molecular_system.MOs)} molecular orbitals (MOs)")
    except ValueError as e:
        print(f"\n⚠️  Orbital calculation failed: {e}")
        print("\nNote: The H2O_H2O tarball has a known issue with coefficient parsing.")
        print("Try using Examples/C2H4_C2H4/ethylene instead for orbital calculations:")
        print("  vis = V.Visualization(input_type='Dalton', input_sub_type='tar',")
        print("                        input_name='Examples/C2H4_C2H4/ethylene')")
        CALC_ORBITALS = False  # Skip visualization

if CALC_ORBITALS:
    # Visualize HOMO (Highest Occupied Molecular Orbital)
    # For water dimer, typically orbital around index 4-5 is interesting
    homo_index = 4  # Example - adjust based on actual system

    print(f"\nVisualizing MO #{homo_index}...")

    if SAVE_MODE:
        print(f"Rendering to {ORBITAL_OUTPUT_FILE}...")

        # Headless rendering
        figure = vis.plot_orbitals_MO(
            orbital_numbers=[homo_index],
            plot_atoms=True,
            atom_scaling=0.5,
            atom_names=True,
            plot_bonds=True,
            contours='80%',  # Show 80% of electron density
            auto_show=False
        )

        vis.mlab.savefig(ORBITAL_OUTPUT_FILE)
        vis.mlab.close()

        print(f"✅ Saved to {ORBITAL_OUTPUT_FILE}")

    else:
        print("\nOpening orbital visualization...")
        print("(Close the Mayavi window when done)\n")

        vis.plot_orbitals_MO(
            orbital_numbers=[homo_index],
            plot_atoms=True,
            atom_scaling=0.5,
            atom_names=True,
            plot_bonds=True,
            contours='80%',  # Show 80% of electron density
            auto_show=True
        )

        print("✅ Orbital visualization complete!")
        print("💡 The colored surfaces show electron density distribution")

# %%
