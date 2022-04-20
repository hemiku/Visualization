# %% 

from visualization.visualization import Visualization

# %%
h2o_h2o = Visualization(input_type='Dalton', input_sub_type='tar', input_name= 'Examples/H2O_H2O')

# %%
h2o_h2o.get_geometry( )

# %%

h2o_h2o.plot_Geometry( )

# %%

h2o_h2o.plot_Geometry( background_color=(0.0, 1.0 ,0.0))

# %%

h2o_h2o.plot_Geometry( background_color=(0.0, 0.0 , 1.0))

# %%

h2o_h2o.plot_Geometry(atom_scaling = 3)
# %%

h2o_h2o.plot_Geometry()

# %%

h2o_h2o.get_orbitals( )
# %%
h2o_h2o.plot_orbitals_MO(orbital_numbers=[6], atom_names=False, plot_atoms=True, atom_scaling=0.2, plot_bonds=False)

# %%
# %%
h2o_h2o.plot_orbitals_MO(orbital_numbers=[8], atom_names=False, plot_atoms=True, atom_scaling=0.2, plot_bonds=False)

# %%
h2o_h2o.plot_grid( )
# %%
h2o_h2o.get_orbitals(R_max_multip=1.0 )

# %%
h2o_h2o.plot_orbitals_MO(orbital_numbers=[8], atom_names=False, plot_atoms=True, atom_scaling=0.2, plot_bonds=False)

# %%
