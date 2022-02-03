# %%

import visualization.visualization as vis

# %%

visualization = vis.Visualization(input_type='Dalton', input_sub_type='tar', input_name='/home/hemik/Hapka/dla_Piotra/ethylene')

# %%


visualization.get_geometry( )

# %%

visualization.get_orbital_data( )

# %%

visualization.orbital_generator.grid.R_max_multip = 5.0
visualization.orbital_generator.grid.x_n= 60 
visualization.orbital_generator.grid.y_n= 60
visualization.orbital_generator.grid.z_n= 60
visualization.orbital_generator.init_grid( )

visualization.orbital_generator.init_grid( )
visualization.orbital_generator.init_AOs()
visualization.generate_AO_orbitals()
visualization.generate_MO_orbitals()

# %% 

visualization.get_geminal_data()
visualization.get_generate_geminals()

# %%
visualization.visualization_data.Atoms_Color["C"] =  (0.5, 0.5, 0.5)
visualization.plot_Geminals( geminal_numbers= [0], plot_bonds=False, atom_names=False, auto_show=False )
visualization.mlab.view(azimuth=0.0, elevation=90.0, distance=18, focalpoint=[0.0, 0.0, 0.0] )
#visualization.mlab.savefig(filename='ethylene_90_1.pdf')
visualization.mlab.show()
# %% 
visualization.plot_Geminals( geminal_numbers= [1], plot_bonds=False, atom_names=False, auto_show=False )
visualization.mlab.view(azimuth=0.0, elevation=90.0, distance=18, focalpoint=[0.0, 0.0, 0.0] )
#visualization.mlab.savefig(filename='ethylene_90_1.pdf')
visualization.mlab.show()
# %%

visualization.plot_Orbitals( orbital_numbers=[5], plot_bonds=False )

