# %%

import visualization.visualization as vis

# %%

visualization = vis.Visualization(input_type='Dalton', input_sub_type='tar', input_name='/SN7500_pool/Aga/dragon/cc2/ab2/start_mol_9.0')

# %%


visualization.get_geometry(get_bonds= False)

# %% 

visualization.molecular_system.atoms_Name.append('C2')
visualization.molecular_system.atoms_Name.append('C2')

visualization.molecular_system.atoms_Name.append('H2')
visualization.molecular_system.atoms_Name.append('H3')
visualization.molecular_system.atoms_Name.append('H2')
visualization.molecular_system.atoms_Name.append('H3')

# %%

visualization.plot_Geometry( Plot_Bonds=0 )



# %%

visualization.Plot_D_AB(atom_names= 1, atom_scaling= 0.3 , bond_scaling = 0.3 , contours = 10, Plot_bonds = False )

# %%

