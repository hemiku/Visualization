# %%

import visualization.visualization as v

# %%

visualization = v.Visualization(input_type='Dalton', input_sub_type='tar', input_name='tests/GVB_AB_AB_restart_1.000')

# %%
visualization.get_geometry()

# %%

#visualization.plot_Geometry()

# %%

visualization.get_orbital_data()
visualization.orbital_generator.init_grid( )

# %%

visualization.orbital_generator.init_AOs()
visualization.generate_AO_orbitals()
visualization.generate_MO_orbitals()

# %% 


visualization.plot_Orbitals(orbital_number=27)








# %%

visualization.orbital_generator.calc_AOs( AO= visualization.orbital_generator.AOs )


# %%

visualization.orbital_generator.AOs

# %% 

print(visualization.molecular_system.nb )



# %%

import visualization.orbital_generator as og

# %%
AO_parameter = og.AOParameters.sqrt15
print(AO_parameter)
# %%
OrbitalsGenrator = og.OrbitalsGenrator()
# %%
visualization.molecular_system.get_Atoms_Name()
# %%
A = visualization.molecular_system.get_Atoms_Name()

# %%
A
# %%
A[1]='ffff'
# %%


import visualization.inputs as inputs

# %%
inputs.INPUT_TYPES
# %%
data_input = inputs.get_input(input_type='Dalton', input_sub_type='tar', input_name='tests/GVB_AB_AB_restart_1.000')
# %%
print(data_input.__class__ )
# %%
data_input.__class__ = inputs.DaltonInput


# %%

visualization.data_input.__class__
# %%
data_input.get_nb()
# %%
