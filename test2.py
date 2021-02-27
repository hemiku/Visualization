# %%

import visualisation.visualisation as v

# %%

visualisation = v.Visualisation(input_type='Dalton', input_sub_type='tar', input_name='tests/GVB_AB_AB_restart_1.000')

# %%
visualisation.get_geometry()

# %%

#visualisation.plot_Geometry()

# %%

visualisation.get_orbital_data()
# %%


import visualisation.orbital_generator as og

 # %%
AO_parameter = og.AOParameters.sqrt15
print(AO_parameter)
# %%
OrbitalsGenrator = og.OrbitalsGenrator()
# %%
visualisation.molecular_system.get_Atoms_Name()
# %%
A = visualisation.molecular_system.get_Atoms_Name()

# %%
A
# %%
A[1]='ffff'
# %%
