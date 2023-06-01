# %%
import visualization.input_molpro as input_molpro
import visualization.input_molpro_sapt as input_sapt
import visualization.visualization as V
import importlib

import visualization.utils as utils

import numpy as np

# %% 


path_g = '/run/user/1000/gvfs/smb-share:server=truenas.local,share=aga/Visualization/SAPT_ELECT/test/'
molpro_output_file = 'water'

sapt_A_g = V.Visualization( input_type='MolproSapt', input_sub_type='output',  input_name= path_g + molpro_output_file   )

sapt_B_g = V.Visualization( input_type='MolproSapt', input_sub_type='output',  input_name= path_g + molpro_output_file   )
sapt_B_g.data_input.monomer = 1 

# %%
sapt_A_g.get_geometry()
sapt_B_g.get_geometry()
# %%
sapt_A_g.molecular_system.atoms_Charge
# %%
sapt_B_g.molecular_system.atoms_Charge
# %%
sapt_A_g.get_orbital_data( )

# %%
sapt_A_g.orbital_generator.grid.R_max_multip = 2
sapt_A_g.orbital_generator.grid.x_n = 30
sapt_A_g.orbital_generator.grid.y_n = 30
sapt_A_g.orbital_generator.grid.z_n = 30
sapt_A_g.orbital_generator.init_grid( )
sapt_A_g.orbital_generator.init_AOs()
sapt_A_g.orbital_generator.spherical  = 1 
# %%
# sapt_A.molecular_system.Coeff[:,:] = fortran.types.read_mo_molpro(  path+'MOLPRO_A.MOPUN', 'CASORB  ', sapt_A.molecular_system.nb).transpose()

# %%
sapt_B_g.get_orbital_data( )

sapt_B_g.orbital_generator.grid.R_max_multip = 2
sapt_B_g.orbital_generator.grid.x_n = 30
sapt_B_g.orbital_generator.grid.y_n = 30
sapt_B_g.orbital_generator.grid.z_n = 30
sapt_B_g.orbital_generator.init_grid( )
sapt_B_g.orbital_generator.init_AOs()
sapt_B_g.orbital_generator.spherical  = True 
# 
# sapt_B.molecular_system.Coeff[:,:] = fortran.types.read_mo_molpro(  path+'MOLPRO_B.MOPUN', 'CASORB  ', sapt_B.molecular_system.nb).transpose()

# %%

ANBasis, BNBasis, NOccupA, NOccupB, A_Occ, B_Occ, ACMO, BCMO, Qmat, QelA, QelB = utils.read_SAPTVIS_ELECTRO( path_g + 'SAPT_ELECT' )


sapt_A_g.molecular_system.Coeff[:,:] = ACMO
sapt_B_g.molecular_system.Coeff[:,:] = BCMO

# %%
# %%

sapt_A_g.generate_AO_orbitals()
sapt_A_g.generate_MO_orbitals()
# %%

sapt_B_g.generate_AO_orbitals()
sapt_B_g.generate_MO_orbitals()



# %%
sapt_electrostatic_A =  np.zeros_like( sapt_A_g.molecular_system.AOs[0] )
sapt_electrostatic_B  =  np.zeros_like( sapt_B_g.molecular_system.AOs[0] )

for i in range(NOccupA):
    sapt_electrostatic_A += QelA[i] * sapt_A_g.molecular_system.MOs[i]**2 

for i in range(NOccupB):
    sapt_electrostatic_B += QelB[i] * sapt_B_g.molecular_system.MOs[i]**2 



# %%

sapt_electrostatic_AB = - ( sapt_electrostatic_A + sapt_electrostatic_B )

# %%


# %%
import numpy as np 

X, Y, Z = sapt_A_g.orbital_generator.grid.return_grid_arrays()
print("dispersion_AB", 1000*np.sum(sapt_electrostatic_AB)*(X[1,0,0]-X[0,0,0]) *( Y[0,1,0]-Y[0,0,0]) * (Z[0,0,1]-Z[0,0,0]))
# %%

dx = (sapt_A_g.orbital_generator.grid.x_max - sapt_A_g.orbital_generator.grid.x_min)/ (sapt_A_g.orbital_generator.grid.x_n +1)
dy = (sapt_A_g.orbital_generator.grid.y_max - sapt_A_g.orbital_generator.grid.y_min)/ (sapt_A_g.orbital_generator.grid.y_n +1 )
dz = (sapt_A_g.orbital_generator.grid.z_max - sapt_A_g.orbital_generator.grid.z_min)/ (sapt_A_g.orbital_generator.grid.z_n +1 )

dv = dx*dy*dz

# %%
X, Y, Z = sapt_A_g.orbital_generator.grid.return_grid_arrays()
plot_figure = sapt_A_g.plot_orbitals_MO(orbital_numbers=[],auto_show=False, plot_bonds=False, atom_names=False, atom_scaling=0.5)
# plot_figure = sapt_A_g.plot_Orbitals_AO(orbital_numbers=[1],auto_show=False, plot_bonds=False, atom_names=False, atom_scaling=0.5, sclalarbar=0)
sapt_A_g.mlab.contour3d( X, Y, Z, sapt_electrostatic_AB,  contours= sapt_A_g.contour_process( ['50%','40%','30%','20%','10%','1%'], sapt_electrostatic_AB) , opacity=0.5)

sapt_A_g.mlab.show()

# %% 

