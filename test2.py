# %%

import visualization.visualization as v
import visualization.dispersion_plot as disp

# %%

visualization = v.Visualization(input_type='Dalton', input_sub_type='tar', input_name='tests/GVB_AB_AB_restart_1.000')
visualization = disp.DispersionPlot(input_type='Dalton', input_sub_type='tar', input_name='tests/GVB_AB_AB_restart_1.000')
visualization.set_GAMMCOR_filename(filename= 'tests/wynik_gammcor_10.txt')
# %%

visualization.get_dispersion_index(  )


# %%

visualization.Plot_D_AB(atom_names= 1, atom_scaling= 0.3 , bond_scaling = 0.3 , contours = 10, sclalarbar=True)

# %%

visualization.orbital_generator.grid.x_n


# %%



visualization.get_geometry()

# %%

visualization.get_orbital_data()
visualization.orbital_generator.init_grid( )

# %%

visualization.orbital_generator.init_AOs()
visualization.molecular_system.grid = visualization.orbital_generator.grid
visualization.generate_AO_orbitals()
visualization.generate_MO_orbitals()

# %% 
#visualization.plot_Orbitals(orbital_numbers=[1])

# %% 
visualization.get_geminal_data()
visualization.get_generate_geminals()

# %% 

#visualization.plot_Geminals( geminal_numbers= [0,10,11] )

# %%

#visualization.plot_Geminals( geminal_numbers= visualization.np.arange(6) )

# %%
# %%
visualization.get_NumberOfFragments()
# %%
Out = visualization.get_GAMMCOR_Output()
# %%
visualization.get_NumberOfFragments()
visualization.get_Fragments()
visualization.get_Epsilons()
# %%
visualization.Monomers
# %%
Monomers_set = set(visualization.Monomers )
visualization.Calc_D( Monomer_A= Monomers_set.pop(), Monomer_B=Monomers_set.pop() )
# %%
visualization.Plot_D_AB(atom_names= 1, atom_scaling= 0.3 , bond_scaling = 0.3 , contours = 10, sclalarbar=True )
# %%


## TODO  ustawić ile geminalu ma być pod powierzchnią, 
# dobrać odpowiedznio legendę  
# Geminal jedna powierzchnia 
# dyspersja kontury
# sprawdzić mayavi latex



