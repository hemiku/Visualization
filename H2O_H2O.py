# %%

import visualization.visualization as v
import visualization.dispersion_plot as disp


# %%

visualization = disp.DispersionPlot(input_type='Dalton', input_sub_type='tar', input_name='/SN7500_pool/Plot_dispersion/H2O_H2O/test')
visualization.set_GAMMCOR_filename(filename= '/SN7500_pool/Plot_dispersion/H2O_H2O/eerpa.out')
# %%

visualization.get_dispersion_index( x_n= 40, y_n= 40, z_n= 40, monomer_A=2, monomer_B=3 )


# %%

figure = visualization.mlab.figure(   "Dispersion", 
                            bgcolor= visualization.visualization_data.background_colors['White'],
                            size=(600, 400) )

visualization.mlab.clf()

figure.scene.camera.position = [-4.4, -12.0, 2.0]
#figure.scene.camera.focal_point = [-3.5, 0.0, 0.0] 
figure.scene.camera.view_angle = 30.0

figure.scene.camera.view_up = [0.0, 0.0, 1.0]

# %%
figure = visualization.Plot_D_AB(atom_names= 0, atom_scaling= 0.1 , bond_scaling = 0.3 , contours = 10, sclalarbar=True, auto_show=False, figure= figure)


# %%


visualization.mlab.show()

# %%


# %%


#  label_text_property.italic = False
#.label_text_property.bold = False 
# module_manager14.scalar_lut_manager.number_of_labels = 14
# module_manager14.scalar_lut_manager.data_range = array([1.47021721e-13, 6.00000000e-04])
# module_manager14.scalar_lut_manager.scalar_bar.unconstrained_font_size = False





# %%
plot_atoms = 1
atom_names = 1
plot_bonds = 1 
atom_scaling = 1.0
bond_scaling = 1.0
contours = 6
background_color = None
sclalarbar = True


_background_color = visualization.visualization_data.background_colors['White']


visualization.mlab.figure(   "Dispersion", 
                            bgcolor=_background_color,
                            size=(1000, 1000) )


if plot_atoms:
    for i in range(visualization.molecular_system.nAtoms):
                visualization.mlab.points3d(    visualization.molecular_system.atoms_R[i, 0],
                                                visualization.molecular_system.atoms_R[i, 1],
                                                visualization.molecular_system.atoms_R[i, 2],
                                                scale_factor=visualization.visualization_data.Atoms_Scale[ visualization.u.letters(visualization.molecular_system.atoms_Name[i])],
                                                resolution=20,
                                                color=visualization.visualization_data.Atoms_Color[ visualization.u.letters(visualization.molecular_system.atoms_Name[i])],
                                                scale_mode='none')

X, Y, Z = visualization.orbital_generator.grid.return_grid_arrays()

D_AB_contur = visualization.mlab.contour3d(X, Y, Z, (visualization.D_AB) ,contours=contours,opacity=0.5) 
    
if sclalarbar:
    _scalarbar = visualization.mlab.scalarbar(object=D_AB_contur, title='D^{AB}' ,orientation='vertical' )
    _scalarbar.scalar_bar.label_text_property.color = (0.0, 0.0, 0.0) 
    _scalarbar.scalar_bar.label_text_property.font_size  = 24
    #_scalarbar.scalar_bar.annotation_text_property.color = ( 0.0, 0.0, 0.0)



visualization.mlab.show()



# %%


dir(D_AB_contur.scene.camera.position.view)

# %%


D_AB_contur.scene.camera.position.view()

# %%

dir(D_AB_contur.scene.camera.pitch)

# %%

D_AB_contur.scene.camera.pitch()

# %%


_scalarbar.scalar_bar.annotation_text_property.color
# %%

dir(_scalarbar.scalar_bar.label_text_property)

# %%
_scalarbar.scalar_bar.label_text_property.font_size 

# %% 
)
_scalarbar.label_text_property = (0.0, 0.0,0.0) 


scalar_bar.number_of_labels


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



