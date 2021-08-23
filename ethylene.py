# %%

import visualization.visualization as v
import visualization.dispersion_plot as disp


# %%
#/SN7500_pool/Plot_dispersion/C2H4_C2H4/Piotr/frozen_4
#test_C2H4_C2H4_7_702.out', eerpa_path='wynik.out')
path = '/SN7500_pool/Plot_dispersion/C2H4_C2H4/Piotr/frozen_4/'
visualization = disp.DispersionPlot(input_type='Dalton', input_sub_type='tar', input_name= path + 'test_C2H4_C2H4_7_702')
visualization.set_GAMMCOR_filename(filename= path + 'wynik.out')
# %%

visualization.get_dispersion_index( x_n= 80, y_n= 80, z_n= 80, monomer_A=2, monomer_B=3 )


# %%

visualization.plot_Geometry(Plot_Bonds=False)

# %%



###3mlab.view(90, 70, 6.2, (-1.3, -2.9, 0.25))

#visualization.mlab.clf()

#figure.scene.camera.position = [-6.4, -12.0, 2.0]

#figure.scene.camera.view_angle = 30.0
#figure.scene.camera.view_up = [0.0, 0.0, 1.0]

# %%
figure = visualization.mlab.figure(   "Dispersion", 
                            bgcolor= visualization.visualization_data.background_colors['White'],
                            size=(600, 400) )


figure = visualization.Plot_D_AB(atom_names= 0, atom_scaling= 0.3 ,plot_bonds=False, bond_scaling = 0.3 , contours = 10, sclalarbar=True, auto_show=False, figure= figure)


view = visualization.mlab.view()
print( view )
# %%
#visualization.mlab.view(azimuth=90.0, elevation=-90.0, distance=11, focalpoint=[-3.8573824125626173, -0.011676371097564697, 0.3456050871419716] )
visualization.mlab.view(azimuth=90.0, elevation=45.0, distance=18, focalpoint=[-3.0, -0.0, -0.32] )
visualization.mlab.savefig(filename='ethylene_1.pdf')
visualization.mlab.show()

# %%


figure = visualization.mlab.figure(   "Dispersion", 
                            bgcolor= visualization.visualization_data.background_colors['White'],
                            size=(600, 400) )

figure = visualization.Plot_D_AB(atom_names= 0, atom_scaling= 0.3 ,plot_bonds=False, bond_scaling = 0.3 , contours = 10, sclalarbar=True, auto_show=False, figure= figure)

#visualization.mlab.view(azimuth=90.0, elevation=-90.0, distance=11, focalpoint=[-3.8573824125626173, -0.011676371097564697, 0.3456050871419716] )
visualization.mlab.view(azimuth=90.0, elevation=135.0, distance=18, focalpoint=[-3.0, -0.0, -0.32] )
visualization.mlab.savefig(filename='ethylene_2.pdf')
visualization.mlab.show()




# %%
view = visualization.mlab.view()

# %%

print( view )


visualization.plot_Geminals( geminal_numbers=[11], atom_names= 0, atom_scaling= 0.3 ,plot_bonds=False )


# %%



## TODO  ustawić ile geminalu ma być pod powierzchnią, 
# dobrać odpowiedznio legendę  
# Geminal jedna powierzchnia 
# dyspersja kontury
# sprawdzić mayavi latex



