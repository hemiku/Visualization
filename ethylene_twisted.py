# %%

import visualization.visualization as v
import visualization.dispersion_plot as disp


# %%

visualization = disp.DispersionPlot(input_type='Dalton', input_sub_type='tar', input_name='/SN7500_pool/Plot_dispersion/C2H4_C2H4/Twist/test')
visualization.set_GAMMCOR_filename(filename= '/SN7500_pool/Plot_dispersion/C2H4_C2H4/Twist/test_erpa.txt')
# %%

visualization.get_dispersion_index( x_n= 120, y_n= 120, z_n= 120, monomer_A=2, monomer_B=3 )


# %%

visualization.plot_Geometry(Plot_Bonds=False)

# %%

# %%
figure = visualization.mlab.figure(   "Dispersion", 
                            bgcolor= visualization.visualization_data.background_colors['White'],
                            size=(600, 400) )


figure = visualization.Plot_D_AB(atom_names= 0, atom_scaling= 0.3 ,plot_bonds=False, bond_scaling = 0.3 , contours = 10, sclalarbar=False, auto_show=False, figure= figure)
#visualization.mlab.view(azimuth=90.0, elevation=-90.0, distance=11, focalpoint=[-3.8573824125626173, -0.011676371097564697, 0.3456050871419716] )
#visualization.mlab.view(azimuth=90.0, elevation=45.0, distance=18, focalpoint=[-3.0, -0.0, -0.32] )
visualization.mlab.view(azimuth=90.0, elevation=135.0, distance=18, focalpoint=[-3.85, 0.0, 0.32] )
visualization.mlab.savefig(filename='ethylene_90_1.pdf')
visualization.mlab.show()

# %%


figure = visualization.mlab.figure(   "Dispersion", 
                            bgcolor= visualization.visualization_data.background_colors['White'],
                            size=(600, 400) )

figure = visualization.Plot_D_AB(atom_names= 0, atom_scaling= 0.3 ,plot_bonds=False, bond_scaling = 0.3 , contours = 10, sclalarbar=False, auto_show=False, figure= figure)

#visualization.mlab.view(azimuth=90.0, elevation=-90.0, distance=11, focalpoint=[-3.8573824125626173, -0.011676371097564697, 0.3456050871419716] )
#visualization.mlab.view(azimuth=90.0, elevation=135.0, distance=18, focalpoint=[-3.0, -0.0, -0.32] )
visualization.mlab.view(azimuth=90.0, elevation=135.0, distance=18, focalpoint=[-3.85, 0.0, 0.0] )
visualization.mlab.savefig(filename='ethylene_90_2.pdf')
visualization.mlab.show()




# %%
view = visualization.mlab.view()

# %%

#visualization.plot_Geminals( geminal_numbers=[11], atom_names= 0, atom_scaling= 0.3 ,plot_bonds=False, contours=[ 0.0014565149649289496] )

# %%


figure = visualization.mlab.figure(  "Geminal", 
                            bgcolor= visualization.visualization_data.background_colors['White'],
                            size=(600, 400) )

visualization.plot_Geminals( geminal_numbers=[11], atom_names= 0, contours=[0.002029361172288052],  atom_scaling= 0.3 ,plot_bonds=False, auto_show=False, figure=figure )
visualization.mlab.view(azimuth=90.0, elevation=135.0, distance=18, focalpoint=[-3.85, -0.0, -0.0] )
visualization.mlab.savefig(filename='ethylene_twisted_geminal.pdf')
visualization.mlab.show()

# %%


figure = visualization.mlab.figure(  "Geminal", 
                            bgcolor= visualization.visualization_data.background_colors['White'],
                            size=(600, 400) )

visualization.plot_Geminals( geminal_numbers=[11], atom_names= 0, contours=[0.002029361172288052],  atom_scaling= 0.3 ,plot_bonds=False, auto_show=False, figure=figure )
visualization.mlab.view(azimuth=90.0, elevation=45.0, distance=18, focalpoint=[-3.85, -0.0, -0.0] )
visualization.mlab.savefig(filename='ethylene_twisted_geminal_2.pdf')
visualization.mlab.show()


# %%



## TODO  ustawić ile geminalu ma być pod powierzchnią, 
# dobrać odpowiedznio legendę  
# Geminal jedna powierzchnia 
# dyspersja kontury
# sprawdzić mayavi latex



