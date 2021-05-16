# %%

import visualization.dispersion_plot as disp

# %%

visualization = disp.DispersionPlot(input_type='Dalton', input_sub_type='tar', input_name='tests/GVB_AB_AB_restart_1.000')
visualization.set_GAMMCOR_filename(filename= 'tests/wynik_gammcor_10.txt')
# %%

visualization.get_dispersion_index(  )


# %%

visualization.Plot_D_AB(atom_names= 1, atom_scaling= 0.3 , bond_scaling = 0.3 , contours = 10 )

# %%

