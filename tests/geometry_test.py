# %%

#import ...visualisation.visualisation as visualisation




# %%


visualisation = my_visualisation_library.visualisation.visualisation.Visualisation(input_type='Dalton', input_sub_type='tar', input_name='GVB_AB_AB_restart_1.000')


# %%

visualisation.data_input.get_nb()
visualisation.data_input.get_Coeff()

# %%

visualisation.data_input.get_nAtoms()
visualisation.data_input.get_Atoms()


# %%

visualisation.plot_Geometry()


# %%

print(visualisation.molecular_system)


# %%

import my_visualisation_library.visualisation.molecular_model as molecular_model

model = molecular_model.Molecular_System()




# %%

Atoms_R, Atoms_Charge,  Atoms_Name  = visualisation.data_input.get_Atoms()
model.set_Atoms( Atoms_R = Atoms_R, Atoms_Charge = Atoms_Charge,  Atoms_Name = Atoms_Name)




# %%

model.Atoms_R

# %%

visualisation.

# %%


Geminale = GG.GemialGenerator(Filename='GVB_AB_AB_restart_1.000.out')
Geminale.set_source(source='tar')
Geminale.get_data()


# %%
Geminale.Calc_AO_init2()
# %%
Geminale.Calc_MO_init()
# %%
Geminale.Calc_Geminals_init()
# %%
Geminale.Geminal_Generator()
# %%
Geminale.Plot_Geminal(Geminal_number=[0, 1, 2, 3, 4, 7], Plot_Bonds=1, Plot_Atoms=1, Plot_Atom_Names=0)
# %%
#import my_visualisation_library.visualisation.Grid
# %%
#grid = Grid.Grid( )
# %%

from my_visualisation_library.src.utils import letters

print(letters('4df322d'))
