# %%
import Kowalski_Geminale as KG
# %%
Geminale = KG.GemialGenerator(Filename='GVB_AB_AB_restart_1.000.out')
#%% 

Geminale.set_source(source= 'tar')

# %%
Geminale.get_data()
# %%
Geminale.Calc_AO_init2()
Geminale.Calc_MO_init()

# %%
Geminale.Calc_Geminals_init()
Geminale.Geminal_Generator()
# %%
Geminale.Plot_Geminal(Geminal_number=[7], Plot_Bonds=1, Plot_Atoms= 1,Plot_Atom_Names=0)
# %%
