# %%

import Kowalski_Geminale as KG


# %%
Geminale = KG.GemialGenerator(Filename='DALTON')


Geminale.set_source(source='MOPUN')

# %%

Geminale.get_data()

# %%
Geminale.set_Init_Grid_param(R_max_multip=2.0, x_n=50, y_n=50, z_n=50)
Geminale.Prepe_Grid()
Geminale.Gen_Grid()

# %%
Geminale.Calc_AO_init2()
Geminale.Calc_MO_init()
Geminale.Calc_Geminals_init()
# %%
%%time
Geminale.Geminal_Generator(x_n=50, y_n=50, z_n=50)
# %%
Geminale.Plot_Geminal(Geminal_number=[6], Plot_Bonds=0)



# %%
