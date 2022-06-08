# %%
from numpy import f2py
# %%
gammcor_source_path = '/home/hemik/GammCor/gammcor/SOURCE/'

# %%
with open(gammcor_source_path + "types_f2py.f90") as sourcefile:
    sourcecode = sourcefile.read()
# %%
with open(gammcor_source_path + "types.f90") as sourcefile:
    sourcecode = sourcefile.read()

# %%
f2py.compile(sourcecode, modulename='types', extension='.f90')
# %%



# %%

