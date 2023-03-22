def letters( input):
    return ''.join(filter(str.isalpha, input))

def strip_input_name( input ):
    return

def read_mo_molpro( filename, text, nbasis ):
    """

    Python version on GAMMCOR function

    temporally in the utils

    """

    import numpy as np

    text_binary = text.encode('ascii') 

    cmo = np.zeros([nbasis,nbasis], dtype=np.float64)

    print("=====")
    with open(filename, "rb") as f:
        f.read(4)
        while (data := f.read(8)):
            
            if data ==  text_binary:
                # print( "I have found" )     
                
                f.read(8)    

                buffer = f.read(4)
                nsym = np.frombuffer(buffer, dtype=np.int32, count=1)[0]

                buffer = f.read(4*nsym)
                nbas = np.frombuffer(buffer, dtype=np.int32, count=nsym)

                buffer = f.read(4*nsym)
                offs = np.frombuffer(buffer, dtype=np.int32, count=nsym)

                ncmo_tmp = np.sum(nbas**2)

                f.read(8)  

                buffer = f.read(8*ncmo_tmp)
                cmo_tmp = np.frombuffer(buffer, dtype=np.float64, count=ncmo_tmp)

                # print(cmo_tmp)

                idx = 0 

                for i_rep in range(nsym):

                    i_off = offs[i_rep]
                    i_nbas = nbas[i_rep]

                    cmo[i_off:i_nbas, i_off:i_nbas] = cmo_tmp[idx:i_nbas*i_nbas].reshape([i_nbas,i_nbas])

                    idx += i_nbas * i_nbas
            break

    return cmo


def read_SAPTVIS( filename):



    from scipy.io import FortranFile
    import numpy as np


    f = FortranFile( filename, 'r')

    ANBasis, BNBasis = f.read_ints(np.int32) # type: ignore
    NOccupA, NOccupB = f.read_ints(np.int32) # type: ignore

    Occ_buff = f.read_reals(np.float64) # type: ignore

    A_Occ = Occ_buff[0:ANBasis]
    B_Occ = Occ_buff[ANBasis:ANBasis+BNBasis]

    ACMO = np.reshape(f.read_reals(np.float64), [ANBasis,ANBasis] ) # type: ignore
    BCMO = np.reshape(f.read_reals(np.float64), [BNBasis,BNBasis] ) # type: ignore
    Qmat = np.reshape(f.read_reals(np.float64), [NOccupA,NOccupB] )  # type: ignore

    return ANBasis, BNBasis, NOccupA, NOccupB, A_Occ, B_Occ, ACMO, BCMO, Qmat
