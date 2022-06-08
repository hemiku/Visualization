def letters( input):
    return ''.join(filter(str.isalpha, input))

def strip_input_name( input ):
    return

def read_mo_molpro( filename, text, nbasis ):
    """

    Python version on GAMMCOR fuction

    temporaly in the utils

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