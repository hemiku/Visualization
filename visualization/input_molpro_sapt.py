
from visualization.input_molpro import MolproInput

class MolproSaptInput(MolproInput):
    """

        Input from molpro

    """
    input_name = 'Molpro'
    monomer:int = 0


    def __init__(self, *args, **kwargs):
        super(MolproSaptInput, self).__init__(*args, **kwargs)

        if 'monomer' in kwargs:
            self.monomer = kwargs['monomer']


    def get_output(self):

        if self.output is not None:

            return self.output

        else:
            with open(self.input_name + ".out", 'r', encoding="utf-8") as f:
                _output = f.read()

            geometry_block_beginning_sentence = '1PROGRAM * SEWARD (Integral evaluation for generally contracted gaussian basis sets)     Author: Roland Lindh, 1990' 

            self.output = _output.split(geometry_block_beginning_sentence)[ 1+ self.monomer ]

            return self.output