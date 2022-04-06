''' packages '''
from ProtFlex_functions import *
import sys
import os
import argparse

''' Input arguments '''

parser = argparse.ArgumentParser(description = "ProtFlex 1.0: python program to obtain flexibility scores for each aminoacid in protein sequences")

required_arg = parser.add_argument_group("Required arguments")

required_arg.add_argument('-i', '--input',
                    dest = "inputfile",
                    action = "store",
                    default = None,
                    required = True,
                    help = "Required argument. It must be a FASTA formatted file")

args = parser.parse_args()

if args.inputfile:
    if os.path.isfile(args.inputfile):
        with open(args.inputfile) as file:
            for line in file:
                if line[0].startswith('>') & (file.read()).count('>') == 0: # filtrem que no sigiu multifasta i tngui format fasta
                    fasta_provided = file.read() # emmagatzemem seq fasta dins de fasta_provided
                else:
                    raise IncorrectInput(args.inputfile)
    else:
        raise IncorrectInput(args.inputfile)

''' 1 fasta a uniprotID'''





''' 3 uniprot a alphafold '''


''' 4 alphafold aconseguir B-factors '''
