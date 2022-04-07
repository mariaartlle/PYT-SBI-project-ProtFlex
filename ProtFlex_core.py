''' packages '''
from ProtFlex_functions import *
from ProtFlex_blastp import *
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
                if line[0].startswith('>'): # filtrem que no sigiu multifasta i tngui format fasta
                    fasta_provided = file.read() # emmagatzemem seq fasta dins de fasta_provided
                else:
                    raise IncorrectInput(args.inputfile)
    else:
        raise IncorrectInput(args.inputfile)

''' 1 fasta a uniprotID'''
print('before blastp')
print(process_FASTA(fasta_provided))




''' 3 uniprot a alphafold '''
# os.mkdir("aplhafold")
#
#
# protein = parsePDB(protein.get_alphafold_structure())

''' 4 Retrieve Normalized Square fluctuactions from GNM model '''
#
# protein_NSF = get_NormSqFluct(protein)


''' Output formatting '''
#
# outputfile = open("")
