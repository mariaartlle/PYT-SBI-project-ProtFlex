import sys
import os
import urllib as urllib
import logging
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from prody import *
import matplotlib.pylab as pl

class IncorrectInput(NameError):
	def __init__(self,input):
		self.__input = input
	def __str__(self):
		return "%s is not a formatted FASTA file" %(self.__input)

class Protein(object):
    def __init__(self, uniprotID, sequence):
        self.__uniprotID = uniprotID
        self.__sequence = sequence

    def __str__(self):
        return '(%s, %s)' %(self.__uniprotID, self.__sequence)

    def get_uniprotID(self):
        return self.__uniprotID

    def get_sequence(self):
        return self.__sequence

    def get_alphafold_structure(self):
        '''
        Retrieve structure with UniprotID from alphafold database
        '''
        url = ('https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v2.pdb' %(self.__uniprotID))
        req = urllib.request.urlretrieve(url, '%s.pdb' %(self.__uniprotID))
        logging.info('Alphafold structure retrieved')
        return '%s.pdb' %(self.__uniprotID)

def FASTA_to_UniprotID(fasta_provided):
    '''
    Performs a BLASTP of the sequence in the provided FASTA against the swissprot database.
    Returns a Protein object with the uniprotID of the hit with lowest e-value and its sequence.
    '''
    logging.info('Starting BLASTP')
    result_handle = NCBIWWW.qblast("blastp", "swissprot", fasta_provided)
    with open("my_blast.xml", "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle = open("my_blast.xml")
    blast_record = NCBIXML.read(result_handle)

    for alignment in blast_record.alignments:
        uniprotID = alignment.title[3:9]
        break
    logging.info('BLASTP finished')
    os.remove("my_blast.xml")
    return Protein(uniprotID, sequence = fasta_provided)


def get_NormSqFluct(protein):
    calphas = protein.select('calpha')
    gnm = GNM()
    gnm.buildKirchhoff(calphas)
    gnm.calcModes(None)
    SqFlucts = calcSqFlucts(gnm[0])
    NormSqFlucts = SqFlucts / (SqFlucts**2 ).sum()**0.5
    return NormSqFlucts


def get_IDs_from_file(uniprot_filename):
    ''' A generator function that reads a file with uniprotIds.
    In each iteration the function returns a protein object'''
    with open(uniprot_filename, 'r') as file:
        for line in file:
            uniprotID = line.strip()
            if len(uniprotID) > 0:
                try:
                    yield Protein(uniprotID)
                except:
                     sys.stderr.write("Incorrect Uniprot ID file input\n")




if __name__ == "__main__":
    print(str(protein.get_alphafold_structure()))
