import sys
import os
import urllib as urllib
import logging
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from prody import *
import matplotlib.pylab as pl
import re

class IncorrectInput(NameError):
	def __init__(self,input):
		self.__input = input
	def __str__(self):
		return "%s is not a formatted FASTA file" %(self.__input)

class Protein(object):
    def __init__(self, uniprotID):
        self.__uniprotID = uniprotID

    def __str__(self):
        return '(%s)' %(self.__uniprotID)

    def get_uniprotID(self):
        return self.__uniprotID

    def get_alphafold_structure(self):
        '''
        Retrieve structure with UniprotID from alphafold database.
        '''
        url = ('https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v2.pdb' %(self.__uniprotID))
        req = urllib.request.urlretrieve(url, '%s.pdb' %(self.__uniprotID))
        logging.info('%s Alphafold structure retrieved' %(self.__uniprotID))
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
    return Protein(uniprotID)


def get_NormSqFluct(protein, graph = None, name = 'Unknown'):
	'''
	Performs a Gaussian Network Model from PDB coordinates and returns the
	normalized Square fluctuations for each alpha carbon of the protein provided.
	The SqFlucts are calculated using the first mode (the slowest one).
	Also saves in the current directory a graphical representation if defined.
	'''
	calphas = protein.select('calpha')
	gnm = GNM(name = name)
	gnm.buildKirchhoff(calphas)
	gnm.calcModes(None)
	hinges = calcHinges(gnm[0])
	SqFlucts = calcSqFlucts(gnm[0])
	NormSqFlucts = SqFlucts / (SqFlucts**2 ).sum()**0.5

	if graph != None:
		pl.figure()
		showNormedSqFlucts(gnm[0], hinges=True)
		pl.savefig('%s_NormSqFlucts.png' %(name))
		logging.info('Graphical representation of protein flexibility in %s_NormSqFlucts.png' %(name))

	return NormSqFlucts


def get_calphaPDB(pdbfile):
	'''
	Generator function to obtain a alpha carbon pdb file.
	When provided a pdb filepath, yields each alpha carbon line.
	'''
	with open(pdbfile, 'r') as fd:
		for line in fd:
			if re.search(r'^ATOM\s+\d+\s+CA\s+', line):
				yield line
