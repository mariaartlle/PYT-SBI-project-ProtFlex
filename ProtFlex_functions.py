import urllib as urllib
import sys
# import ProDy as prody



class IncorrectInput(NameError):
	def __init__(self,input):
		self.input = input
	def __str__(self):
		return "The name %s is not a a formatted FASTA file " %(self.input)

class Protein(object):
    def __init__(self, uniprotID):
        self.uniprotID = uniprotID

    def get_alphafold_structure(self):
        url = ('https://alphafold.ebi.ac.uk/files/AF-%s-F1-model_v2.pdb' %(self.uniprotID))
        req = urllib.request.urlretrieve(url, 'alphafold/%s.pdb' %(self.uniprotID))
        return 'alphafold/%s.pdb' %(self.uniprotID)




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
