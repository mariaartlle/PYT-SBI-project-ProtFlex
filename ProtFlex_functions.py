import urllib as urllib
import sys
import ProDy as prody

i = 0

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


# Faig una llista i hi guardo tots els uniprotIDs del document input com protein objects
input_query = []
for protein in get_IDs_from_file("uniprotIDs.txt"):
    input_query.append(protein.uniprotID)



if __name__ == "__main__":
    print(str(protein.get_alphafold_structure()))
