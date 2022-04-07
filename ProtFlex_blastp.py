import Bio.Blast
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML







            # E_VALUE_THRESH = 0.0001
            # if hsp.expect < E_VALUE_THRESH:
            #     print ("***Best match***")
            #     print("sequence:", alignment.title)
            #     print("length:", alignment.length)
            #     print("e value:", hsp.expect)






#### Lo millor seria que la funció faci return d'un protein object: return Protein(UniprotID)

if __name__ == "__main__":
    process_FASTA('P06401.fa')
