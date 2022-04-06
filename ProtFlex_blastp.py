import Bio.Blast
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML




def process_FASTA(fasta_provided):
    '''
    Performs a BLASTP of the provided FASTA against the swissprot database.
    Returns the uniprotID of the hit with lowest e-value.
    '''
    fasta_string = open(fasta_provided).read()
    result_handle = NCBIWWW.qblast("blastp", "swissprot", fasta_string)
    result_handle = open("my_blast.xml")
    blast_record = NCBIXML.read(result_handle)

    for hsp in alignment.hsps:
        E_VALUE_THRESH = 0.0001
        if hsp.expect < E_VALUE_THRESH:
            print ("***Best match***")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)


# # line code to obtain the best matches with an e-value threshold
#     for hsp in alignment.hsps:
#

if __name__ == "__main__":
    process_FASTA('P06401.fa')
