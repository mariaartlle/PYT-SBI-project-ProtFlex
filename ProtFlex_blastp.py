import Bio.Blast
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML




def process_FASTA(fasta_provided):
    '''
    Performs a BLASTP of the sequence in the provided FASTA against the swissprot database.
    Returns the uniprotID of the hit with lowest e-value.
    '''
    print('lets blastp')
    result_handle = NCBIWWW.qblast("blastp", "swissprot", fasta_provided)
    print('blastp finished')
    with open("my_blast.xml", "w") as out_handle:
        out_handle.write(result_handle.read())
    result_handle = open("my_blast.xml")
    print('reading XML')
    blast_record = NCBIXML.read(result_handle)

    for alignment in blast_record.alignments:
        # for hsp in alignment.hsps:
        print('processing alignments')
        uniprotID = alignment.title[3:9]
        break
    return uniprotID


            # E_VALUE_THRESH = 0.0001
            # if hsp.expect < E_VALUE_THRESH:
            #     print ("***Best match***")
            #     print("sequence:", alignment.title)
            #     print("length:", alignment.length)
            #     print("e value:", hsp.expect)






#### Lo millor seria que la funció faci return d'un protein object: return Protein(UniprotID)

if __name__ == "__main__":
    process_FASTA('P06401.fa')
