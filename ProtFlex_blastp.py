import Bio.Blast
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML


fasta_string = open("prova.fasta").read()
result_handle = NCBIWWW.qblast("blastp", "swissprot", fasta_string)
with open("my_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle = open("my_blast.xml")
blast_record = NCBIXML.read(result_handle)


#line code to obtain the best matches with an e-value threshold
for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        E_VALUE_THRESH = 0.0001
        if hsp.expect < E_VALUE_THRESH:
            print ("***Best match***")
            print("sequence:", alignment.title)
            print("length:", alignment.length)
            print("e value:", hsp.expect)
