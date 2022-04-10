### Modules ###
from ProtFlex_functions import *
# from GaussianNetworkModel import *
import argparse

### Input arguments ###
parser = argparse.ArgumentParser(description = "ProtFlex 1.0: python program to obtain flexibility scores for each aminoacid in protein sequences")

required_arg = parser.add_argument_group("Required arguments")

required_arg.add_argument('-i', '--input',
                    dest = "inputfile",
                    action = "store",
                    default = None,
                    required = True,
                    help = "Required argument. It must be a FASTA formatted file")

args = parser.parse_args()

### Logging file and stdout configuration ###

logging.basicConfig(level=logging.DEBUG,
                    filename="ProtFlex.log",
                    format='%(levelname)s: %(asctime)s %(message)s',
                    filemode='w')

# print info into stdout and logging file
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s: %(message)s')
console.setFormatter(formatter)
logging.getLogger().addHandler(console)

### Process input file ###
logging.info('Processing input file')
if args.inputfile:
    if os.path.isfile(args.inputfile):
        with open(args.inputfile) as file:
            for line in file:
                if line[0].startswith('>'): # filtrem que no sigiu multifasta i tngui format fasta
                    fasta_provided = file.read().strip() # emmagatzemem seq fasta dins de fasta_provided
                else:
                    logging.error('The file provided is not in FASTA format')
                    raise IncorrectInput(args.inputfile)

    else:
        logging.error('A file has not been provided')
        raise IncorrectInput(args.inputfile)


### Retrieve UniprotID ###

# protein = FASTA_to_UniprotID(fasta_provided)
protein = Protein('P06401', 'MTELKAKGPRAPHVAGGPPSPEVGSPLLCRPAAGPFPGSQTSDTLPEVSAIPISLDGLLFPRPCQGQDPSDEKTQDQQSLSDVEGAYSRAEATRGAGGSSSSPPEKDSGLLDSVLDTLLAPSGPGQSQPSPPACEVTSSWCLFGPELPEDPPAAPATQRVLSPLMSRSGCKVGDSSGTAAAHKVLPRGLSPARQLLLPASESPHWSGAPVKPSPQAAAVEVEEEDGSESEESAGPLLKGKPRALGGAAAGGGAAAVPPGAAAGGVALVPKEDSRFSAPRVALVEQDAPMAPGRSPLATTVMDFIHVPILPLNHALLAARTRQLLEDESYDGGAGAASAFAPPRSSPCASSTPVAVGDFPDCAYPPDAEPKDDAYPLYSDFQPPALKIKEEEEGAEASARSPRSYLVAGANPAAFPDFPLGPPPPLPPRATPSRPGEAAVTAAPASASVSSASSSGSTLECILYKAEGAPPQQGPFAPPPCKAPGASGCLLPRDGLPSTSASAAAAGAAPALYPALGLNGLPQLGYQAAVLKEGLPQVYPPYLNYLRPDSEASQSPQYSFESLPQKICLICGDEASGCHYGVLTCGSCKVFFKRAMEGQHNYLCAGRNDCIVDKIRRKNCPACRLRKCCQAGMVLGGRKFKKFNKVRVVRALDAVALPQPVGVPNESQALSQRFTFSPGQDIQLIPPLINLLMSIEPDVIYAGHDNTKPDTSSSLLTSLNQLGERQLLSVVKWSKSLPGFRNLHIDDQITLIQYSWMSLMVFGLGWRSYKHVSGQMLYFAPDLILNEQRMKESSFYSLCLTMWQIPQEFVKLQVSQEEFLCMKVLLLLNTIPLEGLRSQTQFEEMRSSYIRELIKAIGLRQKGVVSSSQRFYQLTKLLDNLHDLVKQLHLYCLNTFIQSRALSVEFPEMMSEVIAAQLPKILAGMVKPLLFHKK')

### Retrieve alphafold structure from UniprotID ###

protein_str = parsePDB(protein.get_alphafold_structure())


### Retrieve Normalized Square fluctuactions from protein using GNM model ###
logging.info('GNM model construction')

protein_NSF = get_NormSqFluct(protein_str, name = str(protein.get_uniprotID()))

### Output formatting ###

with open('out.txt', 'w') as outfile:
    outfile.write('ProtFlex parseable output file\n')
    outfile.write('AT\tAA\tCH\tN\tNormSqFluct\n')
    i = 0
    for element in get_calphaPDB(protein.get_alphafold_structure()):
        element_nice = element.split()[2:6]
        for x in element_nice:
            outfile.write('%s\t' %(x))
        outfile.write('%.5e\t' %(protein_NSF[i]))
        outfile.write('\n')
        i += 1

logging.info('ProtFlex successfully executed')
