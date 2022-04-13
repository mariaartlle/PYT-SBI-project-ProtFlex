### Modules ###
from ProtFlex_functions import *
import argparse

def arguments():
    '''
    Define the input arguments of the program.
    '''
    parser = argparse.ArgumentParser(description = 'ProtFlex 1.0: python program to obtain flexibility scores for each aminoacid in protein sequences.')

    required_arg = parser.add_argument_group('Required arguments')
    optional_arg = parser.add_argument_group('Optional arguments')

    required_arg.add_argument('-i', '--input',
                        dest = 'inputfile',
                        action = 'store',
                        default = None,
                        required = True,
                        help = 'Required argument. It must be a FASTA formatted file.')

    optional_arg.add_argument('-o', '--output',
                        dest = 'outputfile',
                        action = 'store',
                        default = None,
                        required = False,
                        help = 'Optional argument. If defined, the program will overwrite the existing file or create a new one with the name provided.')

    optional_arg.add_argument('-g', '--graph',
                        dest = 'graph',
                        action = 'store_true',
                        default = None,
                        required = False,
                        help = 'Optional argument. If defined, the program will provide a graphical representation of the flexibility scores.')

    return parser.parse_args()

def logger():
    '''
    Logging file and stdout configuration.
    '''
    logging.basicConfig(level=logging.INFO,
                        filename='ProtFlex.log',
                        format='%(levelname)s: %(asctime)s %(message)s',
                        filemode='w')

    # print info into stdout and logging file
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(levelname)s: %(message)s')
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

    return None

def main():

    args = arguments()
    logger()

    ### Process input file ###
    logging.info('Processing input')
    if args.inputfile:
        if os.path.isfile(args.inputfile):
            with open(args.inputfile) as file:
                for line in file:
                    if line[0].startswith('>'):
                        fasta_provided = file.read().strip()
                        logging.info('%s FASTA file provided' %(args.inputfile))
                        protein = FASTA_to_UniprotID(fasta_provided)
                    else:
                        logging.error('The file provided is not in FASTA format')
                        raise IncorrectInput(args.inputfile)
        elif len(args.inputfile) == 6:
            protein = Protein(args.inputfile)
            logging.info('%s UniprotID provided' %(args.inputfile))

        else:
            logging.error('A file or UniprotID has not been provided')
            raise IncorrectInput(args.inputfile)


    ### Retrieve alphafold structure from UniprotID ###
    protein_str = protein.get_alphafold_structure()
    protein_parsed = parsePDB(protein_str)

    ### Retrieve Normalized Square fluctuactions from protein using GNM model ###

    logging.info('GNM model construction')
    if args.graph:
        protein_NSF = get_NormSqFluct(protein_parsed, graph = 1, name = str(protein.get_uniprotID()))
    else:
        protein_NSF = get_NormSqFluct(protein_parsed, name = str(protein.get_uniprotID()))

    ### Output formatting ###

    if args.outputfile:
        fd = args.outputfile
    else:
        fd = 'ProtFlex_%s_out.txt' %(protein.get_uniprotID())

    with open(fd, 'w') as outfile:
        outfile.write('ProtFlex parseable output file\n')
        outfile.write('Analyzed protein: %s\n' %(protein.get_uniprotID()))
        outfile.write('AT\tAA\tCH\tN\tpLDDT\tNormSqFluct\n')
        i = 0
        for element in get_calphaPDB(protein_str):
            element_nice = element.split()[2:6]
            for x in element_nice:
                outfile.write('%s\t' %(x))
            outfile.write('%s\t' %(element.split()[10]))
            outfile.write('%.5e\t' %(protein_NSF[i]))
            outfile.write('\n')
            i += 1
    os.remove('%s' %(protein_str))
    logging.info('Output saved in %s' %(fd))
    logging.info('ProtFlex successfully executed')

if __name__ == '__main__':
    main()
