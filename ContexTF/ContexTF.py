# IMPORTS
import os

from functions import *

import argparse, logging

# ARGUMENTS
def arguments():
    '''
    Define the input arguments of the program.
    '''

    parser = argparse.ArgumentParser(description = 'ContexTF 1.0: functional annotation of prokaryotic regulatory elements')

    required_arg = parser.add_argument_group('Required arguments')
    optional_arg = parser.add_argument_group('Optional arguments')

    required_arg.add_argument('-fna', '--fasta_file',
                    dest = 'fastafile',
                    action = 'store',
                    default = None,
                    required = True,
                    help = 'Required argument. It must be a genomic FASTA formatted file containing the whole genome to study.')

    required_arg.add_argument('-gff3', '--gff3_file',
                    dest = 'gff3file',
                    action = 'store',
                    default = None,
                    required = True,
                    help = 'Required argument. It must be a GFF3 formatted file containing the annotations of the genomic FASTA file.')
    
    optional_arg.add_argument('--tmalign',
                    dest = 'tmalign',
                    action = 'store_true',
                    default = None,
                    required = False,
                    help = 'Optional argument. If defined, ContexTF will perform structural comparisons between the reference TF and the candidates.')
    
    optional_arg.add_argument('--tmalign_min_threshold',
                    dest = 'tmalign_min_threshold',
                    type=float,
                    default=0.5,
                    required = False,
                    help = 'Optional argument. Takes a float between 0 and 1, o.5 as default. If defined, ContexTF will only consider as candidates superimpositions with TM-score above the threshold.')
    
    optional_arg.add_argument('--n_flanking',
                    dest = 'n_flanking',
                    type=int,
                    default=6,
                    required = False,
                    help = 'Optional argument. Number of flanking genes to plot upstream and downstream of each TF candidate, 4 as default value.')
    
    optional_arg.add_argument('--n_flanking3',
                    dest = 'n_flanking3',
                    type=int,
                    default=6,
                    required = False,
                    help = "Optional argument. Number of flanking genes to plot in the 3' site.")
    
    optional_arg.add_argument('--n_flanking5',
                    dest = 'n_flanking5',
                    type=int,
                    default=6,
                    required = False,
                    help = "Optional argument. Number of flanking genes to plot in the 5' site.")
    
    return parser.parse_args()

# LOGGING FILE
def logger():
    '''
    Logging file and stdout configuration.
    '''
    logging.basicConfig(level=logging.INFO,
                        filename='ContexTF.log',
                        format='%(levelname)s|%(asctime)s: %(message)s',
                        datefmt='%Y-%m-%d|%H:%M:%S',
                        filemode='w')

    # print info into stdout and logging file
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter("%(levelname)s|%(asctime)s: %(message)s",
                              "%Y-%m-%d|%H:%M:%S")
    console.setFormatter(formatter)
    logging.getLogger().addHandler(console)

# ERRORS
class IncorrectInput(NameError):
	def __init__(self,input):
		self.__input = input
	def __str__(self):
		return "%s is not an accepted input" %(self.__input)

def HMMER_module(organism_proteome, structural_superimposition=False, current_directory='.', executable_path='.', db_fasta='TF_databases/TF_database_merged.fasta', db_tsv='TF_databases/TF_database_merged.tsv'):
    '''
    When HMMER3.3.2 (http://hmmer.org/) is locally installed and provided with 2 fastas:
        1. Downloads and prepares the Pfam-A database.
        2. hmmscan the db_fasta for PFAM domains
            2.1 Parses the output file and only remain the non-overlapping domains with best e-value.
        3. hmmfetch the domains
        4. hmmsearch of the remaining domains against the input_fasta
            4.1 Parses the hmmsearch output file merging the information about the database tsv file and hmmscan domain information
            4.2 Produces an output file with only the hits which have all the domains found in the reference transcription factor
            4.3 If structural_superimposition == True, maps the target protein identifiers to uniprot
    '''

    if os.path.isdir('HMMER') == False:
        os.mkdir('HMMER')    
    os.chdir('HMMER')

    current_directory = current_directory+'/HMMER/'

    # 1. Prepare the Pfam database
    if os.path.exists(current_directory+'Pfam-A.hmm') == False:
        # Download PFAM database
        logging.info('Downloading Pfam-A database...')
        rec.urlretrieve('https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz', '{}Pfam-A.hmm.gz'.format(current_directory))
        os.system('gunzip {}Pfam-A.hmm.gz'.format(current_directory))
    else:
         logging.info('Pfam-A database found, using this local version to execute hmmscan.')

    # 1.1 Press the PFAM database
    for f in [current_directory+'Pfam-A.hmm.h3'+str(x) for x in ['i', 'p', 'm', 'f']]:
        if os.path.exists(f) == False:
            if glob.glob('*.h3*') != []:
                os.system('rm *.h3*')
            os.system('hmmpress {}Pfam-A.hmm'.format(current_directory))

    # 2. hmmscan of db_fasta
    if os.path.exists(current_directory+'/'+'hmmscan_output_parsed.tsv') == True:
        logging.info('hmmscan file found (version: {}), proceding with hmmsearch'.format(time.strftime('%m/%d/%Y|%H:%M:%S', time.gmtime(os.path.getmtime(current_directory+'/'+'hmmscan_output_parsed.tsv')))))
        hmmscan_file = current_directory+'/'+'hmmscan_output_parsed.tsv'
    else: 
        logging.info('Starting hmmscan...')
        os.system('hmmscan -E 0.001 --domtblout {}hmmscan_out.tsv {}Pfam-A.hmm {} > hmmscan.out'.format(current_directory, current_directory, executable_path+db_fasta))
        logging.info('hmmscan of {} finished.'.format(db_fasta.split('/')[-1]))

        # 2.1 Parse the output hmmscan file and remove overlapping domains with lower e-values
        hmmscan_file = non_overlapping_parser('hmmscan_out.tsv')
        logging.info('Parsing of hmmscan finished.')

    # 3. hmmfetch of the domains
    hmmfetch_file = retrieve_domains_from_hmmscan(hmmscan_file)
    if glob.glob('*.h3m.ssi') == False:
        os.system('hmmfetch --index Pfam-A.hmm')
    os.system('hmmfetch -o hmmfetch_out -f Pfam-A.hmm {}'.format(hmmfetch_file))

    # 4. Perform hmmsearch of the domains found in db_fasta into the target proteome
    if os.path.exists(current_directory+'/'+'hmmsearch_output_parsed.tsv') == True:
        logging.info('hmmsearch file found (version: {}), proceding with it'.format(time.strftime('%m/%d/%Y|%H:%M:%S', time.gmtime(os.path.getmtime(current_directory+'/'+'hmmscan_output_parsed.tsv')))))
        hmmsearch_file = current_directory+'/'+'hmmsearch_output_parsed.tsv'
    else:
        logging.info('Starting hmmsearch...')
        os.system('hmmsearch --domtblout hmmsearch_out.tsv -E 0.001 hmmfetch_out {} > hmmsearch.out'.format(organism_proteome))
        logging.info('hmmsearch of domains retrieved from {} into {} done'.format(db_fasta.split('/')[-1], organism_proteome.split('/')[-1]))
        hmmsearch_file = hmmsearch_parser('hmmsearch_out.tsv', hmmscan_file, executable_path+db_tsv)

    hmmer_df = only_hits_with_all_domains(hmmsearch_file)

    if structural_superimposition == True:
        # 1. Maps the target protein identifiers in the hmmsearch_file to UniprotIDs
        target_names = set(hmmer_df['Target_name'].to_list())
        mapped_dict = mapAccession2Uniprot(target_names)

        hmmer_df['Target_uniprot'] = ''
        for i in range(len(hmmer_df.axes[0])):
            accession = hmmer_df.at[i, 'Target_name']
            if accession in mapped_dict.keys():
                hmmer_df.at[i, 'Target_uniprot'] = mapped_dict[accession]
            else:
                hmmer_df.at[i, 'Target_uniprot'] = '-'

    hmmer_df.to_csv('hmmsearch_hits_all_domains_parsed.tsv', sep='\t', index=False)
    logging.info('Parsing of hmmsearch finished.')

    os.chdir('../')

    return hmmer_df

def TMAlign_module(hmmer_df, TMAlign_threshold=None):
    '''
    When TM-Align (https://zhanggroup.org/TM-align/) is installed on your local machine:
        1. Downloads the AlphaFold predicted structure of every protein
        2. Performs pair-wise alignments between the reference protein and every hit in the target proteome
        3. Produces a parseable output summary of the alignments with structural superimposition results
    '''
    if os.path.isdir('TMalign') == False:
        os.mkdir('TMalign')    
    os.chdir('TMalign')
    pwd = os.getcwd()
    # 1. Downloads the AlphaFold predicted structure of every protein
    superimposition_dict = get_alphafolds(hmmer_df, pwd=pwd)

    if 'best_alignments.tsv' in os.listdir(pwd):
        TMalign_results_file = pwd+'/best_alignments.tsv'
        logging.info('TM-Align results file found (version: {}), proceding with it'.format(time.strftime('%m/%d/%Y|%H:%M:%S', time.gmtime(os.path.getmtime(TMalign_results_file)))))
    else:
        # 2. Performs pair-wise alignments between the reference protein and every hit in the target proteome
        superimposition_dict = TM_align(superimposition_dict, TMAlign_folder=pwd)
        # 3. Produces a parseable output summary of the alignments with structural superimposition results
        TMalign_results_file = TM_align_parser(superimposition_dict, TMAlign_threshold=TMAlign_threshold)
        logging.info('Parsing of TMalign results done.')
    os.chdir('../')

    return TMalign_results_file

def GC_module(hmmer_df, path2gff=None, path2fna=None, TM_Align_file=None, n_flanking5=None, n_flanking3=None, executable_path='.'):
    '''
    GC_module can be executed with hmmer or TM-Align output (TM_Align_file)
        1. Downloads the genomic files for every protein 
        2. Retrieves the flanking genes upstream and downstream of every TF (n_flanking3 and n_flanking5)
    
    :return: a genomic context figure for each reference TF and all of its candidates
    '''

    if os.path.isdir('GC') == False:
        os.mkdir('GC')    
    os.chdir('GC')
    pwd = os.getcwd()
    # 1. Get the target dict with targets and gff and fna path for all targets
    targets_dict = retrieve_targets_dict(hmmer_df, path2gff=path2gff, path2fna=path2fna, TM_Align_file=TM_Align_file, executable_path=executable_path)

    # 2. Retrieve the genomic context for each target and candidates

    plasmids = ['AAP70493.1','ADO85569.1','AEM66515.1','AAA21920.1','BAD11074.1','AAP83141.1','BAB32408.1','AFI98560.1','ACO06750.1','AAB36583.1','AAK38101.1', 'ANB66399.1', 'AAG00065.1', 'AAA25445.1', 'AAA68939.2', 'CAC87048.1', 'BAA36282.1']
    
    if os.path.isdir('GC_figures') == False:
        os.mkdir('GC_figures')
    figures_dir = pwd+'/GC_figures/'

    already_figures = [str(x).strip('.png') for x in os.listdir(pwd+'/GC_figures/') if x.endswith('.png')]
 
    # print(len(plasmids))
    # print(len(targets_dict.keys()))
    # print(len(already_figures))

    for key in targets_dict.keys():
        if key not in plasmids and targets_dict[key]['TF'] not in already_figures:
            logging.info('GC module: plotting genome context for {}'.format(targets_dict[key]['TF']))
            all_syntenies = {}
            target_dict = {}
            target_dict[key] = targets_dict[key]

            all_syntenies = genomic_context(target_dict, n_flanking5=n_flanking5, n_flanking3=n_flanking3)

            # 3. Execute MMseqs2 o BLASTp /UniRep? to find families
            all_syntenies, protein_families_summary = find_and_add_protein_families(all_syntenies, out_label=targets_dict[key]['TF'], method='mmseqs')

            # 4. Draw the genomic context and save a figure
            draw_genomic_context(all_syntenies, protein_families_summary, targets_dict, save_dir=figures_dir)
            # Clean files 
            os.system('rm *.fa')
        elif key in plasmids:
            logging.info('No genomic context for {}: {} only sequence reference in NCBI databases is a plasmid/biosynthetic gene cluster ({}) it is not annotated in any species.'.format(targets_dict[key]['TF'], key, targets_dict[key]['fna'].split('/')[-1].strip('_genomic.fna ')))
    

def main():
    '''
    Execute the ContexTF pipeline. There's 3 different modules that can be executed independently
    when providing the needed files.
    '''

    # SETUP
    args = arguments()
    logger()
    logging.info('Time of start')

    # get path of the executable 
    executable_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))+'/'

    # establish current directory
    current_directory = os.getcwd()

    # Process arguments

    fastafile = current_directory + '/' + args.fastafile
    gff3file = current_directory + '/' + args.gff3file

    n_flanking = args.n_flanking
    n_flanking5 = args.n_flanking5
    n_flanking3 = args.n_flanking3
    if n_flanking3 == n_flanking5 and n_flanking != n_flanking3:
        if n_flanking != 4:
            n_flanking3 = n_flanking
            n_flanking5 = n_flanking


    # Execute pipeline
    if os.path.isfile(fastafile):
        if os.path.isfile(gff3file):

            # 1. Obtain proteome
            organism_proteome = obtain_cds_proteins_from_genome(gff3file, fastafile, current_directory=current_directory)

            if args.tmalign: 
                tmalign_min_threshold = args.tmalign_min_threshold

                # 2. HMMER
                hmmer_df = HMMER_module(organism_proteome, structural_superimposition=True, current_directory=current_directory, executable_path=executable_path)
                # 3. TM-Align
                TMalign_results_file = TMAlign_module(hmmer_df, TMAlign_threshold=tmalign_min_threshold)
                # 4. Genomic Context plotting of the TM_align candidates
                GC_module(hmmer_df, path2gff=gff3file, path2fna=fastafile, TM_Align_file=TMalign_results_file, n_flanking5=n_flanking5, n_flanking3=n_flanking3, executable_path=executable_path)
            else: 
                # 2. HMMER
                hmmer_df = HMMER_module(organism_proteome, current_directory=current_directory)
                # 3. Genomic Context plotting of the HMMER candidates
                GC_module(hmmer_df, path2gff=gff3file, path2fna=fastafile, n_flanking5=n_flanking5, n_flanking3=n_flanking3, executable_path=executable_path)
        else:
            logging.error('Please provide a GFF3 formatted file.')
            raise IncorrectInput(args.gff3file)
    else:
        logging.error('Please provide a genomic FASTA formatted file.')
        raise IncorrectInput(args.fastafile)

if __name__ == '__main__':
    main()
