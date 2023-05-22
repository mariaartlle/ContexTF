# IMPORTS
import os, argparse, logging, re, time, json, zlib, operator, logging, glob

import urllib.request as rec
import pandas as pd
import numpy as np
import subprocess as sp
import soothsayer_utils as so
import matplotlib as plt
import matplotlib.cm as cm

from functools import reduce
from operator import itemgetter
from itertools import groupby, chain
from collections import Counter
from Bio import SeqIO, Seq
from Bio.SeqRecord import SeqRecord
from BCBio import GFF
from Bio.Blast import NCBIXML
from scipy.cluster import hierarchy
from scipy.spatial import distance
from pygenomeviz import GenomeViz
from matplotlib.lines import Line2D
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry

pd.options.mode.chained_assignment = None  # default='warn'

######################################################## PIPELINE FUNCTIONS ############################################################

### 0. UTILITIES ###
def get_subdirs_list(pwd='.'):
    folders_list = []
    os.chdir(pwd)
    for TMAlign_folder, dirs, files in os.walk(pwd):
        for subdir in dirs:
            folders_list.append(subdir)
    return folders_list

def divide_chunks(l, n):
    # looping till length l
    for i in range(0, len(l), n):
        yield l[i:i + n]

def breakby(biglist, sep, delim=None):
    for item in biglist:
        p = item.split(sep)
        if p[0] != '':
            yield p[0]
        if len(p) > 1:
            yield delim
            if p[1] != '':
                yield p[1]

def divide_list_consecutive_numbers(data_list):
    '''
    Divides a list of integers in smaller lists of consecutive integers.
    :param data_list: list of integers to divide
    :return: yields the list of consecutive integers
    '''
    for k, g in groupby(enumerate(data_list), lambda i_x: i_x[0] - i_x[1]):
            yield list(map(itemgetter(1), g))

def extract_fasta_from_Uniprot(UniprotID_list, output_file='out_file.fasta'):
    '''
    Given a list of UniprotIDs creates a multifasta output file. The name of the
    output file can be specified when calling the function.
    '''
    with open(output_file, 'w') as out_file:
        for uniprot in UniprotID_list:
            try:
                f = rec.urlopen('https://rest.uniprot.org/uniprotkb/{}.fasta'.format(uniprot))
                fasta = f.read().decode("utf-8", "ignore")
                out_file.write('{}\n'.format(fasta))
            except:
                print('Error with {}'.format(uniprot))

def retrieve_seq_from_uniprot(uniprotID):
    '''
    Given a UniprotID yields its sequence.
    '''
    f = rec.urlopen('https://rest.uniprot.org/uniprotkb/%s.fasta' % (uniprotID))
    fasta = f.read().decode("utf-8", "ignore")
    seq = ''
    for element in fasta.split('\n'):
        if '>' not in element:
            seq += element.strip()

    return seq

def column2list(csvfile, column_name, sep='\t'):
    '''
    Given a csv/tab file and a column name, returns all the column entries
    as an ordened list.
    '''

    data = pd.read_csv(csvfile, sep = sep)
    df = pd.DataFrame(data)

    my_list = df[str(column_name)].tolist()

    return my_list

def mapAccession2Uniprot(id_list, to_db='UniProtKB'):
    '''
    Given a list of protein identifiers, maps them to UniprotIDs.
    :param id_list: list of identifiers to map
    :param from_db: database from which are the identifiers
    :param to_db: database to which identifiers are mapped
    :return mapping_dict: dictionary with initial identifiers as keys and uniprotIDs as values
    '''
    # Uniprot mapping setup
    # Extracted from https://www.uniprot.org/help/id_mapping
    POLLING_INTERVAL = 3
    API_URL = "https://rest.uniprot.org"

    retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
    session = requests.Session()
    session.mount("https://", HTTPAdapter(max_retries=retries))

    def check_response(response):
        try:
            response.raise_for_status()
        except requests.HTTPError:
            print(response.json())
            raise

    def submit_id_mapping(from_db, to_db, ids):
        request = requests.post(
            f"{API_URL}/idmapping/run",
            data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
        )
        check_response(request)
        return request.json()["jobId"]

    def get_next_link(headers):
        re_next_link = re.compile(r'<(.+)>; rel="next"')
        if "Link" in headers:
            match = re_next_link.match(headers["Link"])
            if match:
                return match.group(1)

    def check_id_mapping_results_ready(job_id):
        while True:
            request = session.get(f"{API_URL}/idmapping/status/{job_id}")
            check_response(request)
            j = request.json()
            if "jobStatus" in j:
                if j["jobStatus"] == "RUNNING":
                    print(f"Retrying in {POLLING_INTERVAL}s")
                    time.sleep(POLLING_INTERVAL)
                else:
                    raise Exception(j["jobStatus"])
            else:
                return bool(j["results"] or j["failedIds"])

    def get_batch(batch_response, file_format, compressed):
        batch_url = get_next_link(batch_response.headers)
        while batch_url:
            batch_response = session.get(batch_url)
            batch_response.raise_for_status()
            yield decode_results(batch_response, file_format, compressed)
            batch_url = get_next_link(batch_response.headers)

    def combine_batches(all_results, batch_results, file_format):
        if file_format == "json":
            for key in ("results", "failedIds"):
                if key in batch_results and batch_results[key]:
                    all_results[key] += batch_results[key]
        elif file_format == "tsv":
            return all_results + batch_results[1:]
        else:
            return all_results + batch_results
        return all_results

    def get_id_mapping_results_link(job_id):
        url = f"{API_URL}/idmapping/details/{job_id}"
        request = session.get(url)
        check_response(request)
        return request.json()["redirectURL"]

    def decode_results(response, file_format, compressed):
        if compressed:
            decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
            if file_format == "json":
                j = json.loads(decompressed.decode("utf-8"))
                return j
            elif file_format == "tsv":
                return [line for line in decompressed.decode("utf-8").split("\n") if line]
            elif file_format == "xlsx":
                return [decompressed]
            elif file_format == "xml":
                return [decompressed.decode("utf-8")]
            else:
                return decompressed.decode("utf-8")
        elif file_format == "json":
            return response.json()
        elif file_format == "tsv":
            return [line for line in response.text.split("\n") if line]
        elif file_format == "xlsx":
            return [response.content]
        elif file_format == "xml":
            return [response.text]
        return response.text

    def get_xml_namespace(element):
        m = re.match(r"\{(.*)\}", element.tag)
        return m.groups()[0] if m else ""

    def merge_xml_results(xml_results):
        merged_root = ElementTree.fromstring(xml_results[0])
        for result in xml_results[1:]:
            root = ElementTree.fromstring(result)
            for child in root.findall("{http://uniprot.org/uniprot}entry"):
                merged_root.insert(-1, child)
        ElementTree.register_namespace("", get_xml_namespace(merged_root[0]))
        return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)

    def print_progress_batches(batch_index, size, total):
        n_fetched = min((batch_index + 1) * size, total)
        print(f"Mapped: {n_fetched} / {total} accession numbers to Uniprot")

    def get_id_mapping_results_search(url):
        parsed = urlparse(url)
        query = parse_qs(parsed.query)
        file_format = query["format"][0] if "format" in query else "json"
        if "size" in query:
            size = int(query["size"][0])
        else:
            size = 500
            query["size"] = size
        compressed = (
            query["compressed"][0].lower() == "true" if "compressed" in query else False
        )
        parsed = parsed._replace(query=urlencode(query, doseq=True))
        url = parsed.geturl()
        request = session.get(url)
        check_response(request)
        results = decode_results(request, file_format, compressed)
        total = int(request.headers["x-total-results"])
        print_progress_batches(0, size, total)
        for i, batch in enumerate(get_batch(request, file_format, compressed), 1):
            results = combine_batches(results, batch, file_format)
            print_progress_batches(i, size, total)
        if file_format == "xml":
            return merge_xml_results(results)
        return results

    logging.info('...Mapping protein identifiers to UniprotIDs')
    for from_db in ['RefSeq_Protein','EMBL-GenBank-DDBJ_CDS']:
        job_id = submit_id_mapping(
            from_db=from_db, to_db=to_db, ids=id_list
        )

        if check_id_mapping_results_ready(job_id):
            link = get_id_mapping_results_link(job_id)
            results = get_id_mapping_results_search(link)
            if results['results'] != []:
                break
    i = 0
    mapping_dict = {}
    while i <= len(id_list)-1:
        accession_id = results['results'][int(i)]['from']
        uniprot = results['results'][int(i)]['to']['primaryAccession']
        mapping_dict[accession_id] = uniprot
        i += 1
    return mapping_dict

### 1. OBTAIN PROTEOME FILE ###
def gff_yielder(in_handle, ref_recs):
    """
    Generate SeqRecord and SeqFeatures for gff instances
    """
    for rec in GFF.parse(in_handle, target_lines=1, base_dict=ref_recs):
        yield rec

def obtain_cds_proteins_from_genome(gff_file, genome, current_directory='.'):
    '''
    Creates an output file with all the CDS defined in the gff from the assembly file.
    We obtain a fasta file with all the CDS translated and with the identifier as headers.
    '''
    def protein_recs(gff_file, ref_recs, translation_table='Bacterial'):
        """
        Generate protein records for every CDS in the assembly provided.
        """
        var = False
        i = 7
        with open(gff_file) as in_handle:
            for rec in gff_yielder(in_handle, ref_recs):
                if rec.features != []:
                    for feature in rec.features:
                        if feature.type == 'CDS':
                            if 'pseudo' not in feature.qualifiers.keys():
                                seq_exons = rec.seq[feature.location.nofuzzy_start:feature.location.nofuzzy_end]
                                cds_id = feature.id
                                gene_seq = Seq.Seq(str(reduce(operator.add, seq_exons, "")))
                                if feature.strand == -1:
                                    gene_seq = gene_seq.reverse_complement()
                                protein_seq = pad_seq(gene_seq).translate(table=translation_table)
                                if protein_seq != '':
                                    yield SeqRecord(protein_seq, cds_id + '|' + feature.qualifiers["ID"][0], "", "")

    def pad_seq(seq):
        """
        Pad sequence to multiple of 3 with Ns to avoid conflicts when translating
        """
        return seq + ['', 'NN', 'N'][len(seq) % 3]

    with open(genome) as in_handle:
        ref_recs = SeqIO.to_dict(SeqIO.parse(in_handle, "fasta"))

    base, ext = os.path.splitext(gff_file)
    out_file = current_directory+"/%s-proteins.fa" % base.split('/')[-1]
    with open(out_file, "w") as out_handle:
        SeqIO.write(protein_recs(gff_file, ref_recs), out_handle, "fasta")

    return out_file

### 2. HMMER MODULE FUNCTIONS ###
def extract_df_hmmer(file, hmmerprogram=None, hmmerformat='domtblout'):
    '''
    File converter into a dataframe. By default processes tabulated files, if indicated also processes
    tabulated HMMER output.
    '''
    if hmmerprogram != None:
        if hmmerformat != None:
            df = so.read_hmmer(file, program=hmmerprogram, format=hmmerformat)
            df.columns = ['_'.join(col) for col in df.columns.values]
            return df

def non_overlapping_parser(hmmer_file, hmmerprogram='hmmscan'):
    '''
    When a hmmer output 'domtblout' file is provided, returns a parsed file with only the most significant non-overlapping domains.
    2 core functions:
        - check_overlapping() returns true when there is overlapping domains in a dataframe
        - remove_overlapping() iterates through a dataframe and removes the overlapping domains with lower i-Evalue
    '''
    def check_overlapping(subset):
        # reset df indexing and create and empty dictionary with coordinate numbers
        subset.reset_index(drop=True, inplace=True)
        coord_dict = dict.fromkeys(range(5000+1), [])

        for i in range(len(subset.axes[0])):
            domain = subset.at[i, 'identifier_target_name']
            from_coord = subset.at[i, 'ali_coord_from']
            to_coord = subset.at[i, 'ali_coord_to']
            coords_list = list(range(int(from_coord), int(to_coord) + 1))

            for coord in coords_list:
                if coord_dict[coord] == []:
                    coord_dict[coord] = [domain]
                else:
                    coord_dict[coord].append(domain)

        # iterate in dictionary and return True if overlapping exists
        for key, values in coord_dict.items():
            if len(values) > 1:
                return True

    def remove_overlapping(subset):

        subset_new = pd.DataFrame()
        for i in range(len(subset.axes[0])):
            try:
                subset_not_overlap = pd.DataFrame()
                subset_overlap = pd.DataFrame()

                from_coord = subset.at[i, 'ali_coord_from']
                to_coord = subset.at[i, 'ali_coord_to']
                i_value = subset.at[i, 'this_domain_i-value']  # we'll differentiate same domains for their i-Evalue
                x = set(range(int(from_coord), int(to_coord) + 1))
                row = subset.loc[(subset['this_domain_i-value'] == i_value) & (subset['ali_coord_from'] == from_coord) & (subset['ali_coord_to'] == to_coord)]
                for j in range(len(subset.axes[0])):
                    if i != j:
                        from_coord1 = subset.at[j, 'ali_coord_from']
                        to_coord1 = subset.at[j, 'ali_coord_to']
                        i_value1 = subset.at[j, 'this_domain_i-value']
                        y = range(int(from_coord1), int(to_coord1))
                        inter = x.intersection(y)
                        row1 = subset.loc[(subset['this_domain_i-value'] == i_value1) & (subset['ali_coord_from'] == from_coord1) & (subset['ali_coord_to'] == to_coord1)]
                        if inter != set():
                            subset_overlap = pd.concat([subset_overlap, row1])
                            subset_overlap = pd.concat([subset_overlap, row]).drop_duplicates()

                        else:
                            subset_not_overlap = pd.concat([subset_not_overlap, row]).drop_duplicates()

                if subset_overlap.empty == False:
                    e_value_list = [float(x) for x in subset_overlap['this_domain_i-value']]
                    e_value = min(e_value_list)
                    if e_value_list.count(e_value) > 1:
                        subset1 = subset_overlap.loc[subset_overlap['this_domain_i-value'] == str(e_value)]
                        score_list = [float(x) for x in subset1['this_domain_score']]
                        max_score = max(score_list)
                        min_row = subset1.loc[(subset1['this_domain_i-value'] == str(e_value)) & (subset1['this_domain_score'] == str(max_score))]
                    else:
                        min_row = subset_overlap.loc[subset_overlap['this_domain_i-value'] == str(e_value)]
                    subset_new = pd.concat([subset_new, min_row]).drop_duplicates()

                flt = ['this_domain_i-value', 'ali_coord_from', 'ali_coord_to']
                if ((subset_not_overlap[flt] == row[flt]).all(axis=1)).sum() != 0 and subset_overlap.empty == True:
                    subset.drop(i, inplace=True)
                    subset_new = pd.concat([subset_new, subset_not_overlap])

                subset.reset_index(drop=True, inplace=True)
            except:
                break
        return subset_new

    df = extract_df_hmmer(hmmer_file, hmmerprogram=hmmerprogram, hmmerformat='domtblout')
    hmmer_parsed = pd.DataFrame()

    # make a unique ordered list with all the queries for every TF
    queries = set()
    queries_list = [x for x in df['identifier_query_name']]
    queries_list = [x for x in queries_list if x not in queries and queries.add(x) is None]

    for query in queries_list:
        subset = df.loc[df['identifier_query_name'] == query]  # make subsets for every TF (process whole hmmer file
        subset.reset_index(drop=True, inplace=True)

        while check_overlapping(subset) == True:
            subset = remove_overlapping(subset)

        hmmer_parsed = pd.concat([hmmer_parsed, subset])

    hmmer_parsed.reset_index(drop=True, inplace=True)
    hmmer_parsed.to_csv('{}_output_parsed.tsv'.format(hmmerprogram), sep='\t', index=False)

    return '{}_output_parsed.tsv'.format(hmmerprogram)

def retrieve_domains_from_hmmscan(file, outfile='hmmfetch_keyfile.txt'):
    '''
    Retrieve the domains key from a hmmscan parsed output file.
    Produces a input file for hmmfetch without repeated domains.
    '''

    # df = extract_df_hmmer(file, hmmerprogram='hmmscan', hmmerformat='domtblout')

    domains = set()
    # domains_list = [x for x in df['identifier_target_name']]
    domains_list = column2list(file, 'identifier_target_name')
    domains_list = [x for x in domains_list if x not in domains and domains.add(x) is None]

    out_file = open(outfile, 'w')
    for i in domains_list:
        out_file.write('%s\n' %(i))
    out_file.close()

    return outfile

def hmmsearch_parser(hmmsearch_file, hmmscan_file, db_tsv):
    '''
    When provided with a database .tsv file and output hmmscan and hmmsearch files, parses them and merges them
    into a unique file with all the information.
    '''

    # input file processing
    tsv_datafile = pd.read_csv(db_tsv, sep='\t')
    db_df = pd.DataFrame(tsv_datafile)

    scan_df = pd.read_csv(hmmscan_file, sep='\t')
    search_df = extract_df_hmmer(hmmsearch_file, hmmerprogram='hmmsearch')

    # Dict with alias and uniprot
    alias_dict = {}
    for alias in db_df['Alias']:
        s = db_df.loc[db_df['Alias'] == alias, 'UniprotID']
        alias_dict[alias] = s.to_string(index=False).strip(' ')

    # Dict with uniprotID and PFAM domains + their e-value
    uni_pfam_dict = {}
    for line in scan_df['identifier_query_name']:
        uniprot = line.split('|')[1]
        uni_pfam_dict.setdefault(uniprot, [])

        domain_list = list(scan_df.loc[scan_df['identifier_query_name'] == line, 'identifier_target_name'])
        e_value_list = list(scan_df.loc[scan_df['identifier_query_name'] == line, 'this_domain_i-value'])
        zipped = list(zip(domain_list, e_value_list))
        uni_pfam_dict[uniprot] = zipped

    # Dict with domains and their description
    domain_description = {}
    for domain in scan_df['identifier_target_name']:
        s = scan_df.loc[scan_df['identifier_target_name'] == domain, 'identifier_query_description']
        domain_description[domain] = ((s.to_string(index=False)).split('\n')[0]).strip(' ')

    # Dict with pfam domain and each hit of acnes + stats + description
    pfam_acnes_dict = {}
    for domain in search_df['identifier_query_name']:
        pfam_acnes_dict.setdefault(str(domain).strip(), [])
        acnes_targets = search_df.loc[search_df['identifier_query_name'] == domain, 'identifier_target_name']
        acnes_targets = [str(x.split('|')[0]).strip('cds-') for x in acnes_targets]
        coords1 = search_df.loc[search_df['identifier_query_name'] == domain, 'ali_coord_from']
        coords2 = search_df.loc[search_df['identifier_query_name'] == domain, 'ali_coord_to']
        coordinates = list(zip(coords1, coords2))
        zipped0 = list(zip(acnes_targets, coordinates)) # target + coordinates

        d = search_df.loc[search_df['identifier_query_name'] == domain]
        acnes_stats = (d.to_string(index=False, header=False)).split('\n')
        acnes_stats = ['{},{}'.format(x.split()[6], x.split()[11:13]) for x in acnes_stats]
        acnes_description = list(search_df.loc[search_df['identifier_query_name'] == domain, 'identifier_query_description'])
        acnes_description = [phrase.translate({ord(i):'\t' for i in ','}) for phrase in acnes_description] # remove the ',' from the description
        zipped1 = list(zip(acnes_stats, acnes_description)) # stats + description
        zipped2 = list(zip(zipped0, zipped1))

        pfam_acnes_dict[domain] = zipped2

    domain_not_found = []
    alias_wo_domain = []
    # Formatting the output file
    with open('hmmsearch_output_parsed.tsv', 'w') as out_file:
        out_file.write('Alias\tUniprotID\tPFAM_domain\tE-value_domain\tDomain_description\tTarget_name\tFrom_coordinate\tTo_coordinate\tE-value_full_sequence\tc-Evalue_domain\ti-Evalue_domain\tDescription\n')

        for alias, uniprot in alias_dict.items():
            if uniprot not in uni_pfam_dict:
                alias_wo_domain.append((alias, uniprot.split('\n')[0]))
            else:
                a = uni_pfam_dict[uniprot]
                for i in a:
                    domain = str(i[0]).translate({ord(e): None for e in "''[]"})
                    if domain not in pfam_acnes_dict:
                        domain_not_found.append(domain)
                    else:
                        for value in pfam_acnes_dict[domain]:
                            out_file.write('%s\t' %(alias))
                            out_file.write('%s\t' %(uniprot))
                            out_file.write('%s\t' %(domain))
                            out_file.write('%s\t' %(str(i[1]).translate({ord(e):None for e in "''[]"})))
                            out_file.write('%s\t' %(domain_description[domain]))
                            out_file.write('%s\n'%((str(value).translate({ord(i):None for i in "(())[]\"\"''"})).translate({ord(i):'\t' for i in ','}))) # avoid the dict formatting and separate in tab
    logging.info('hmmscan did not find any domain in the following transcription factors: {}'.format(set(alias_wo_domain)))
    logging.info('hmmsearch did not find any of the following domains in the provided proteome: {}'.format(set(domain_not_found)))
    return 'hmmsearch_output_parsed.tsv'

def only_hits_with_all_domains(hmmsearch_file):
    '''
    Returns a dataframe with only the hmmsearch targeted proteins that contain all the domains of the initial TF.
    :param hmmsearch_file:
    :return dataframe of the hits with all domains in the reference TF:
    '''

    # process the hmmer file
    data = pd.read_csv(hmmsearch_file, sep='\t')
    hmmer_df = pd.DataFrame(data)

    # there can be alias shared between different uniprots so we need a more distinctive identifier
    hmmer_df['Alias'] = hmmer_df['Alias'].str.strip() + '|' + hmmer_df['UniprotID'].str.strip()

    # make a unique ordered list with all the queries for every TF
    queries = set()
    queries_list = [x for x in hmmer_df['Alias']]
    queries_list = [x for x in queries_list if x not in queries and queries.add(x) is None]

    subset_new = pd.DataFrame()

    for query in queries_list:
        subset = hmmer_df.loc[hmmer_df['Alias'] == query]
        alias_dict = {}
        for alias in subset['Alias']:
            domains = set(list(subset.loc[hmmer_df['Alias'] == alias, 'PFAM_domain']))
            hits = [(domain, set(subset.loc[hmmer_df['PFAM_domain'] == domain, 'Target_name'])) for domain in domains]
            alias_dict.setdefault(alias, hits)

        alias_keys = [x for x in alias_dict.keys()]

        for alias in alias_keys:
            count_list = []
            n_domains = len(alias_dict[alias])
            for i in range(n_domains):
                count_list.append(list(alias_dict[alias][i][1]))

            counts = Counter(list(chain.from_iterable(count_list)))

            for uniprot, n in counts.items():
                if n == n_domains:
                    row = subset.loc[subset['Target_name'] == uniprot]
                    subset_new = pd.concat([subset_new, row])
    subset_new.reset_index(inplace=True)
    subset_new.drop(labels='index', axis=1, inplace=True)
    return subset_new

### 3. TM-ALIGN MODULE FUNCTIONS ###
def yield_PDB(pdbfile, CA = False):
    '''
    Generator function to yield every pdb file line.
    :param pdbfile: pdb filepath
    :param CA: when flagged true, only yields the carbon alpha lines
    '''

    if CA == True:
        with open(pdbfile, 'r') as fd:
            for line in fd:
                if re.search(r'^ATOM\s+\d+\s+CA\s+', line):
                    yield line
    else:
        with open(pdbfile, 'r') as fd:
            for line in fd:
                yield line

def get_alphafold_structure(uniprotID, outfilename = None):
    '''
    Retrieve structure with UniprotID from alphafold database.
    '''

    url = ('https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb'.format(uniprotID))
    if outfilename == None:
        outfilename = uniprotID
    try:
        req = rec.urlretrieve(url, 'AF-{}.pdb'.format(outfilename))
        return 'AF-{}.pdb'.format(outfilename)
    except:
        logging.info('No structure available for {}'.format(uniprotID))
        return None

def check_pLDDT(pdbfile, protein_percentage = 10):
    pdbCA = []
    CA_dict = {}
    for element in yield_PDB(pdbfile, CA = True):
        pdbCA.append(element.split()[5])
        if float(element.split()[10]) < 70:
            CA_dict[int(element.split()[5])] = { 'aa': element.split()[3],
                                                'pLDDT': element.split()[10]
                                                }
    # Define the threshold for considering long coils (in base of the protein percentage provided)
    threshold = (protein_percentage * len(pdbCA)) / 100
    if len(str(threshold).split('.')) > 1:
        if str(threshold)[-1] >= '5':
            threshold = int(threshold) + 1
        else:
            threshold = int(threshold)

    # identify the protein fragments with low pLDDT that represent the same or higher percentage of the protein
    remove_list = []
    for fragment in divide_list_consecutive_numbers(CA_dict.keys()):
        if len(fragment) >= threshold:
            remove_list.append(fragment)

    remove_list_flat = list(chain(*remove_list))

    # rewrite the pdb file
    with open('{}'.format(pdbfile.strip('AF-')), 'w') as outfile:
        for line in yield_PDB(pdbfile):
            if re.search(r'^ATOM', line):
                if int(line.split()[5]) not in remove_list_flat:
                    outfile.write(line)
            else:
                outfile.write(line)

def get_alphafolds(hmmer_df, pwd='.'):
    
    logging.info('...Downloading the AlphaFold structures from UniprotIDs')
    # set working directory the folder in which all the TF folders with the alphafolds will be created
    os.chdir(pwd)

    # make dict with the reference TF and the putative proteins to be superimposed 
    superimposition_dict = {}
    aliases = set(hmmer_df['Alias'].to_list())
    for alias in aliases: 
        superimposition_dict[alias.split('|')[1]] = {
            'targets': list(hmmer_df.loc[hmmer_df['Alias'] == alias, 'Target_uniprot']),
            'alias': alias
            }

    # download all structures to be superimposed
    ref_uniprots = set(hmmer_df['UniprotID'].to_list())
    if os.path.isdir('reference_alphafolds') == False:
        os.mkdir('reference_alphafolds')
    os.chdir('reference_alphafolds')
    # continue where it was left off
    already_uniprots = [f.strip('AF-ref').strip('.pdb') for f in os.listdir('.') if f.endswith('.pdb')]

    for uniprot in ref_uniprots:
        if uniprot not in already_uniprots: 
            refpdb = get_alphafold_structure(uniprot, outfilename='ref{}'.format(uniprot))
            if refpdb == None: 
                # remove from the superimposition_dict the TF without structure available
                superimposition_dict.pop(uniprot)
            else:
                check_pLDDT(refpdb)
    os.system('rm AF-*')
    os.chdir('../')

    target_uniprots = set(hmmer_df['Target_uniprot'].to_list())
    if os.path.isdir('target_alphafolds') == False:
        os.mkdir('target_alphafolds')
    os.chdir('target_alphafolds')
    # continue where it was left off
    already_uniprots = [f.strip('AF-').strip('.pdb') for f in os.listdir('.') if f.endswith('.pdb')]
    for uniprot in target_uniprots:
        if uniprot not in already_uniprots: 
            pdb = get_alphafold_structure(uniprot)
            if pdb == None:
                superimposition_dict1 = superimposition_dict.copy()
                for key, values in superimposition_dict1.items():
                    if uniprot in values['targets']:
                        # remove from the superimposition_dict the target protein without structure available
                        temp = list(set(values['targets']))
                        temp.remove(uniprot)
                        superimposition_dict[key] = {
                            'targets': temp,
                            'alias': values['alias']
                        }
            else:
                check_pLDDT(pdb)
    os.system('rm AF-*')
    os.chdir('../')
    logging.info('All AlphaFold structures downloaded')
    return superimposition_dict

def TM_align(superimposition_dict, TMAlign_folder='.'):
    '''
    Relies in a local installation of TM-Align (https://zhanggroup.org/TM-align/).
    When provided a list with the names of the folders created after executing the get_alphafolds() function.
    In each folder the TM-Align pairwise comparisons are performed for each hit and the reference TF protein.
    The output of all comparisons along their TM-scores are merged in a file '{TF_name}.out'.
    '''
    logging.info('...Performing pair-wise structural superimpositions')
    os.chdir(TMAlign_folder)
    for refuniprot, values in superimposition_dict.items():
        if '{}.txt'.format(refuniprot) not in os.listdir('.'):
            for target in values['targets']:
                os.system('TMalign reference_alphafolds/ref{}.pdb target_alphafolds/{}.pdb -o {}.sup > {}.out'.format(refuniprot, target, target, target))
            os.system('cat *.out > {}.txt'.format(refuniprot))
            os.system('rm *sup*')
            os.system('rm *out')

    logging.info('All {} structural superimpositions performed.'.format(len(superimposition_dict)))

    return superimposition_dict

def TM_align_parser(superimposition_dict, TMAlign_threshold=0.5, summary_each_TF=False):
    '''
    Filters all pair-wise comparisons made by TM-Align based on a threshold for the TM-score (default=0.5).
    Provides a summary tsv with the best alignments and statistics for each TF.
    '''
    with open('best_alignments.tsv', 'w') as outfile:
        outfile.write('TF\tTarget_uniprot\tTM_score\tRMSD\tRatio_aaIdentical/aaAligned\n')
    
    for uniprotTF, values in superimposition_dict.items():
        count_alignments = 0
        with open('{}.txt'.format(uniprotTF), 'r') as f:
            lines = f.read().splitlines()
            tm_dict = {}
            for i in range(len(lines)):
                if lines[i].startswith(' *') and lines[i+6].startswith('Name of Chain_1:'):
                    uniprot_acnes = str(lines[i+7].split(' ')[3]).strip('.pdb').strip('target_alphafolds/')
                    RMSD = str(lines[i+11].split(',')[1]).strip(' RMSD=').strip()
                    ratio_identical_aligned = lines[i+11].split(',')[2].strip(' Seq_ID=n_identical/n_aligned= ').strip()
                    TM_score = lines[i+12].split(' ')[1]
                    tm_dict[uniprot_acnes] = { 'TM_score': float(TM_score),
                                                'RMSD': RMSD,
                                                'ratio': ratio_identical_aligned,
                                             }

                    count_alignments += 1
        better_alignments = [k.strip('target_alphafolds/') for k, d in tm_dict.items() if d['TM_score'] >= float(TMAlign_threshold)]
        worst_alignments = [a for a in tm_dict.keys() if a not in better_alignments]

        if summary_each_TF==True:
            with open('{}_score_summary.out'.format(values['alias']), 'w') as outfile:
                outfile.write('Summary of pair-wise alignments of {}\n\n'.format(values['alias']))
                outfile.write('Best pair-wise alignments: {}\n'.format(len(better_alignments)))

                if len(better_alignments) != 0:
                    outfile.write('Uniprot_acnes\tTM_score\tRMSD\tRatio_aaIdentical/aaAligned\n')
                for element in better_alignments:
                    outfile.write('{}\t'.format(element))
                    outfile.write('{}\t{}\t{}\t'.format(tm_dict[element]['TM_score'], tm_dict[element]['RMSD'], tm_dict[element]['ratio']))
                    outfile.write('\n')

                outfile.write('\nOther pair-wise alignments: {}\n'.format(len(worst_alignments)))
                if len(worst_alignments) != 0:
                    outfile.write('Target_uniprot\tTM_score\tRMSD\tRatio_aaIdentical/aaAligned\n')
                for element in worst_alignments:
                    outfile.write('{}\t'.format(element))
                    outfile.write('{}\t{}\t{}\t'.format(tm_dict[element]['TM_score'], tm_dict[element]['RMSD'], tm_dict[element]['ratio']))
                    outfile.write('\n')
                outfile.write('\nNumber of pair-wise alignments performed: {}'.format(count_alignments))

        with open('best_alignments.tsv', 'a') as outfile:
            for element in better_alignments:
                outfile.write('{}\t'.format(values['alias']))
                outfile.write('{}\t'.format(element))
                outfile.write('{}\t{}\t{}\t'.format(tm_dict[element]['TM_score'], tm_dict[element]['RMSD'],tm_dict[element]['ratio']))
                outfile.write('\n')
    os.system('rm *.txt')
    return os.getcwd()+'/best_alignments.tsv'


### 4. GENOMIC CONTEXT MODULE FUNCTIONS ###
# 4.1 Obtain a dictionary with all the targets and their correspondant genomic files
def retrieve_targets_dict(hmmer_df, path2gff= '.', path2fna='.', TM_Align_file=None, executable_path='.'):
    '''
    Parses the provided file into a dictionary with the identifiers of the target proteins and
    the reference proteins to which have mapped.
    :param targets_file: the file of targets to be parsed into the dictionary
    :param TM_Align: if True, takes the output file of executing the TM_Align_module, if set False
                    takes the dataframe of the HMMER_module
    '''
    dbdf = pd.read_csv(executable_path+'TF_databases/TF_database_merged.tsv', sep='\t')

    def download_genomic_files(dict2download):
        logging.info('GC module: downloading genomic files to plot genomic context...')
        downloaded_dict = {}
        if os.path.isdir('downloaded_files'):
            os.chdir('downloaded_files')
        else:
            os.mkdir('downloaded_files')
            os.chdir('downloaded_files')
        currentdir = os.getcwd()
        for key, values in dict2download.items():
            if 'https' in values['link']:
                link = values['link'] +'/'+ values['link'].split('/')[-1]
                if link.split('/')[-1]+'_genomic.gff' in os.listdir():
                    if os.path.exists(currentdir+'/'+link.split('/')[-1]+'_genomic.fna'):
                        fna = currentdir+'/'+link.split('/')[-1]+'_genomic.fna'
                        faa = None
                    pass
                else:
                    url = link+'_genomic.gff.gz'
                    req = rec.urlretrieve(url, '{}'.format(url.split('/')[-1]))
                    url = link+'_genomic.fna.gz'
                    req = rec.urlretrieve(url, '{}'.format(url.split('/')[-1]))
                    faa = None
                    fna = currentdir+'/'+link.split('/')[-1]+'_genomic.fna'
                    os.system('gunzip *.gz')

                downloaded_dict[key] = {
                'uniprot': values['uniprot'],
                'gff': currentdir+'/'+link.split('/')[-1]+'_genomic.gff',
                'fna': fna,
                'faa': faa
                }
            else:
                downloaded_dict[key] = {
                'uniprot': values['uniprot'],
                'gff': executable_path+values['link'].split('|')[0].strip('../'),
                'fna': executable_path+values['link'].split('|')[1].strip('../'),
                'faa': None
                }
        os.chdir('../')
        return downloaded_dict

    dict2download = {}
    targets_dict = {}
    tmp = {}
    if TM_Align_file != None:
        data = pd.read_csv(TM_Align_file, sep='\t', index_col=False)
        df = pd.DataFrame(data)

        df['UniprotID'] = df['TF'].apply(lambda x: x.split('|')[1])
        uniprots = set(df['UniprotID'].to_list())
        for uniprot in uniprots:
            line = dbdf.loc[dbdf['UniprotID'] == uniprot, ['Accession_mapped','Assembly_link']]
            accession = list(line['Accession_mapped'])[0]
            link = list(line['Assembly_link'])[0]
            # retrieve the links to download the genomic files for reference TF
            dict2download[accession] = {
            'uniprot': uniprot,
            'link': link
            }
            # dictionary of the genomic contexts to compare
            target_uniprots = list(df.loc[df['UniprotID'] == uniprot, 'Target_uniprot'])
            TF = list(df.loc[df['UniprotID'] == uniprot, 'TF'])[0]
            mapper_dict = dict(zip(hmmer_df.Target_uniprot, hmmer_df.Target_name))
            target_accessions = []
            for uni in target_uniprots:
                acc = mapper_dict[uni]
                target_accessions.append(acc)
            tmp[(accession,TF)] = target_accessions

    else:
        uniprots = set(hmmer_df['UniprotID'].to_list())
        for uniprot in uniprots:
            line = dbdf.loc[dbdf['UniprotID'] == uniprot, ['Accession_mapped','Assembly_link']]
            accession = list(line['Accession_mapped'])[0]
            link = list(line['Assembly_link'])[0]
            # retrieve the links to download the genomic files for reference TF
            dict2download[accession] = {
                'uniprot': uniprot,
                'link': link
            }
            # dictionary of the genomic contexts to compare 
            target_accessions = list(hmmer_df.loc[hmmer_df['UniprotID'] == uniprot, 'Target_name'])
            TF = list(hmmer_df.loc[hmmer_df['UniprotID'] == uniprot, 'Alias'])[0]
            tmp[(accession,TF)] = target_accessions
        
    downloaded_dict = download_genomic_files(dict2download)

    targets_dict1 = tmp.copy()
    for key, values in targets_dict1.items():
        accession, TF = key
        tmp_dict = {}
        for value in values:
            tmp_dict[value] = {
            'gff': path2gff,
            'fna': path2fna
            }

        targets_dict[accession] = {
        'TF': TF,
        'gff': downloaded_dict[accession]['gff'],
        'fna': downloaded_dict[accession]['fna'],
        'faa': downloaded_dict[accession]['faa'],
        'targets': tmp_dict
        }
    del downloaded_dict, targets_dict1, tmp, tmp_dict

    return targets_dict

# 4.2 Routines to obtain the genomic contexts from the local gneomic files
# Functions modified from GCsnap (https://github.com/JoanaMPereira/GCsnap) - developed by Joana Pereira, Structural Computational Biology, Biozentrum University of Basel, Basel Switzerland
def parse_local_assembly(assembly_file, get_chunk=False, chunk_size=0, target=None):
    assembly = {'ncbi_codes': [], 'starts': [], 'ends': [], 'directions': [], 'names': [], 'scaffolds': []}
    curr_scaffold = 0
    found_chunk = False
    chunk_timer = 0

    with open(assembly_file, 'r') as in_assembly:
        for line in in_assembly:
            if not line.startswith('#'):
                line_data = line.split('\t')
                if line_data != ['\n']:
                    if line_data[2] == 'CDS' and 'pseudo=' not in line:

                        start = int(line_data[3])
                        end = int(line_data[4])
                        direction = line_data[6]

                        if 'cds-' in line:
                            ncbi_code = line_data[8].split('ID=cds-')[1].split(';')[0]
                        else:
                            if 'Name=' in line:
                                ncbi_code = line_data[8].split('Name=')[1].split(';')[0]
                            else:
                                ncbi_code = 'unk'

                        if 'pseudo=' not in line and 'product=' in line and 'fragment' not in line:
                            prot_name = line_data[8].split('product=')[1].split(';')[0]
                        else:
                            prot_name = 'pseudogene'

                        if ncbi_code in assembly['ncbi_codes'] and assembly['ncbi_codes'][
                            -1] == ncbi_code:  # it means that this is some kind of fragmented gene (has introns?...) and so we have to collect the largest interval encompassing it
                            if start < assembly['starts'][-1]:
                                assembly['starts'][-1] = start
                            if end > assembly['ends'][-1]:
                                assembly['ends'][-1] = end
                        else:
                            if '|' in ncbi_code:
                                ncbi_code = ncbi_code.replace('|', '_')

                            assembly['ncbi_codes'].append(ncbi_code)
                            assembly['names'].append(prot_name)
                            assembly['scaffolds'].append(curr_scaffold)
                            assembly['starts'].append(start)
                            assembly['ends'].append(end)
                            assembly['directions'].append(line_data[6])

                        if get_chunk:
                            if ncbi_code == target:
                                chunk_timer = 1
                            elif chunk_timer > 0:
                                chunk_timer += 1

                            if chunk_timer == round(chunk_size / 2):
                                break

            elif line.startswith('##sequence-region'):
                curr_scaffold += 1

    if get_chunk:
        chunked_assembly = {}
        for key in assembly:
            chunked_assembly[key] = assembly[key][-chunk_size:]
        assembly = chunked_assembly

    return assembly

def get_n_flanking_genes(target_ncbi_code, assembly, n_5 = None, n_3 = None):

    flanking_genes = {'relative_starts': [], 'relative_ends': []}

    if target_ncbi_code in assembly['ncbi_codes']:
        index_of_target = assembly['ncbi_codes'].index(target_ncbi_code)
        direction_of_target = assembly['directions'][index_of_target]

        if direction_of_target == '+':
            genomic_context_block = [index_of_target-n_5, index_of_target+n_3]
        else:
            genomic_context_block = [index_of_target-n_3, index_of_target+n_5]

        for i in range(genomic_context_block[0], genomic_context_block[1]+1):
            if i>= 0 and i < len(assembly['scaffolds']):
                if assembly['scaffolds'][i] == assembly['scaffolds'][index_of_target]:
                    for key in assembly.keys():
                        if key != 'scaffolds':
                            if key not in flanking_genes:
                                flanking_genes[key] = []
                            flanking_genes[key].append(assembly[key][i])

                            if key == 'starts' or key == 'ends':
                                flanking_genes['relative_{}'.format(key)].append(assembly[key][i] - assembly['starts'][index_of_target] + 1)

        print(' ...Found {} flanking genes for {}'.format(len(flanking_genes['ncbi_codes'])-1, target_ncbi_code))

        if direction_of_target == '-':
            index_of_target_in_flanking = flanking_genes['ncbi_codes'].index(target_ncbi_code)
            current_relative_starts = flanking_genes['relative_starts']
            current_relative_ends = flanking_genes['relative_ends']

            flanking_genes['relative_starts'] = [value*(-1)+current_relative_ends[index_of_target_in_flanking]+1 for value in current_relative_ends]
            flanking_genes['relative_ends'] = [value*(-1)+current_relative_ends[index_of_target_in_flanking]+1 for value in current_relative_starts]

            for key in flanking_genes:
                flanking_genes[key] = flanking_genes[key][::-1]
                if key == 'directions':
                    flanking_genes[key] = list(''.join(flanking_genes[key]).replace('+','p'))
                    flanking_genes[key] = list(''.join(flanking_genes[key]).replace('-','+'))
                    flanking_genes[key] = list(''.join(flanking_genes[key]).replace('p','-'))

        return flanking_genes

    else:
        print(' ...{} was not found in assembly'.format(target_ncbi_code))
        return 'nan'

def add_sequences_to_flanking_genes_from_local_assembly(flanking_genes, target_ncbi_code, protein_file):

    flanking_genes['sequences'] = []
    flanking_genes['species'] = []
    flanking_genes['taxID'] = []

    # obtain a dictionary with protein identifiers (ncbi_codes) as keys and their sequences as values
    proteins_dict = {}
    for seq_record in SeqIO.parse(protein_file, 'fasta'):
        if 'cds' in seq_record.id:
            proteins_dict[str(seq_record.id.split('|')[0].strip('cds-'))] = str(seq_record.seq).strip('*')
        elif '.prt' in protein_file:
            proteins_dict[str(seq_record.id.split()[0])] = str(seq_record.seq).strip('*')
        else:
            proteins_dict[str(seq_record.id.split()[0].strip('>'))] = str(seq_record.seq).strip('*')

    for ncbi_code in flanking_genes['ncbi_codes']:
        seq = ''
        species = ''
        uniprot_code = ''
        swiss_model_link = ''
        taxid = ''
        # retrieve the sequence from the CDS_proteins file with the ncbi_code (AAT81777.1)
        ### MODIFICATION 18/04 ###

        flanking_genes['sequences'].append(proteins_dict[ncbi_code])

        if species != '':
            flanking_genes['species'] = species
        if taxid != '':
            flanking_genes['taxID'] = taxid

    return flanking_genes

# 4.3 Function to retrieve the genomic contexts of all target proteins in target_dict
def genomic_context(target_dict, n_flanking5=4, n_flanking3=4):

    all_syntenies = {}
    for reference in target_dict.keys():
        all_syntenies[reference] = {}

        if target_dict[reference]['faa'] != None:
            proteome = target_dict[reference]['faa']
        else:
            # obtain proteome from gff and fna
            proteome = obtain_cds_proteins_from_genome(target_dict[reference]['gff'], target_dict[reference]['fna'])
        # parse gff
        assembly = parse_local_assembly(target_dict[reference]['gff'])
        # extract flanking genes from gff
        flanking_genes = get_n_flanking_genes(reference, assembly, n_5=n_flanking5, n_3=n_flanking3)
        # add protein sequences to flanking genes from obtained proteome
        genomic_context = add_sequences_to_flanking_genes_from_local_assembly(flanking_genes, reference, proteome)

        all_syntenies[reference]['flanking_genes'] = genomic_context
        all_syntenies[reference]['assembly_id'] = reference

        tmp_dict = target_dict[reference]['targets']
        for target in tmp_dict.keys():
            all_syntenies[target] = {}
            # obtain proteome from gff and fna
            proteome = obtain_cds_proteins_from_genome(tmp_dict[target]['gff'], tmp_dict[target]['fna'])
            # parse gff
            assembly = parse_local_assembly(tmp_dict[target]['gff'])
            # extract flanking genes from gff
            flanking_genes = get_n_flanking_genes(target, assembly, n_5=n_flanking5, n_3=n_flanking3)
            # add protein sequences to flanking genes from obtained proteome
            genomic_context = add_sequences_to_flanking_genes_from_local_assembly(flanking_genes, target, proteome)

            all_syntenies[target]['flanking_genes'] = genomic_context
            all_syntenies[target]['assembly_id'] = target

    return all_syntenies

# 4.4 Routines to find and add the protein families to every gene found
# Functions adapted from GCsnap (https://github.com/JoanaMPereira/GCsnap) - developed by Joana Pereira, Structural Computational Biology, Biozentrum University of Basel, Basel Switzerland
def compute_all_against_all_distance_matrix(in_syntenies, out_label = None, num_threads = None, num_alignments = None, max_evalue = None, num_iterations = None, blast = None, mmseqs = None, default_base = None, tmp_folder = None, mode = 'flanking_sequences', method = 'psiblast', min_coverage=None):

    def write_flanking_sequences_to_fasta(all_syntenies, out_dir, out_label, exclude_pseudogenes = False, mode = 'flanking_sequences'):

        out_fasta = '{}/{}_{}.fasta'.format(out_dir, out_label, mode)
        seqs_lens = {}

        with open(out_fasta, 'w') as outfst:
            for target in all_syntenies.keys():

                if mode == 'flanking_sequences':
                    try:
                        flanking_genes = all_syntenies[target]['flanking_genes']
                        for i, ncbi_code in enumerate(flanking_genes['ncbi_codes']):

                            if flanking_genes['names'][i] != 'pseudogene' or not exclude_pseudogenes:
                                outfst.write('>{}|{}\n{}\n'.format(ncbi_code, flanking_genes['names'][i], flanking_genes['sequences'][i]))
                                seqs_lens[ncbi_code] = len(flanking_genes['sequences'][i])
                    except:
                        logging.info('No flanking genes found for {}'.format(target))

                if mode == 'operon':
                    outfst.write('>{}\n'.format(target))
                    for sequence in all_syntenies[target]['flanking_genes']['sequences']:
                        outfst.write(sequence)
                    outfst.write('\n')

        return out_fasta, seqs_lens

    ### Routines for executing local psiblast to find the families
    def make_blast_database_from_fasta(infasta, blast = None):

        blastDB = '{}_blastDB'.format(infasta[:-6])

        make_blastDB = sp.Popen(['makeblastdb','-in', infasta, '-dbtype', 'prot', '-out', blastDB], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = make_blastDB.communicate()

        if len(stderr) > 0:
            print(stderr)
            if 'BLAST engine error' not in stderr:
                make_blastDB = sp.Popen(['{}makeblastdb'.format(blast),'-in', infasta, '-dbtype', 'prot', '-out', blastDB], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
                stdout, stderr = make_blastDB.communicate()
        return blastDB

    def run_blast_for_flanking_sequences(seq_fasta, database, num_threads = None, num_alignments = None, max_evalue = None, num_iterations = None, blast = None, tmp_folder = None):

        print(' ... ...Running BLAST')
        blast_outfile = '{}/{}_{}.xml'.format(tmp_folder, seq_fasta.split('/')[-1][:-6], max_evalue)

        run_blast = sp.Popen(['psiblast', '-query', seq_fasta, '-db', database, '-out', blast_outfile, '-num_threads', str(num_threads), '-evalue', str(max_evalue), '-inclusion_ethresh', str(max_evalue), '-num_iterations', str(num_iterations),'-num_alignments', str(num_alignments), '-outfmt', '5'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_blast.communicate()

        if len(stderr) > 0:
            print(stderr)
            if 'BLAST engine error' not in stderr:
                run_blast = sp.Popen(['{}psiblast'.format(blast), '-query', seq_fasta, '-db', database, '-out', blast_outfile, '-num_threads', str(num_threads), '-evalue', str(max_evalue), '-inclusion_ethresh', str(max_evalue), '-num_iterations', str(num_iterations),'-num_alignments', str(num_alignments), '-outfmt', '5'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
                stdout, stderr = run_blast.communicate()

        return blast_outfile

    def extract_distance_matrix_from_blast_output(blast_results, default_base = None, min_coverage = 70, sequences_lengths = {}):

        def merge_intervals(intervals):

            intervals = np.array(intervals)
            intervals = intervals[intervals[:, 0].argsort()]

            new_intervals = []
            for interval in intervals:
                if len(new_intervals) == 0:
                    new_intervals.append(interval)
                else:
                    previous_interval = list(new_intervals[-1])
                    if interval[0] <= previous_interval[1]:
                        overlap_range = interval+previous_interval
                        new_intervals[-1] = [min(overlap_range), max(overlap_range)]
                    else:
                        new_intervals.append(interval)

            return new_intervals

        result_handle = open(blast_results)
        blast_records = NCBIXML.parse(result_handle)

        all_queries = sorted(list(set([record.query.split('|')[0] for record in blast_records])))
        distance_matrix = [[default_base for query in all_queries] for query in all_queries]

        result_handle = open(blast_results)
        blast_records = NCBIXML.parse(result_handle)

        for record in blast_records:
            query_name = record.query.split('|')[0]
            query_index = all_queries.index(query_name)

            for alignment in record.alignments:
                target_name = alignment.title.split('|')[2].split()[1]
                target_index = all_queries.index(target_name)

                if query_name != target_name:
                    query_intervals = []
                    sbjct_intervals = []

                    for hsp in alignment.hsps:
                        query_intervals.append(np.array([hsp.query_start, hsp.query_end]))
                        sbjct_intervals.append(np.array([hsp.sbjct_start, hsp.sbjct_end]))

                    query_intervals  = merge_intervals(query_intervals)
                    target_intervals = merge_intervals(sbjct_intervals)

                    if query_name in sequences_lengths and target_name in sequences_lengths:
                        query_length  = sequences_lengths[query_name]
                        target_lenght = sequences_lengths[target_name]

                        query_coverage  = sum([i[-1]-i[0] for i in query_intervals])*100.0/float(query_length)
                        target_coverage = sum([i[-1]-i[0] for i in target_intervals])*100.0/float(target_lenght)

                        if query_coverage >= min_coverage and target_coverage >= min_coverage:
                            distance_matrix[query_index][target_index] = 0
                            distance_matrix[target_index][query_index] = 0

            distance_matrix[query_index][query_index] = 0

        return np.array(distance_matrix), all_queries

    ### Routines for executing local MMseqs2 to find the families

    def get_queries_labels(seq_fasta):

        all_queries = []
        with open(seq_fasta, 'r') as infasta:
            for line in infasta:
                if line.startswith('>'):
                    query = line.split('|')[0].strip('>')
                    all_queries.append(query)
        return all_queries

    def run_mmseqs_for_flanking_sequences(seq_fasta, num_threads = None, max_evalue = None, num_iterations = None, mmseqs = None, tmp_folder = None, sensitivity=7.5, min_coverage=None):

        print(' ... ...Running MMseqs')
        mmseqs_outfile = '{}/{}_{}.mmseqs'.format(tmp_folder, seq_fasta.split('/')[-1][:-6], max_evalue)

        run_mmseqs = sp.Popen(['mmseqs', 'easy-search', seq_fasta, seq_fasta, mmseqs_outfile, tmp_folder, '-e', str(max_evalue), '-s', str(sensitivity), '-c', str(min_coverage), '--num-iterations', str(num_iterations), '--threads', str(num_threads), '--format-output', 'query,target,evalue'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
        stdout, stderr = run_mmseqs.communicate()

        if len(stderr) > 0:
            try:
                run_mmseqs = sp.Popen(['{}mmseqs'.format(mmseqs), 'search', seq_fasta, seq_fasta, mmseqs_outfile, tmp_folder, '-e', str(max_evalue), '-s', str(sensitivity), '-c', str(min_coverage), '--num-iterations', str(num_iterations), '--slice-search', '--threads', str(num_threads), '--format-output', 'query,target,evalue'], stdout=sp.PIPE, stderr=sp.PIPE, stdin=sp.PIPE)
                stdout, stderr = run_mmseqs.communicate()
            except:
                print(stderr)
                print("\n--> ERROR:  There's no MMseqs installation")
                exit()

        return mmseqs_outfile

    def extract_distance_matrix_from_mmseqs_output(mmseqs_results, all_queries, default_base = None):

        distance_matrix = [[default_base if i!=j else 0 for i in all_queries] for j in all_queries]

        all_queries = {query: i for i, query in enumerate(all_queries)}

        with open(mmseqs_results, 'r') as mmseqs_records:
            for hsp in mmseqs_records:
                hsp = hsp.split()
                if len(hsp) > 0:
                    query = hsp[0].split('|')[0]
                    query_index = all_queries[query]

                    target = hsp[1].split('|')[0]
                    if target != query:
                        target_index = all_queries[target]
                        distance_matrix[query_index][target_index] = 0
                        distance_matrix[target_index][query_index] = 0
        
        return np.array(distance_matrix)
    
    out_label = out_label.replace('|', '_')
    out_dir = '{}/{}_all_against_all_searches'.format(os.getcwd(), out_label)
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
 
    flanking_fasta, sequences_len = write_flanking_sequences_to_fasta(in_syntenies, out_dir, out_label, mode = mode)

    if method == 'psiblast':
        sequences_database = make_blast_database_from_fasta(flanking_fasta, blast = blast)
        blast_results = run_blast_for_flanking_sequences(flanking_fasta, sequences_database, num_threads = num_threads, num_alignments = num_alignments, max_evalue = max_evalue, num_iterations = num_iterations, blast = blast, tmp_folder = tmp_folder)
        distance_matrix, queries_labels = extract_distance_matrix_from_blast_output(blast_results, default_base = default_base, sequences_lengths = sequences_len, min_coverage = min_coverage)

    elif method == 'mmseqs':
        queries_labels = get_queries_labels(flanking_fasta)
        mmseqs_results = run_mmseqs_for_flanking_sequences(flanking_fasta, num_threads = num_threads, max_evalue = max_evalue, num_iterations = num_iterations, min_coverage = min_coverage/100, mmseqs = mmseqs, tmp_folder = tmp_folder)
        distance_matrix = extract_distance_matrix_from_mmseqs_output(mmseqs_results, queries_labels, default_base = default_base)

    return distance_matrix, queries_labels

def find_clusters_in_distance_matrix(distance_matrix, t = 0):

    distance_matrix = distance.squareform(distance_matrix)
    linkage = hierarchy.linkage(distance_matrix, method = 'single')
    clusters = hierarchy.fcluster(linkage, t, criterion = 'distance')
    clusters = [int(i) for i in clusters]

    return clusters

def mask_singleton_clusters(clusters_list, mask = 0):

    new_clusters_list = []

    for value in clusters_list:
        if list(clusters_list).count(value) == 1:
            new_clusters_list.append(mask)
        else:
            new_clusters_list.append(value)

    return new_clusters_list

def get_protein_families_summary(in_syntenies, write_to_file = True, out_label = None):

    families = {}

    for target in in_syntenies:
        for i, family in enumerate(in_syntenies[target]['flanking_genes']['families']):
            curr_ncbi_code = in_syntenies[target]['flanking_genes']['ncbi_codes'][i]

            if family not in families:
                families[family] = {'name': [], 'members': [], 'all_names': []}

            if curr_ncbi_code not in families[family]['members']:

                if family != 0:
                    families[family]['all_names'].append(in_syntenies[target]['flanking_genes']['names'][i])
                else:
                    families[family]['all_names'] = ['Non-conserved']

                families[family]['members'].append(curr_ncbi_code)

    for family in families:
        if len(set(families[family]['all_names'])) > 1 and 'hypothetical protein' in set(families[family]['all_names']):
            families[family]['name'] = [name for name in families[family]['all_names'] if name != 'hypothetical protein']
        else:
            families[family]['name'] = families[family]['all_names']

        try:
            families[family]['name'] = statistics.mode(families[family]['name'])
        except:
            families[family]['name'] = families[family]['name'][0]

    if 10000 in families:
        n_pseudogenes = len(families[10000]['members'])
    else:
        n_pseudogenes = 0

    if 0 in families:
        n_nonconserved = len(families[0]['members'])
    else:
        n_nonconserved = 0

    print(' ... Found {} conserved protein families, {} pseudogenes and {} non-conserved protein coding regions'.format(len([i for i in families if i != 0 and i != 10000]), n_pseudogenes, n_nonconserved))

    if write_to_file:
        out_file = '{}_protein_families_summary.txt'.format(out_label)
        with open(out_file, 'w') as outf:
            for family in sorted(list(families.keys())):
                if families[family]['name'] != 'Non-conserved':
                    outf.write('\n ### Family: {} -> {}\n\n'.format(family, families[family]['name']))
                    for i, member in enumerate(families[family]['members']):
                        outf.write('     {}\t{}\n'.format(member, families[family]['all_names'][i]))

    return families

def find_and_add_protein_families(in_syntenies, out_label = 'default', num_threads = 1, num_alignments = 50000, max_evalue = 1e-3, num_iterations = 1, psiblast = 'psiblast', mmseqs = 'mmseqs', min_coverage = 70, default_base = 10, tmp_folder = '/tmp', method = 'psiblast'):
    ordered_ncbi_codes = []

    print(' ... Doing all against all searches with {}'.format(method))

    distance_matrix, ordered_ncbi_codes = compute_all_against_all_distance_matrix(in_syntenies, out_label = out_label, num_threads = num_threads, num_alignments = num_alignments, max_evalue = max_evalue, num_iterations = num_iterations, min_coverage = min_coverage, method = method, mmseqs = mmseqs, blast = psiblast, default_base = default_base, tmp_folder = tmp_folder)
    protein_clusters = find_clusters_in_distance_matrix(distance_matrix)
    protein_clusters = mask_singleton_clusters(protein_clusters)

    # define the correct numbering of the families so that the reference family is the largest, pseudogenes are > 10000 and all the others
    # are continuous
    curr_numbers = []

    for target in in_syntenies:
        in_syntenies[target]['flanking_genes']['families'] = []
        for i, ncbi_code in enumerate(in_syntenies[target]['flanking_genes']['ncbi_codes']):

            protein_name = in_syntenies[target]['flanking_genes']['names'][i]

            protein_family = protein_clusters[ordered_ncbi_codes.index(ncbi_code)]

            if protein_name == 'pseudogene':
                protein_family = 10000
            if ncbi_code == target:
                protein_family = max(protein_clusters)+1

            in_syntenies[target]['flanking_genes']['families'].append(protein_family)
            curr_numbers.append(protein_family)

            if ncbi_code == in_syntenies[target]['assembly_id'][0]:
                in_syntenies[target]['target_family'] = protein_family

    curr_numbers = sorted(list(set(curr_numbers)))
    for target in in_syntenies:
        for i, ncbi_code in enumerate(in_syntenies[target]['flanking_genes']['ncbi_codes']):
            protein_name = in_syntenies[target]['flanking_genes']['names'][i]
            protein_family = in_syntenies[target]['flanking_genes']['families'][i]

            if protein_family <= max(protein_clusters):
                if protein_family != range(min(protein_clusters), max(protein_clusters)+1)[curr_numbers.index(protein_family)]:
                    protein_family = range(min(protein_clusters), max(protein_clusters)+1)[curr_numbers.index(protein_family)]
                    in_syntenies[target]['flanking_genes']['families'][i] = protein_family

        # in_syntenies[target]['target_family'] = in_syntenies[target]['flanking_genes']['families'][in_syntenies[target]['flanking_genes']['ncbi_codes'].index(in_syntenies[target]['assembly_id'][0])]
        in_syntenies[target]['target_family'] = in_syntenies[target]['flanking_genes']['families'][in_syntenies[target]['flanking_genes']['ncbi_codes'].index(in_syntenies[target]['assembly_id'])]

    protein_families = get_protein_families_summary(in_syntenies, write_to_file = False, out_label = out_label)

    os.system('rm -r *_all_against_all_searches')

    return in_syntenies, protein_families

# 4.5 Draw the genomic context and save it as a .png
def draw_genomic_context(all_syntenies, protein_families_summary, targets_dict, save_dir='.'):
    families_flat = {}
    color_dict = {}
    cmap = cm.get_cmap('tab20', len(protein_families_summary.keys()))
    color_list = [plt.colors.rgb2hex(cmap(i)[:3]) for i in range(cmap.N)]

    i = 0
    for key, values in protein_families_summary.items():
        for member in values['members']:
            if member in all_syntenies.keys():
                families_flat[member] = {
                'family': values['name'],
                'color': '#38A9AB'
                }
                color_dict['Target proteins'] = '#38A9AB'
            elif values['name'] == 'Non-conserved':
                families_flat[member] = {
                'family': 'Non-conserved',
                'color': '#E8EDE9'
                }
                color_dict['Non-conserved'] = '#E8EDE9'
            else:
                if values['name'] in color_dict.keys():
                    families_flat[member] = {
                    'family': values['name'],
                    'color': color_dict[values['name']]
                    }
                    pass 
                else:
                    families_flat[member] = {
                    'family': values['name'],
                    'color': color_list[i]
                    }
                    color_dict[values['name']] = color_list[i]
        i += 1

    operon_list = []
    for key, values in all_syntenies.items():
        stringdict = {}
        figure_size = len(range(values['flanking_genes']['relative_starts'][0], values['flanking_genes']['relative_ends'][-1]))
        labels = values['flanking_genes']['ncbi_codes']
        descriptions = values['flanking_genes']['names']

        cds_list = []
        for i in range(len(values['flanking_genes']['relative_starts'])-1):
            if values['flanking_genes']['directions'][i] == '+':
                sense = 1
            else:
                sense = -1
            cds = (values['flanking_genes']['relative_starts'][i], values['flanking_genes']['relative_ends'][i], sense)
            cds_list.append(cds)

        stringdict = {
        'name': key,
        'size': figure_size,
        'cds_list': cds_list,
        'labels': labels,
        'descriptions': descriptions
        }
        operon_list.append(stringdict)

    gv = GenomeViz(tick_style="axis")
    for operon in operon_list:
        uniprot = ''
        name, size, cds_list = operon["name"], operon["size"], operon["cds_list"]
        if name in targets_dict.keys():
            figname = targets_dict[name]['TF']
            uniprot = targets_dict[name]['TF'].split('|')[1]
        relative_start = min(cds_list)[0]
        correction = 0 - relative_start

        if uniprot != '':
            name = '{}\n{}'.format(name, uniprot)
        track = gv.add_feature_track(name, size)
        for idx, cds in enumerate(cds_list, 0):
            start, end, strand = cds
            start = start + correction
            end = end + correction
            label = operon['labels'][idx]
            try:
                track.add_feature(start, end, strand, label=label, linewidth=1, labelrotation=30, labelvpos="top", labelhpos="center", labelha="center", facecolor=families_flat[label]['color'])
            except:
                continue

    fig = gv.plotfig()

    markers = [Line2D([], [], marker=">", color=color, label=family) for family, color in color_dict.items()]
    
    fig.legend(handles=markers, loc='center left', bbox_to_anchor=(1, 0.5), title="Families")

    fig.savefig("{}.png".format(save_dir+figname))


######################################################## PIPELINE CORE ############################################################

# ARGUMENTS
def arguments():
    '''
    Define the input arguments of the program.
    '''

    parser = argparse.ArgumentParser(description = 'ContexTF 1.0: functional annotation of prokaryotic regulatory elements')

    required_arg = parser.add_argument_group('Required arguments')
    optional_arg = parser.add_argument_group('Optional arguments')

    # Required arguments
    required_arg.add_argument('-fna', '--fasta_file', dest = 'fastafile', action = 'store', default = None, required = True, help = 'Required argument. It must be a genomic FASTA formatted file containing the whole genome to study.')
    required_arg.add_argument('-gff3', '--gff3_file', dest = 'gff3file', action = 'store', default = None, required = True, help = 'Required argument. It must be a GFF3 formatted file containing the annotations of the genomic FASTA file.')
    
    # Optional arguments 
    optional_arg.add_argument('--tmalign', dest = 'tmalign', action = 'store_true', default = None, required = False, help = 'Optional argument. If defined, ContexTF will perform structural comparisons between the reference TF and the candidates.')    
    optional_arg.add_argument('--tmalign_min_threshold', dest = 'tmalign_min_threshold', type=float, default=0.5, required = False, help = 'Optional argument. Takes a float between 0 and 1, o.5 as default. If defined, ContexTF will only consider as candidates superimpositions with TM-score above the threshold.')
    optional_arg.add_argument('--n_flanking', dest = 'n_flanking', type=int, default=6, required = False, help = 'Optional argument. Number of flanking genes to plot upstream and downstream of each TF candidate, 4 as default value.')    
    optional_arg.add_argument('--n_flanking3', dest = 'n_flanking3', type=int, default=6, required = False, help = "Optional argument. Number of flanking genes to plot in the 3' site.")
    optional_arg.add_argument('--n_flanking5', dest = 'n_flanking5', type=int, default=6, required = False, help = "Optional argument. Number of flanking genes to plot in the 5' site.")
    optional_arg.add_argument('--GC_method', dest='GC_method', type=str, default='psiblast', choices=['mmseqs', 'psiblast'], required = False, help ='Optional argument. Select MMseqs2 or PSIblast to perform the all-against-all searches to assign protein families and plot the genomic context, default: psiblast.')
    optional_arg.add_argument('--mmseqs_path', dest='mmseqs_path', type=str, default=None, required = False, help ='Optional argument. Provide the path of your mmseqs2 local installation, if not in path.')
    optional_arg.add_argument('--blast_path', dest='blast_path', type=str, default=None, required = False, help = 'Optional argument. Provide the path of your blast local installation, if not in path.')
    
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

# PIPELINE CORE MODULES 
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

    if os.path.isdir('1_HMMER') == False:
        os.mkdir('1_HMMER')    
    os.chdir('1_HMMER')
    current_directory = os.getcwd()+'/'

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
    if os.path.isdir('2_TMalign') == False:
        os.mkdir('2_TMalign')    
    os.chdir('2_TMalign')
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

def GC_module(hmmer_df, path2gff=None, path2fna=None, TM_Align_file=None, n_flanking5=None, n_flanking3=None, executable_path='.', n_module=2, method='psiblast', mmseqs=None, psiblast=None):
    '''
    GC_module can be executed with hmmer or TM-Align output (TM_Align_file)
        1. Downloads the genomic files for every protein 
        2. Retrieves the flanking genes upstream and downstream of every TF (n_flanking3 and n_flanking5)
        3. Uses MMseqs2 or PSIblast to perform an all-against-all sequence search to cluster proteins of the same family
        4. Generates a genomic context figure for each reference TF and its candidates, with the proteins colored by the family clustering
    '''

    if os.path.isdir('{}_GC'.format(n_module)) == False:
        os.mkdir('{}_GC'.format(n_module))    
    os.chdir('{}_GC'.format(n_module))
    pwd = os.getcwd()
    # 1. Get the target dict with targets and gff and fna path for all targets
    targets_dict = retrieve_targets_dict(hmmer_df, path2gff=path2gff, path2fna=path2fna, TM_Align_file=TM_Align_file, executable_path=executable_path)

    plasmids = ['AAP70493.1','ADO85569.1','AEM66515.1','AAA21920.1','BAD11074.1','AAP83141.1','BAB32408.1','AFI98560.1','ACO06750.1','AAB36583.1','AAK38101.1', 'ANB66399.1', 'AAG00065.1', 'AAA25445.1', 'AAA68939.2', 'CAC87048.1', 'BAA36282.1']
    
    if os.path.isdir('GC_figures') == False:
        os.mkdir('GC_figures')
    figures_dir = pwd+'/GC_figures/'

    # Routine to resume from where it left off 
    already_figures = [str(x).strip('.png') for x in os.listdir(pwd+'/GC_figures/') if x.endswith('.png')]
    if len(already_figures) > 0: 
        logging.info('Found {}/{} GC figures. Plotting the {} remaining GC figures.'.format(len(already_figures), len(targets_dict.keys()), len(targets_dict.keys())-len(already_figures)))
    else: 
        logging.info('Plotting the genomic context for {} reference TF and their candidates.'.format(len(targets_dict.keys())))
    
    # 2. Retrieve the genomic context for each target and candidates
    for key in targets_dict.keys():
        if key not in plasmids and targets_dict[key]['TF'] not in already_figures:
            logging.info('GC module: plotting genome context for {}'.format(targets_dict[key]['TF']))
            all_syntenies = {}
            target_dict = {}
            target_dict[key] = targets_dict[key]

            all_syntenies = genomic_context(target_dict, n_flanking5=n_flanking5, n_flanking3=n_flanking3)

            # 3. Execute MMseqs2 o BLASTp to find families
            all_syntenies, protein_families_summary = find_and_add_protein_families(all_syntenies, out_label=targets_dict[key]['TF'], method=method, mmseqs=mmseqs, psiblast=psiblast)

            # 4. Draw the genomic context and save a figure
            draw_genomic_context(all_syntenies, protein_families_summary, targets_dict, save_dir=figures_dir)

            # Clean files 
            os.system('rm *.fa')
        elif key in plasmids:
            logging.info('No genomic context for {}: {} only sequence reference in NCBI databases is a plasmid/biosynthetic gene cluster ({}) it is not annotated in any species.'.format(targets_dict[key]['TF'], key, targets_dict[key]['fna'].split('/')[-1].strip('_genomic.fna ')))
    
# PIPELINE EXECUTION
def main():
    '''
    Execute the ContexTF pipeline:
        1) Sequence homology: execute HMMER module.
        2) Structural homology: if flagged, executes the TMAlign module.
        3) Genomic context homology: executes the GC module with the output of HMMER or TMAlign.
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

    GC_method = args.GC_method
    mmseqs_path = args.mmseqs_path
    blast_path = args.blast_path

    # Execute modules
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
                GC_module(hmmer_df, path2gff=gff3file, path2fna=fastafile, TM_Align_file=TMalign_results_file, n_flanking5=n_flanking5, n_flanking3=n_flanking3, executable_path=executable_path, n_module=3, method=GC_method, mmseqs=mmseqs_path, psiblast=blast_path)
            else: 
                # 2. HMMER
                hmmer_df = HMMER_module(organism_proteome, current_directory=current_directory)
                # 3. Genomic Context plotting of the HMMER candidates
                GC_module(hmmer_df, path2gff=gff3file, path2fna=fastafile, n_flanking5=n_flanking5, n_flanking3=n_flanking3, executable_path=executable_path, method=GC_method, mmseqs=mmseqs_path, psiblast=blast_path)
        else:
            logging.error('Please provide a GFF3 formatted file.')
            raise IncorrectInput(args.gff3file)
    else:
        logging.error('Please provide a genomic FASTA formatted file.')
        raise IncorrectInput(args.fastafile)
    
    logging.info('ContexTF execution finished :)')

if __name__ == '__main__':
    main()
