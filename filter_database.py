import database
import os
import argparse
from hypedsearch.src.objects import Spectrum
from hypedsearch.src.preprocessing import spectra_filtering
from pyteomics import mzml

# comet_results = sys.argv[1]
# prot_path = sys.argv[2]
# prot_dir = sys.argv[3]
parser = argparse.ArgumentParser(description='Tool for making a database from Comet results')
parser.add_argument('--Comet_results', dest='comet_results', nargs='*', default=[''], type = str)
parser.add_argument('--prot_path', dest='prot_path', type=str)
parser.add_argument('--prot_dir', dest='prot_dir', type=str)
parser.add_argument('--auto-includes', dest='auto_includes', nargs='*', default=[''], type=str)
parser.add_argument('--digest-left', dest='digest_left', nargs='*', default=[''], type = str, help='The Amino Acid for which the digest cuts left of. Default=None')
parser.add_argument('--digest-right', dest='digest_right', nargs='*', default=[''], type = str, help='The Amino Acid for which the digest cuts right of. Default=None')
parser.add_argument('--peak-filter', dest='peak_filter', type=int, default=0, help='The number of peaks to take from a spectrum. The most abundant peaks will be taken. Leave blank if you want no filter or to use relative abundance filter. Defualt=0')
parser.add_argument('--abundance-filter', dest='rel_abund_filter', type=float, default=0.0, help='Take only peaks from a spectrum where the abundance of the peak is >= the percentage give. Leave blank if you want no filter or to use peak filter. Default=0.0')
parser.add_argument('--spectra-folder', dest='spectra_folder', type=str, default="", help='Directory to the folder containing all the spectra files')
parser.add_argument('--ppm-tolerance', dest='ppm_tol', type=float, default=0.0, help='Tolerance in ppm to identify peaks')


args = parser.parse_args()
# print(args)

comet_results = args.comet_results
prot_path = args.prot_path
prot_dir = args.prot_dir
auto_includes = args.auto_includes
left_digest = args.digest_left
right_digest = args.digest_right
spectra_folder = args.spectra_folder
ppm_tolerance = args.ppm_tol
peak_filter = args.peak_filter
relative_abundance_filter = args.rel_abund_filter

# comet_results = ["/home/naco3124/snakemake/output/comet_run_1/BMEM_AspN_Fxn5.txt"]
# prot_path = "/home/naco3124/snakemake/database/UniProt_mouse.fasta"
# prot_dir = "/home/naco3124/snakemake/database"
# auto_includes = ["DLQTLAL-EVE", "DLQTLAL-NAAR", "DLQTLAL-WSRM", "DPGVAQLELGG-EVEDPQVAQLELGGGPGAG"]
# left_digest = ["D"]
# right_digest = [""]
# spectra_folder = "/home/naco3124/snakemake/spectra"
# ppm_tolerance = 10
# peak_filter = 25
# relative_abundance_filter = .01

# print(comet_results, prot_path, prot_dir)

def file_exists(file_name: str) -> bool:
    '''Determine if a file exists

    :param file_name: Path to the file in question
    :type file_name: str
    
    :returns: True if the file exists
    :rtype: bool
    '''
    
    return os.path.isfile(file_name)

def read(filename: str, peak_filter=0, relative_abundance_filter=0) -> list:
    if not file_exists(filename):
        print('File {} not found. Please make sure that this file exists'.format(filename))
        return

    spectra = []
    
    filecontents = mzml.read(filename)

    content: dict

    for spec_num, content in enumerate(filecontents):

        masses = list(content['m/z array'])
        abundances = list(content['intensity array'])

        # peak filter if the number is > 0
        if peak_filter > 0:
            masses, abundances = spectra_filtering.peak_filtering(masses, abundances, peak_filter)

        # if peak filter is not set and relative abundances is, filter by that
        elif relative_abundance_filter > 0:

            # if an integer is given, make it a float in the range (0, 1)
            while relative_abundance_filter > 1:
                relative_abundance_filter /= 100

            masses, abundances = spectra_filtering.relative_abundance_filtering(masses, abundances, relative_abundance_filter)

        # get the total intensity
        ti = sum(abundances)

        # get the precursor and its charge
        # we will assume its the first entry in the list
        precursor = None
        precursor_charge = 0

        if not len(content['precursorList']['precursor']) or not len(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon']):
            precursor = max(masses)
            precursor_charge = 1

        else:
            precursor = float(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['selected ion m/z'])
            precursor_abundance = float(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['peak intensity'])
            precursor_charge = int(content['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]['charge state'])

        # get the id
        id = content.get('id', '')

        retention_time = content['scanList']['scan'][0]['scan start time']
        
        spectra.append(
            Spectrum(
            spec_num,
            masses,
            abundances,
            precursor,
            precursor_charge,
            filename, 
            id,
            retention_time,
            precursor_abundance
        ))

    return spectra


def load(filename: str, peak_filter: int = 0, relative_abundance_filter: float = 0.0) -> list:
    if not file_exists(filename):
        print(f'File {filename} not found. Please make sure that this file exists')

    # return based on file type
    ext = filename.split('.')[-1]

    # if ext.lower() == 'mzxml':
    #     foo = mzXML.read(filename, peak_filter, relative_abundance_filter)
    #     return foo

    if ext.lower() == 'mzml':
        return read(filename, peak_filter, relative_abundance_filter)

    else:
        print(f'File {filename} is not of supported types (mzML, mzXML)')
        return []

def load_spectra(
    spectra_file, ppm_tol: int, peak_filter: int = 0, relative_abundance_filter: float = 0.0):
    linear_spectra = []
    all_spectra = []
    these_spectra = load(
        spectra_file, 
        peak_filter=peak_filter, 
        relative_abundance_filter=relative_abundance_filter
    )
    all_spectra += these_spectra
    # leave next 2 lines commented; uncomment only to test just specific indices
    # index_list = [5049] #For using a condensed database
    # these_spectra, all_spectra = reduce_database(all_spectra, these_spectra, index_list)
    linear_spectra += list(set([
        x for spectrum in these_spectra for x in spectrum.mz_values
    ]))
    linear_spectra.sort()
    return (all_spectra)


def get_spectra_files(spectra_folder):
    spectra_files = []
    for (root, _, filenames) in os.walk(spectra_folder):
        for fname in filenames:
            spectra_files.append(os.path.join(root, fname))
    return spectra_files

proteins = database.build(prot_path)
spectra_files = get_spectra_files(spectra_folder)

def assign_abundances(prot_unique, spectrum_mapping):
    #output a dictionary mapping proteins to their abundances
    protein_abundances = dict()
    for label, spec_id, prec_abundance in spectrum_mapping: # (filename, spec_id, precursor_abundance) -> all proteins that share this spectrum
        #If separate proteins share a peptide sequence, the count goes to the protein with the highest number of unique spectra
        current_max = 0
        for prot in spectrum_mapping[(label, spec_id, prec_abundance)]:
            if len(prot_unique[prot]) > current_max:
                current_max = len(prot_unique[prot])
                most_unique = prot
        if most_unique not in protein_abundances.keys():
            protein_abundances[most_unique] = 0
        protein_abundances[most_unique] += prec_abundance 
    
    return protein_abundances

def get_parents(comet_results, ppm_tol, peak_filter, relative_abundance_filter):
    # As input we have each Comet file and then as output we have the top X most abundant proteins
    
    spectrum_mapping = dict()#Dictionary holding all data for proteins
    prot_unique = dict() #Dictionary mapping proteins to the number of unique spectra matched to that protein
    prot_abundances = dict() #Dictionary mapping proteins to their abundance
    
    for file in comet_results: #Only take the top hit
        prev_num = -1
        with open(file, "r") as f:
            for i, line in enumerate(f):
                A = line.split("\t")
                if i == 0:
                    filepath = A[1]
                    input_spectra = load_spectra(filepath+".mzML", ppm_tol, peak_filter, relative_abundance_filter)
                    continue
                elif i == 1:
                    continue
                else:
                    spec_id = int(A[0])
                    if spec_id != prev_num: #take the top hit
                        prev_num = spec_id
                        comet_seq = A[11]
                        comet_score = A[9]
                        if int(comet_score) > len(comet_seq)/2: #Only add things where at least 1/8 of the ions we'd expect to see matched
                            target_spectrum = input_spectra[spec_id]
                            precursor_abundance = target_spectrum.precursor_abundance
                            id = (filepath, spec_id, precursor_abundance)
                            parent_seqs = A[15].split(",")
                            for parent in parent_seqs:
                                if id not in spectrum_mapping.keys():
                                    spectrum_mapping[id] = []
                                spectrum_mapping[id].append(parent)
                                
                                if parent not in prot_unique.keys():
                                    prot_unique[parent] = []
                                prot_unique[parent].append(spec_id)
    
    prot_abundances = assign_abundances(prot_unique, spectrum_mapping)
    
    parents = sorted(prot_abundances, key = lambda x: prot_abundances[x], reverse=True)
    return parents
    

# def get_parents(file_list):
#     parents = dict()
#     # for file in result_list:
#     for file in file_list:
#         with open(file, "r") as f:
#             next(f)
#             next(f)
#             for line in f:
#                 A = line.split("\t")
#                 comet_seq = A[11]
#                 comet_score = A[9]
#                 if int(comet_score) > len(comet_seq)/2: #Only add things where at least 1/8 of the ions we'd expect to see matched
#                     parent_seqs = A[15].split(",")
#                     for seq in parent_seqs:
#                         prot_name = seq.split("|")[2]
#                         if prot_name not in parents.keys():
#                             parents[prot_name] = 0
#                         if "DECOY" not in prot_name:
#                             parents[prot_name] += 1
#     sorted_parents = sorted(parents, key=lambda x: parents[x], reverse=True)
#     for i in range(0, len(sorted_parents)):
#         if ("INS1_MOUSE" in sorted_parents[:i]) and ("INS2_MOUSE" in sorted_parents[:i]) and ("CMGA_MOUSE" in sorted_parents[:i]):
#             print("The 3 good proteins in sorted parents up to:", i)
#             return sorted_parents[:i]

parents = get_parents(comet_results, ppm_tolerance, peak_filter, relative_abundance_filter)
parents = parents[:50]

def build_small_db(parent_list, proteins, auto_includes, digest_left, digest_right):
    new_proteins = []
    for parent in parent_list:
        found = False
        if "Hybrid" not in parent and "DECOY" not in parent:
            for protein in proteins:
                description = protein[0]
                if parent in description:
                    found = True
                    new_proteins.append(protein)
            if found == False:
                print("investigate")
        else:
            print(parent, "not found in proteins")
            
    #now adding in the custom proteins
    for seq in auto_includes:
        left_seq = seq.split("-")[0]
        seq_found = False
        for protein in new_proteins:
            protein_seq =  protein[1]
            if left_seq in protein_seq:
                seq_found = True
                break
        
        if seq_found == False:
            for protein in proteins:
                prot_seq = protein[1]
                if left_seq[0] in digest_left:
                    if left_seq in prot_seq:
                        new_proteins.append(protein)
                        break
                        
    for seq in auto_includes:
        break_flag = False
        right_seq = seq.split("-")[1]
        seq_found = False
        for protein in new_proteins:
            protein_seq =  protein[1]
            for char in digest_left:
                if right_seq + char in protein_seq:
                    seq_found = True
                    break_flag = True
                    break
                
                if break_flag:
                    break_flag = False
                    break
        
        if seq_found == False:
            for protein in proteins:
                prot_seq = protein[1]
                if right_seq[-1] in digest_right:
                    if right_seq in prot_seq:
                        new_proteins.append(protein)
                        break
                else:
                    for char in digest_left:
                        if right_seq + char in prot_seq:
                            new_proteins.append(protein)
                            break_flag = True
                            break
                    
                    if break_flag:
                        break_flag = False
                        break
        
    new_proteins = set(new_proteins)
    return new_proteins

all_proteins = build_small_db(parents, proteins.proteins, auto_includes, left_digest, right_digest)

with open(os.path.join(prot_dir, "Comet_filtered_db.fasta"), 'w') as d:
    for protein in all_proteins:
        d.write('>' + protein[0] + '\n')
        prot_seq = protein[1]
        prot_seq = '\n'.join(prot_seq[i:i+70] for i in range(0, len(prot_seq), 70))
        d.write(prot_seq + "\n\n")
