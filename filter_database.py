import database
import glob
import sys
import os
import argparse

# comet_results = sys.argv[1]
# prot_path = sys.argv[2]
# prot_dir = sys.argv[3]
parser = argparse.ArgumentParser(description='Tool for making a database from Comet results')
parser.add_argument('--Comet_results', dest='comet_results', nargs='*', default=[''], type = str)
parser.add_argument('--prot_path', dest='prot_path', type=str)
parser.add_argument('--prot_dir', dest='prot_dir', type=str)
args = parser.parse_args()
# print(args)

comet_results = args.comet_results
prot_path = args.prot_path
prot_dir = args.prot_dir

# print(comet_results, prot_path, prot_dir)

proteins = database.build(prot_path)

# def build_table(directory):
#     for filepath in directory+"/*.txt":
#         sequences, ids = [], []
#         with open(filepath, 'r') as t:
#             for line in t:
#                 A = line.split('\t')
#                 seq = A[11]
#                 id = A[0].strip('\"')
#                 sequences.append(seq)
#                 ids.append(id)
#     return sequences, ids

def get_parents(file_list):
    parents = set()
    # for file in result_list:
    for file in file_list:
        with open(file, "r") as f:
            next(f)
            next(f)
            prev_spec_num = -1
            for line in f:
                A = line.split("\t")
                spec_num = A[0]
                if spec_num != prev_spec_num:
                    parent_seqs = A[15].split(",")
                    parent_seq = parent_seqs[0]
                    if "DECOY" not in parent_seq:
                        parent = parent_seq.split("|")[2]
                        parents.add(parent)
                    prev_spec_num = spec_num
    return parents

parents = get_parents(comet_results)

def build_small_db(parent_list, proteins):
    new_proteins = []
    for parent in parent_list:
        found = False
        if "Hybrid" not in parent:
            for protein in proteins:
                description = protein[0]
                if parent in description:
                    found = True
                    new_proteins.append(protein)
            if found == False:
                print("investigate")
        else:
            print(parent, "not found in proteins")
        
    new_proteins = set(new_proteins)
    return new_proteins

all_proteins = build_small_db(parents, proteins.proteins)

with open(os.path.join(prot_dir, "Comet_filtered_db.fasta"), 'w') as d:
    for protein in all_proteins:
        d.write('>' + protein[0] + '\n')
        prot_seq = protein[1]
        prot_seq = '\n'.join(prot_seq[i:i+70] for i in range(0, len(prot_seq), 70))
        d.write(prot_seq + "\n\n")
