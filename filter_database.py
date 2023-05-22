import database
import glob
import sys
import os

comet_results = sys.argv[1]
prot_path = sys.argv[2]
prot_dir = sys.argv[3]

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

def get_parents(file):
    parents = []
    # for file in result_list:
    with open(file, "r") as f: #This will break when we have to run on full sample!
        next(f)
        for line in f:
            A = line.split("\t")
            parent = A[15]
            parents.append(parent)
    return parents
    
parents = get_parents(comet_results)[1:]

def build_small_db(parent_list, proteins):
    new_proteins = []
    for parent in parent_list:
        found = False
        line_arr = parent.split(",")
        # print("line_arr:", line_arr)
        first_parent = line_arr[0]
        # print("first_parent:", first_parent)
        split_parent = first_parent.split("|")
        # print("split_parent:", split_parent)
        tag = split_parent[0]
        # print("tag:", tag)
        if "DECOY" not in tag:
            name = split_parent[2]
            if "Hybrid" not in name:
                for protein in proteins:
                    description = protein[0]
                    if name in description:
                        found = True
                        new_proteins.append(protein)
            else:
                sections = name.split("--")
                target_section = sections[-1]
                left_parent, right_parent = target_section.split("+")
                for protein in proteins: #TODO: Make this a dictionary
                    description = protein[0]
                    if left_parent in description or right_parent in description:
                        found = True
                        new_proteins.append(protein)
        
        if not found:
            if "DECOY" not in tag:
                print("Manually add:", parent)
        
    new_proteins = set(new_proteins)
    return new_proteins

all_proteins = build_small_db(parents, proteins.proteins)

with open(os.path.join(prot_dir, "Comet_filtered_db.fasta"), 'w') as d:
    for protein in all_proteins:
        d.write('>' + protein[0] + '\n')
        prot_seq = protein[1]
        prot_seq = '\n'.join(prot_seq[i:i+70] for i in range(0, len(prot_seq), 70))
        d.write(prot_seq + "\n\n")
