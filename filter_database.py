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
parser.add_argument('--auto-includes', dest='auto_includes', nargs='*', default=[''], type=str)
parser.add_argument('--digest-left', dest='digest_left', nargs='*', default=[''], type = str, help='The Amino Acid for which the digest cuts left of. Default=None')
parser.add_argument('--digest-right', dest='digest_right', nargs='*', default=[''], type = str, help='The Amino Acid for which the digest cuts right of. Default=None')

args = parser.parse_args()
# print(args)

comet_results = args.comet_results
prot_path = args.prot_path
prot_dir = args.prot_dir
auto_includes = args.auto_includes
left_digest = args.digest_left
right_digest = args.digest_right

# comet_results = ["/home/naco3124/snakemake/output/comet_run_1/BMEM_AspN_Fxn5.txt"]
# prot_path = "/home/naco3124/snakemake/database/UniProt_mouse.fasta"
# prot_dir = "/home/naco3124/snakemake/database"
# auto_includes = ["DLQTLAL-EVE", "DLQTLAL-NAAR", "DLQTLAL-WSRM", "DPGVAQLELGG-EVEDPQVAQLELGGGPGAG"]
# left_digest = ["D"]
# right_digest = [""]

# print(comet_results, prot_path, prot_dir)

proteins = database.build(prot_path)

def get_parents(file_list):
    parents = dict()
    # for file in result_list:
    for file in file_list:
        with open(file, "r") as f:
            next(f)
            next(f)
            for line in f:
                A = line.split("\t")
                comet_seq = A[11]
                comet_score = A[9]
                if int(comet_score) > len(comet_seq)/2: #Only add things where at least 1/8 of the ions we'd expect to see matched
                    parent_seqs = A[15].split(",")
                    for seq in parent_seqs:
                        prot_name = seq.split("|")[2]
                        if prot_name not in parents.keys():
                            parents[prot_name] = 0
                        if "DECOY" not in prot_name:
                            parents[prot_name] += 1
    sorted_parents = sorted(parents, key=lambda x: parents[x], reverse=True)
    for i in range(0, len(sorted_parents)):
        if ("INS1_MOUSE" in sorted_parents[:i]) and ("INS2_MOUSE" in sorted_parents[:i]) and ("CMGA_MOUSE" in sorted_parents[:i]):
            print("The 3 good proteins in sorted parents up to:", i)
            return sorted_parents[:i]

parents = get_parents(comet_results)

def build_small_db(parent_list, proteins, auto_includes, digest_left, digest_right):
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
