import argparse
import glob
import database
import os

CLI=argparse.ArgumentParser()
CLI.add_argument("--initial-db", type=str, default='')
CLI.add_argument("--HS-outputs", nargs="*", type=str, default=[])

args = CLI.parse_args()

def get_hybrids(output_dir):
    hybrids = []
    for file in output_dir:
        with open(file, "r") as g:
            #next(g)
            for line in g:
                A = line.split("\t")
                hybrid_label = A[1]
                if hybrid_label == "Hybrid":
                    hybrid_seq = A[2]
                    left_parent = A[5]
                    right_parent = A[6]
                    protein_name = "sp|P01325|Hypedsearch_Hybrid_" + hybrid_seq + " " + left_parent + "-" + right_parent + " " + hybrid_seq
                    full_protein = (protein_name, hybrid_seq.replace("-", ""))
                    hybrids.append(full_protein)
    return hybrids

def append_to_existing_db(initial_db_filepath, hybrid_proteins):
    UniProt = database.build(initial_db_filepath)
    with open(os.path.join(os.path.dirname(initial_db_filepath), "HS_hybrid_database.fasta"), "w") as d:
        for protein in UniProt.proteins:
            d.write('>' + protein[0] + '\n')
            prot_seq = protein[1]
            prot_seq = '\n'.join(prot_seq[i:i+70] for i in range(0, len(prot_seq), 70))
            d.write(prot_seq + "\n\n")
        for protein in hybrid_proteins:
            d.write('>' + protein[0] + '\n')
            prot_seq = protein[1]
            prot_seq = '\n'.join(prot_seq[i:i+70] for i in range(0, len(prot_seq), 70))
            d.write(prot_seq + "\n\n")
            
HS_hybrids = get_hybrids(args.HS_outputs)
append_to_existing_db(args.initial_db, HS_hybrids)