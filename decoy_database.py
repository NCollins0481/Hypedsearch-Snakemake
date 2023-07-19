import sys
import os
import database
import random

# database_dir = sys.argv[1]
# prot_dir = sys.argv[2]
database_dir = '/home/naco3124/snakemake/database/Comet_filtered_db.fasta'
prot_dir = '/home/naco3124/snakemake/database'
print(prot_dir)

def invert_database(protein_list):
    inverted_proteins = []
    for protein in protein_list:
        headers = protein[0].split("|")
        decoy_header = ">" + headers[0] + "|" + headers[1] + "|DECOY_" + headers[2]
        
        prot_seq = protein[1]
        l = list(prot_seq)
        random.shuffle(l)
        shuffled = ''.join(l)

        inverted_prot = (decoy_header, shuffled)
        inverted_proteins.append(inverted_prot)
    
    return inverted_proteins

def write_inverted(inverted, prot_dir):
    with open(os.path.join(prot_dir, "Decoy_Comet_filtered_db.fasta"), 'w') as d:
        
        for tup in inverted:
            d.write(tup[0] + "\n")
            prot_seq = tup[1]
            prot_seq = '\n'.join(prot_seq[i:i+70] for i in range(0, len(prot_seq), 70))
            d.write(prot_seq + "\n\n")
            
proteins = database.build(database_dir)
inverted = invert_database(proteins.proteins)
write_inverted(inverted, prot_dir)