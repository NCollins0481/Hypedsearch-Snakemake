import database
import sys

comet_results = sys.argv[1]
prot_path = sys.argv[2]

proteins = database.build(prot_path)

def build_table(directory):
    [print(x) for x in directory+"/*.txt"]
    for filepath in directory+"/*.txt":
        sequences, ids = [], []
        with open(filepath, 'r') as t:
            for line in t:
                A = line.split('\t')
                seq = A[21]
                id = A[14].strip('\"')
                sequences.append(seq)
                ids.append(id)
    return sequences, ids
    
specmill_sequences, seq_ids = build_table(comet_results)
with open("filter_database_testing.txt", 'w') as f:
    [f.write(x + "\t" + y + "\n") for x in specmill_sequences for y in seq_ids]

def get_proteins(protein_list, seq_ids, sequences):
    proteins = []
    for i, id in enumerate(seq_ids):
        print("On i=", i)
        found = False
        for protein in protein_list:
            description = protein[0]
            if id in description:
                if sequences[i] in protein[1]:
                    proteins.append(protein)
                    found = True
                    continue
        if 'HYBRID' in id:
            print(id) #these need to be added manually
            continue
        elif not found:
            print("Not found at", i)

    print(len(seq_ids), len(proteins))
    new_proteins = set(proteins)
    return new_proteins
            
            

all_proteins = get_proteins(proteins.proteins, seq_ids, specmill_sequences)

with open(prot_path, 'w') as d:
    for protein in all_proteins:
        d.write('>' + protein[0] + '\n')
        prot_seq = protein[1]
        seq_len = len(prot_seq)
        curr_write = 0
        prot_seq = '\n'.join(prot_seq[i:i+70] for i in range(0, len(prot_seq), 70))
        d.write(prot_seq + "\n\n")
