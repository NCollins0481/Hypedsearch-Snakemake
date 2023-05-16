import os

spectra_files = []
for (root, _, filenames) in os.walk('/home/naco3124/snakemake/spectra'):
    for fname in filenames:
        spectra_files.append(os.path.join(root, fname))

