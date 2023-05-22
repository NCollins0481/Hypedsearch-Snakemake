# Snakemake pipeline for Hypedsearch

Automates a run of Hypedsearch and compares the output to Comet.
Made for Linux machines. Windows and MacOS compatability in progress.

# Installation
1. Install [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
2. Install snakemake 
```bash
$> conda install -n base -c conda-forge mamba
```
3. Copy this repository
```bash
$> git clone https://github.com/NCollins0481/Hypedsearch-Snakemake
```
4. Inside *conf*, edit the *config.yaml* and *comet_params* files with filepaths to the target directories and specify the resources allocated
5. If running the small example, go to *hypedsearch/src/preprocessing/preprocessing_utils* and uncomment the two lines to reduce to a few spectra
6. Inside the *Hypedsearch-Snakemake*, run the pipeline with:
```bash
$> conda activate snakemake
$> snakemake -j1 --use-conda
```
