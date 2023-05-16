from types import SimpleNamespace
import os
configfile: 'conf/config.yaml'

config = SimpleNamespace(**config)

spectra_files = []
for (root, _, filenames) in os.walk(config.spectra_dir):
    for fname in filenames:
        spectra_files.append(os.path.join(root, fname))

spectra_files = [spectra_files[0]]

database_file = config.database_file
output_dir = config.output_dir
bin_directory = config.bin_direc

base_files = []
for file in spectra_files:
    e = os.path.splitext(file)[0]
    b = os.path.basename(e)
    base_files.append(b)

rule All:
    input:
        expand(f'{output_dir}/comet_run_1/{{dataset}}.txt', dataset=base_files)

rule GetComet:
    output:
        comet = f'{bin_directory}/comet'
    shell:
        """
        wget -O {output.comet} https://github.com/UWPR/Comet/releases/download/v2023.01.2/comet.linux.exe
        chmod +x {output.comet}
        """

rule RunComet:
    input:
        comet = rules.GetComet.output.comet,
        mzML_files = expand("{spectra_file}", spectra_file=spectra_files),
        database_file = database_file,
        param_file = 'conf/comet.params'
    output:
        output_texts = expand(f'{output_dir}/comet_run_1/{{dataset}}.txt', dataset=base_files),
        output_xml = expand(f'{output_dir}/comet_run_1/{{dataset}}.pep.xml', dataset=base_files)
    shell:
        f"""
        {{input.comet}} -P{{input.param_file}} -D{{input.database_file}} {{input.mzML_files}}
        mkdir -p {output_dir}/comet_run_1
        bash mover.sh {config.spectra_dir} {output_dir}/comet_run_1
        """
