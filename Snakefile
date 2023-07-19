from types import SimpleNamespace
import os
configfile: 'conf/config.yaml'

config = SimpleNamespace(**config)
digest_left, digest_right = config.digest_left, config.digest_right

spectra_files = []
for (root, _, filenames) in os.walk(config.spectra_dir):
    for fname in filenames:
        spectra_files.append(os.path.join(root, fname))

database_file = config.database_file
database_dir = os.path.dirname(database_file)
print(database_dir)
output_dir = config.output_dir
bin_directory = config.bin_direc
environment_directory = os.path.join(os.path.dirname(database_dir), "environments")

base_files = []
for file in spectra_files:
    e = os.path.splitext(file)[0]
    b = os.path.basename(e)
    base_files.append(b)

rule All:
    input:
        output_texts = expand(f'{output_dir}/comet_run_2/{{dataset}}.txt', dataset=base_files),
        output_xml = expand(f'{output_dir}/comet_run_2/{{dataset}}.pep.xml', dataset=base_files),
        decoy_database = expand(f'{output_dir}/Hypedsearch_outputs_Decoy/HS_{{dataset}}.txt', dataset=base_files)

rule GetComet:
    output:
        comet = f'{bin_directory}/comet.exe'
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

rule CondenseDatabase:
    input:
        # dependencies = rules.GetHypedsearchDependencies.output.dependencies,
        output_texts = rules.RunComet.output.output_texts
    output:
        filtered_db = f'{database_dir}/Comet_filtered_db.fasta'
    conda:
        f'{environment_directory}/Hypedsearch.yaml'
    shell:
        f"""
        python3 -m filter_database --Comet_results {{input.output_texts}} --prot_path {{config.database_file}} --prot_dir {database_dir}
        """

        # echo {{input.output_texts}}
        # echo {{config.database_file}}
        # echo {database_dir}


rule CreateInvertedDB:
    input:
        filtered_database = f'{database_dir}/Comet_filtered_db.fasta'
    output:
        inv_filtered_db = f'{database_dir}/Decoy_Comet_filtered_db.fasta'
    conda:
        f'{environment_directory}/Hypedsearch.yaml'
    shell:
        f"""
        python3 -m decoy_database {{input.filtered_database}} {database_dir}
        """

rule RunHypedsearch:
    input:
        filtered_database = f'{database_dir}/Comet_filtered_db.fasta'
    output:
        Hypedsearch_outputs = expand(f'{output_dir}/Hypedsearch_outputs/HS_{{dataset}}.txt', dataset=base_files)
    conda:
        f'{environment_directory}/Hypedsearch.yaml'
    shell:
        f"""
        cd hypedsearch/src
        python3 -m main --spectra-folder {{config.spectra_dir}}\
        --database-file {{input.filtered_database}}\
        --output-dir {{config.output_dir}}/Hypedsearch_outputs\
        --no-config\
        --min-peptide-len {{config.min_peptide_len}}\
        --max-peptide-len {{config.max_peptide_len}}\
        --tolerance {{config.ppm_tolerance}}\
        --precursor-tolerance {{config.precursor_tolerance}}\
        --peak-filter {{config.num_peaks}}\
        --verbose {{config.verbose}}\
        --cores {{config.cores}}\
        --new-database {{config.new_db}}\
        --digest-left {{config.digest_left}}\
        --digest-right {{config.digest_right}}\
        --n {{config.top_results}}\
        """

rule RunHypedsearchDecoy:
    input:
        filtered_database = f'{database_dir}/Decoy_Comet_filtered_db.fasta'
    output:
        Hypedsearch_outputs = expand(f'{output_dir}/Hypedsearch_outputs_Decoy/HS_{{dataset}}.txt', dataset=base_files)
    conda:
        f'{environment_directory}/Hypedsearch.yaml'
    shell:
        f"""
        cd hypedsearch/src
        python3 -m main --spectra-folder {{config.spectra_dir}}\
        --database-file {{input.filtered_database}}\
        --output-dir {{config.output_dir}}/Hypedsearch_outputs_Decoy\
        --no-config\
        --min-peptide-len {{config.min_peptide_len}}\
        --max-peptide-len {{config.max_peptide_len}}\
        --tolerance {{config.ppm_tolerance}}\
        --precursor-tolerance {{config.precursor_tolerance}}\
        --peak-filter {{config.num_peaks}}\
        --verbose {{config.verbose}}\
        --cores {{config.cores}}\
        --new-database {{config.new_db}}\
        --digest-left {{config.digest_left}}\
        --digest-right {{config.digest_right}}\
        --n {{config.top_results}}\
        """
        # --digest {{config.digest}}\

rule BuildHybridDatabase:
    input:
        # dependencies = rules.GetHypedsearchDependencies.output.dependencies,
        Hypedsearch_outputs = expand(f'{output_dir}/Hypedsearch_outputs/HS_{{dataset}}.txt', dataset=base_files)
    output:
        hybrid_database = f'{database_dir}/HS_hybrid_database.fasta'
    conda:
        f'{environment_directory}/Hypedsearch.yaml'
    shell:
        f"""
        python3 -m build_hybrid_database --initial-db {database_file} --HS-outputs {{input.Hypedsearch_outputs}}
        """

rule RunCometAgain:
    input:
        comet = rules.GetComet.output.comet,
        mzML_files = expand("{spectra_file}", spectra_file=spectra_files),
        param_file = 'conf/comet.params',
        database_file = f'{database_dir}/HS_hybrid_database.fasta'
    output:
        output_texts = expand(f'{output_dir}/comet_run_2/{{dataset}}.txt', dataset=base_files),
        output_xml = expand(f'{output_dir}/comet_run_2/{{dataset}}.pep.xml', dataset=base_files)
    shell:
        f"""
        {{input.comet}} -P{{input.param_file}} -D{{input.database_file}} {{input.mzML_files}}
        mkdir -p {output_dir}/comet_run_2
        bash mover.sh {config.spectra_dir} {output_dir}/comet_run_2
        """

rule CompareHypedsearch:
    input:
        comet1_texts = expand(f'{output_dir}/comet_run_1/{{dataset}}.txt', dataset=base_files),
        comet2_texts = expand(f'{output_dir}/comet_run_2/{{dataset}}.txt', dataset=base_files),
        Hypedsearch_outputs = expand(f'{output_dir}/Hypedsearch_outputs/HS_{{dataset}}.txt', dataset=base_files)
