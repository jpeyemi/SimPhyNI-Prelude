###################################
# SNAKEFILE FOR Fasta to SimPhyNI #
###################################

import pandas as pd
import os

# --- Configuration ---
# Options: "poppunk", "raxml", "user"
# You can override this via command line: snakemake --config tree_method=raxml
TREE_METHOD = config.get("tree_method", "poppunk")

# Set to False to skip genomic distance calculation: --config calc_distances=False
CALC_DISTANCES = config.get("calc_distances", True)

# Path to Bakta annotation database.
# Option 1 — point to an existing shared DB on your HPC:
#   snakemake --config bakta_db=/shared/databases/bakta/db
# Option 2 — download the light DB locally first by running:
#   snakemake bakta_db_download
BAKTA_DB = config.get("bakta_db", "resources/bakta_db/db-light")

# --- Load Metadata ---
df = pd.read_csv("inputs/master_metadata.csv")

# Lists all analyses names corresponding to entries in `master_metadata.csv`
ANALYSES = df['analysis_id'].unique().tolist() # If not all specify which analyses
 
# --- Helper Functions ---

def get_samples_by_analysis(analysis_name):
    return df[df['analysis_id'] == analysis_name]['sample_id'].tolist()

def get_fasta(wildcards):
    return df[df['sample_id'] == wildcards.sample]['fasta_path'].values[0]

def get_fastas_by_analysis(analysis_name):
    return df[df['analysis_id'] == analysis_name]['fasta_path'].tolist()

def get_optional_pheno(wildcards):
    pheno_path = f"inputs/{wildcards.analysis}_phenotype.csv"
    if os.path.exists(pheno_path):
        return pheno_path
    return []

def get_simphyni_r_flag(wildcards):
    pheno_path = f"inputs/{wildcards.analysis}_phenotype.csv"
    if os.path.exists(pheno_path):
        try:
            df_pheno = pd.read_csv(pheno_path, nrows=0)
            num_phenotypes = len(df_pheno.columns) - 1
            if num_phenotypes > 0:
                indices = ",".join([str(i) for i in range(num_phenotypes)])
                return f"-r {indices}"
        except Exception as e:
            print(f"Warning reading phenotype file: {e}")
    return "-r ALL"

# --- Tree Logic ---
def get_tree_file(wildcards):
    """
    Determines which tree to use based on the TREE_METHOD config.
    """
    if TREE_METHOD == "user":
        # 1. User provided tree
        return f"inputs/{wildcards.analysis}_tree.nwk"
    
    elif TREE_METHOD == "raxml":
        # 2. RAxML-NG tree (High Resolution)
        # return f"results/{wildcards.analysis}/raxml/{wildcards.analysis}.raxml.bestTree"
        return f"results/{wildcards.analysis}/raxml/{wildcards.analysis}_rooted.nwk"
    
    elif TREE_METHOD == "poppunk":
        # 3. PopPUNK tree (Fast/Low Res)
        # return f"results/{wildcards.analysis}/poppunk_tree/poppunk_tree_core_NJ.nwk"
        return f"results/{wildcards.analysis}/poppunk_tree/{wildcards.analysis}_rooted.nwk"
    
    else:
        raise ValueError(f"Invalid tree_method: {TREE_METHOD}. Choose 'raxml', 'poppunk', or 'user'.")


# --- Rules ---

rule all:
    input:
        expand("results/{analysis}/simphyni/{analysis}/results_AA.csv", analysis=ANALYSES),

        # for genomic distances between pairs — skip with --config calc_distances=False
        expand("results/{analysis}/simphyni/{analysis}/distances.csv", analysis=ANALYSES) if CALC_DISTANCES else [],

# --- Bakta ---
rule bakta:
    input:
        get_fasta
    output:
        gff = "results/{analysis}/bakta/{sample}/{sample}.gff3",
        faa = "results/{analysis}/bakta/{sample}/{sample}.faa"
    params:
        outdir = "results/{analysis}/bakta/{sample}",
        prefix = "{sample}",
        db = BAKTA_DB
    log: "logs/{analysis}/bakta_{sample}.log"
    conda: "envs/bakta.yaml"
    threads: 8
    shell:
        "bakta --db {params.db} --output {params.outdir} --prefix {params.prefix} "
        "--threads {threads} --force {input} > {log} 2>&1"

rule bakta_db_download:
    output:
        db = directory("resources/bakta_db/db-light")
    log: "logs/bakta_db_download.log"
    conda: "envs/bakta.yaml"
    shell:
        "bakta_db download --output resources/bakta_db --type light > {log} 2>&1"

# --- Panaroo (Conditional Alignment) ---
rule panaroo:
    input:
        gffs = lambda wildcards: expand(
            "results/{analysis}/bakta/{sample}/{sample}.gff3",
            analysis=wildcards.analysis,
            sample=get_samples_by_analysis(wildcards.analysis)
        )
    output:
        folder = directory("results/{analysis}/panaroo_results"),
        csv = "results/{analysis}/panaroo_results/gene_presence_absence.csv",
        ref = "results/{analysis}/panaroo_results/pan_genome_reference.fa",
        # Only expect an alignment file if we are running RAxML
        aln = "results/{analysis}/panaroo_results/core_gene_alignment_filtered.aln" if TREE_METHOD == "raxml" else []
    params:
        # Dynamically add alignment flags only if needed
        extra_args = "-a core --aligner mafft --core_threshold 0.95" if TREE_METHOD == "raxml" else ""
    log: "logs/{analysis}/panaroo.log"
    conda: "envs/panaroo.yaml"
    threads: 32
    shell:
        "panaroo -i {input.gffs} -o {output.folder} -t {threads} --clean-mode strict --merge_paralogs --remove-invalid-genes {params.extra_args} > {log} 2>&1"

rule to_faa:
    input: "results/{analysis}/panaroo_results/pan_genome_reference.fa"
    output: "results/{analysis}/panaroo_results/pan_genome_reference.faa"
    log: "logs/{analysis}/to_faa.log"
    conda: "envs/panaroo.yaml"
    shell: "python scripts/fa_to_faa.py {input} {output} > {log} 2>&1"

# --- RAxML-NG ---
rule raxml_ng:
    input:
        aln = "results/{analysis}/panaroo_results/core_gene_alignment_filtered.aln"
    output:
        tree = "results/{analysis}/raxml/{analysis}.raxml.bestTree"
    params:
        prefix = "results/{analysis}/raxml/{analysis}",
        dir = "results/{analysis}/raxml"
    log: "logs/{analysis}/raxml_ng.log"
    conda: "envs/raxml.yaml"
    threads: 16
    shell:
        """
        mkdir -p {params.dir}
        raxml-ng --msa {input.aln} --model GTR+G --prefix {params.prefix} --threads {threads} --seed 12345 > {log} 2>&1
        """

# --- PopPUNK ---
rule create_poppunk_list:
    input:
        fastas = lambda wildcards: get_fastas_by_analysis(wildcards.analysis)
    output:
        list_file = "results/{analysis}/poppunk/{analysis}_assemblies.txt"
    run:
        samples = get_samples_by_analysis(wildcards.analysis)
        paths = get_fastas_by_analysis(wildcards.analysis)
        with open(output.list_file, "w") as f:
            for s, p in zip(samples, paths):
                f.write(f"{s}\t{p}\n")

rule poppunk_db_fit:
    input:
        r_files = "results/{analysis}/poppunk/{analysis}_assemblies.txt"
    output:
        db = directory("results/{analysis}/poppunk/db"),
        model = directory("results/{analysis}/poppunk/model")
    threads: 64
    resources:
        mem_mb = 32000 # Increased memory for large sketches
    log: "logs/{analysis}/poppunk_db_fit.log"
    conda: "envs/poppunk.yaml"
    shell:
        """
        rm -rf {output.db} {output.model}
        mkdir -p results/{wildcards.analysis}/poppunk
        poppunk --create-db --r-files {input.r_files} --output {output.db} --threads {threads} --no-plot > {log} 2>&1
        poppunk --fit-model dbscan --ref-db {output.db} --output {output.model} --threads {threads} >> {log} 2>&1
        """

rule poppunk_tree:
    input:
        db = "results/{analysis}/poppunk/db",
        model = "results/{analysis}/poppunk/model"
    output:
        tree = "results/{analysis}/poppunk_tree/poppunk_tree_core_NJ.nwk"
    log: "logs/{analysis}/poppunk_tree.log"
    conda: "envs/poppunk.yaml"
    threads: 16
    resources:
        mem_mb = 16000 # Increased memory for tree building
    shell:
        """
        poppunk_visualise \
            --ref-db {input.db} \
            --model-dir {input.model} \
            --output results/{wildcards.analysis}/poppunk_tree \
            --tree nj \
            --threads {threads} \
            --microreact > {log} 2>&1
        """

# --- Add the Midpoint Rooting Rule ---
rule midpoint_root:
    input:
        # This dynamic input allows the rule to work for both RAxML and PopPUNK
        tree = (f"results/{{analysis}}/raxml/{{analysis}}.raxml.bestTree"
                if TREE_METHOD == "raxml" else
                f"results/{{analysis}}/poppunk_tree/poppunk_tree_core_NJ.nwk")
    output:
        rooted = (f"results/{{analysis}}/raxml/{{analysis}}_rooted.nwk"
                  if TREE_METHOD == "raxml" else
                  f"results/{{analysis}}/poppunk_tree/{{analysis}}_rooted.nwk")
    log: "logs/{analysis}/midpoint_root.log"
    conda: "envs/panaroo.yaml"
    script:
        "scripts/midpoint_root.py"

# --- SimPhyNI Prep & Execution ---

rule pan_to_sim:
    input:
        pan_csv = "results/{analysis}/panaroo_results/gene_presence_absence.csv",
        pheno_file = get_optional_pheno
    output:
        piv = "results/{analysis}/simphyni/gene_piv.csv"
    params:
        pheno_arg = lambda wildcards, input: f"--pheno {input.pheno_file}" if input.pheno_file else ""
    log: "logs/{analysis}/pan_to_sim.log"
    conda: "envs/panaroo.yaml"
    shell:
        "python scripts/panaroo_to_simphyni.py "
        "--panaroo {input.pan_csv} "
        "{params.pheno_arg} "
        "--output {output.piv} > {log} 2>&1"

rule simphyni:
    input:
        csv = "results/{analysis}/simphyni/gene_piv.csv",
        tree = get_tree_file # This now auto-switches based on config
    output:
        "results/{analysis}/simphyni/{analysis}/simphyni_results.csv"
    params:
        r_flag = get_simphyni_r_flag
    log: "logs/{analysis}/simphyni.log"
    conda: "envs/simphyni.yaml"
    threads: 64
    shell:
        "simphyni run -T {input.tree} -t {input.csv} -s {wildcards.analysis} {params.r_flag} --outdir results/{wildcards.analysis}/simphyni > {log} 2>&1"

rule add_aa:
    input:
        sim_csv = "results/{analysis}/simphyni/{analysis}/simphyni_results.csv",
        faa = "results/{analysis}/panaroo_results/pan_genome_reference.faa",
        gene_piv = "results/{analysis}/simphyni/gene_piv.csv"
    output:
        outfile = "results/{analysis}/simphyni/{analysis}/results_AA.csv"
    log: "logs/{analysis}/add_aa.log"
    conda: "envs/panaroo.yaml"
    script:
        "scripts/add_aa_sequences.py"

rule calc_genomic_distance:
    input:
        sim_csv = "results/{analysis}/simphyni/{analysis}/results_AA.csv",
        panaroo = "results/{analysis}/panaroo_results/gene_presence_absence.csv",
    output:
        outfile = "results/{analysis}/simphyni/{analysis}/distances.csv"
    params:
        gff_dir = "results/{analysis}/bakta",
        pval_col = config.get("dist_pval_col", "pval_bh"),
        alpha = config.get("dist_alpha", 0.05)
    log: "logs/{analysis}/calc_genomic_distance.log"
    threads: 32
    conda: "envs/panaroo.yaml"
    shell:
        "python -u scripts/cluster_distance_summary.py {input.sim_csv} {input.panaroo} {params.gff_dir} {output.outfile} "
        "--threads {threads} --pval-col {params.pval_col} --alpha {params.alpha} > {log} 2>&1"




