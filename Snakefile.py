###################################
# SNAKEFILE FOR Fasta to SimPhyNI #
###################################

import pandas as pd
import os

# --- Configuration ---
# Options: "poppunk", "raxml", "user"
# You can override this via command line: snakemake --config tree_method=raxml
TREE_METHOD = config.get("tree_method", "poppunk")

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

        # for genomic distances between pairs if wanted
        expand("results/{analysis}/simphyni/{analysis}/distances.csv", analysis=ANALYSES),

# --- Prokka ---
rule prokka:
    input:
        get_fasta
    output:
        gff = "results/{analysis}/prokka/{sample}/{sample}.gff",
        faa = "results/{analysis}/prokka/{sample}/{sample}.faa"
    params:
        outdir = "results/{analysis}/prokka/{sample}",
        prefix = "{sample}"
    conda: "envs/prokka.yaml"
    threads: 8
    shell:
        "rm -rf {params.outdir} && "
        "prokka --outdir {params.outdir} --prefix {params.prefix} --cpus {threads} {input}"

# --- Panaroo (Conditional Alignment) ---
rule panaroo:
    input:
        gffs = lambda wildcards: expand(
            "results/{analysis}/prokka/{sample}/{sample}.gff",
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
    conda: "envs/panaroo.yaml"
    threads: 32
    shell:
        "panaroo -i {input.gffs} -o {output.folder} -t {threads} --clean-mode strict --merge_paralogs --remove-invalid-genes {params.extra_args}"

rule to_faa:
    input: "results/{analysis}/panaroo_results/pan_genome_reference.fa"
    output: "results/{analysis}/panaroo_results/pan_genome_reference.faa"
    conda: "envs/panaroo.yaml"
    shell: "python scripts/fa_to_faa.py {input} {output}"

# --- RAxML-NG ---
rule raxml_ng:
    input:
        aln = "results/{analysis}/panaroo_results/core_gene_alignment_filtered.aln"
    output:
        tree = "results/{analysis}/raxml/{analysis}.raxml.bestTree"
    params:
        prefix = "results/{analysis}/raxml/{analysis}",
        dir = "results/{analysis}/raxml"
    conda: "envs/raxml.yaml" 
    shell:
        """
        mkdir -p {params.dir}
        raxml-ng --msa {input.aln} --model GTR+G --prefix {params.prefix} --threads {threads} --seed 12345
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
    conda: "envs/poppunk.yaml"
    shell:
        """
        rm -rf {output.db} {output.model}
        mkdir -p results/{wildcards.analysis}/poppunk
        poppunk --create-db --r-files {input.r_files} --output {output.db} --threads {threads} --no-plot
        poppunk --fit-model dbscan --ref-db {output.db} --output {output.model} --threads {threads}
        """

rule poppunk_tree:
    input:
        db = "results/{analysis}/poppunk/db",
        model = "results/{analysis}/poppunk/model"
    output:
        tree = "results/{analysis}/poppunk_tree/poppunk_tree_core_NJ.nwk"
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
            --microreact
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
    run:
        from Bio import Phylo
        tree = Phylo.read(input.tree, "newick")
        tree.root_at_midpoint()
        Phylo.write(tree, output.rooted, "newick")

# --- SimPhyNI Prep & Execution ---

rule pan_to_sim:
    input:
        pan_csv = "results/{analysis}/panaroo_results/gene_presence_absence.csv",
        pheno_file = get_optional_pheno
    output:
        piv = "results/{analysis}/simphyni/gene_piv.csv"
    params:
        pheno_arg = lambda wildcards, input: f"--pheno {input.pheno_file}" if input.pheno_file else ""
    conda: "envs/panaroo.yaml"
    shell:
        "python scripts/panaroo_to_simphyni.py "
        "--panaroo {input.pan_csv} "
        "{params.pheno_arg} "
        "--output {output.piv}"

rule simphyni:
    input:
        csv = "results/{analysis}/simphyni/gene_piv.csv",
        tree = get_tree_file # This now auto-switches based on config
    output:
        "results/{analysis}/simphyni/{analysis}/simphyni_results.csv"
    params:
        r_flag = get_simphyni_r_flag
    conda: "envs/simphyni.yaml"
    shell:
        "simphyni run -T {input.tree} -t {input.csv} -s {wildcards.analysis} {params.r_flag} --outdir results/{wildcards.analysis}/simphyni"

rule add_aa:
    input:
        sim_csv = "results/{analysis}/simphyni/{analysis}/simphyni_results.csv",
        faa = "results/{analysis}/panaroo_results/pan_genome_reference.faa",
        gene_piv = "results/{analysis}/simphyni/gene_piv.csv"
    output:
        outfile = "results/{analysis}/simphyni/{analysis}/results_AA.csv"
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
        gff_dir = "results/{analysis}/prokka"
    threads: 32
    conda: "envs/panaroo.yaml",
    shell:
        "python -u scripts/cluster_distance_summary.py {input.sim_csv} {input.panaroo} {params.gff_dir} {output.outfile} --threads {threads}"




