# SimPhyNI-Prelude

A Snakemake workflow to run **SimPhyNI** starting from genome assemblies (FASTA).

This pipline takes genome assemblies (FASTA) as inputs then performs annotation (Prokka), pangenome inference (Panaroo), and phylogenetic tree construction (PopPUNK or RAxML) to prepare the necessary inputs for SimPhyNI.

> **Note:** This is a wrapper pipeline. For documentation on the core tool and interpreting results, please visit the **[SimPhyNI GitHub Repository](https://github.com/jpeyemi/SimPhyNI)**.

---

## Workflow Overview

1. **Annotation:** Annotates FASTA files using `Prokka`.
2. **Pangenome:** Generates a pangenome and core-gene alignment using `Panaroo`.
3. **Phylogeny:** Builds a rooted phylogenetic tree using either:
* **PopPUNK** (Fast, NJ tree based on core distances)
* **RAxML-NG** (High-resolution ML tree from core alignment)
* **User-Provided** (Skip tree building)


4. **SimPhyNI:** Formats inputs and executes phylogenetically-aware GWAS/correlated evolution analysis.

---

## Installation & Requirements

### 1. Clone Repository

```bash
git clone https://github.com/jpeyemi/SimPhyNI-Prelude.git
cd simphyni-prelude

```

### 2. Dependencies

You need **Conda** and **Snakemake v9.16** installed. Additionally, because this pipeline performs some lightweight Python operations on the head node, your base environment must include `pandas` and `biopython`.


#### Option 1: Lightweight Launcher (Recommended)

Create a small, dedicated environment just for running this pipeline.

```bash
conda create -n simphyni_prelude -c conda-forge -c bioconda \
    snakemake=9.16 pandas biopython

conda activate simphyni_prelude

```

#### Option 2: Use the Full SimPhyNI Environment

If you have already installed the core **SimPhyNI** package, you can simply activate that environment, as it includes the necessary dependencies.

Please refer to the [SimPhyNI GitHub](https://github.com/jpeyemi/SimPhyNI) for installation instructions.

```bash
# Example activation (if you named your env 'simphyni')
conda activate simphyni

```

*All other heavy dependencies (Prokka, Panaroo, RAxML) are handled automatically by Snakemake via Conda environments.*

---

## Input Formatting

The pipeline relies on a specific directory structure and metadata file.

### 1. Master Metadata (`inputs/master_metadata.csv`)

This file controls which samples are processed. It **must** contain the following columns:

| sample_id | analysis_id | fasta_path |
| --- | --- | --- |
| Sample_A | Analysis1 | /path/to/genomes/sampleA.fasta |
| Sample_B | Analysis1 | /path/to/genomes/sampleB.fasta |
| Sample_C | Analysis2 | /path/to/genomes/sampleC.fasta |
| ... | ... | ... |

* **sample_id**: Unique name for the isolate.
* **analysis_id**: A grouping ID (e.g., species or clade). All samples with the same `analysis_id` will be processed together.
* **fasta_path**: Absolute path to the assembly file.

### 2. Phenotypes (Optional)

If you have phenotype data, create a CSV file named `inputs/{analysis_id}_phenotype.csv` for each phenotype analysis.

* Rows should match `sample_id`.
* Columns are binary traits (0/1).

| sample_id | Phenotype1 | Phenotype2 | ... |
| --- | --- | --- | --- |
| Sample_A | 1 | 1 | ... |
| Sample_B | 1 | 0 | ... |
| ... | ... | ...| ... | 

> **Note:** This pipeline does not handle missing phenotype annotations. They will be treated as absences. If this behavior is not desired, its recommended that you exclude samples with missing phenotypes

---

## Configuration

You can control the pipeline behavior by editing `config.yaml` or passing arguments via the command line.

### Specify Analyses

SimPhyNI-Prelude supports simultaneous analyses using the `analysis_id` column in `master_metadata.csv`. By default all analyses will be run, but the `ANALYSES` variable can be set in `Snakefile.py` to specify a subset.

### Tree Methods

You can choose how the phylogenetic tree is generated:

1. **`poppunk`** (Default): Very fast. Uses PopPUNK core distances to build a Neighbor-Joining tree.
2. **`raxml`** (Recommended): Slower but more accurate. Uses Panaroo core gene alignment to build a Maximum Likelihood tree.
3. **`user`**: Uses a pre-existing tree file located at `inputs/{analysis_id}_tree.nwk`.

---

## Usage

### 1. Local Execution

You can run this pipeline locally using standard snakemake arguments

**Basic local run:**

```bash
# Run locally with 8 cores
snakemake -s Snakefile.py --cores 8 --use-conda --config tree_method=user

```


### 2. Running on HPC (Recommended)

To run this pipeline on a High Performance Computing (HPC) cluster, you must install the appropriate Snakemake executor plugin and submit the master Snakemake process as a job.

#### Step 1: Install Executor Plugin

You need to install the executor plugin that matches your cluster's scheduler. See the [Snakemake Plugin Catalog](https://snakemake.github.io/snakemake-plugin-catalog/index.html) for a full list.

For **SLURM**, run:

```bash
pip install snakemake-executor-plugin-slurm

```

#### Step 2: Configure Cluster Profile

The repository includes a default profile in the `cluster_profile/` directory. You may need to edit `cluster_profile/config.yaml` to update the **partition**, **account**, or **resource limits** for your specific cluster environment.

#### Step 3: Submit Pipeline

To avoid timeouts on login or interactive nodes, edit the script `run_simphyni-prelude.sh` for your computing cluster submit it to the scheduler. This script launches the rule: all Snakemake process, which will then submit individual tasks (Prokka, Panaroo, etc.) to the cluster automatically.

**Submit the job:**

```bash
sbatch run_simphyni-prelude.sh

```

---

## Output Structure

Results are organized by `analysis_id` inside the `results/` folder.

```text
results/
├── {analysis_id}/
│   ├── prokka/               # Prokka annotations per sample
│   ├── panaroo_results/      # Panaroo graph and pangenome files
│   ├── poppunk/              # (If selected) PopPUNK db and model
│   ├── poppunk_tree/         # (If selected) Rooted NJ tree
│   ├── raxml/                # (If selected) RAxML best tree
│   └── simphyni/
│       ├── gene_piv.csv      # Formatted input for SimPhyNI
│       └── {analysis_id}/
│           ├── simphyni_results.csv  # MAIN OUTPUT: All association test results
│           └── results_AA.csv        # Results annotated with AA sequences

```

---

## Troubleshooting

* **Biopython Error:** If you see `ModuleNotFoundError: No module named 'Bio'`, ensure you installed `biopython` in your active conda environment (see Installation section).
* **Missing Inputs:** Ensure `inputs/master_metadata.csv` exists and paths in the `fasta_path` column are absolute, not relative.
* **Cluster Timeouts:** If jobs fail due to time limits, edit `cluster_profile/config.yaml` to increase the `runtime` for the specific rule (e.g., `raxml_ng`).
    * Note that `raxml_ng` specifically has builtin checkpointing so it will pick up where it left off

---

## Contact

For questions, please open an issue or contact Ishaq Balogun at https://github.com/jpeyemi.

---

## References

If you use this pipeline, please cite the underlying tools:

* **SimPhyNI:** [Balogun I. (2025)](https://www.google.com/search?q=https://github.com/viba1/SimPhyNI)
* **Prokka:** [Seemann T. (2014)](https://doi.org/10.1093/bioinformatics/btu153)
* **Panaroo:** [Tonkin-Hill G, et al. (2020)](https://doi.org/10.1186/s13059-020-02090-4)
* **RAxML-NG:** [Kozlov AM, et al. (2019)](https://doi.org/10.1093/bioinformatics/btz305)
* **PopPUNK:** [Lees JA, et al. (2019)](https://doi.org/10.1101/gr.241455.118)
* **Snakemake:** [Mölder F, et al. (2021)](https://doi.org/10.12688/f1000research.29032.2)
