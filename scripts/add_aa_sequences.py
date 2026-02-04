import pandas as pd
import re
from collections import defaultdict
import sys

# Inputs from Snakemake
simphyni_csv = snakemake.input.sim_csv
panaroo_faa = snakemake.input.faa
ibd_ecoli_csv = snakemake.input.gene_piv
outfile = snakemake.output.outfile

# --- Helper functions ---

def read_fasta(filepath):
    """Reads FASTA into a dictionary: {HeaderID: Sequence}"""
    seqs = {}
    with open(filepath) as f:
        seq_id = None
        seq = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id:
                    seqs[seq_id] = "".join(seq)
                # Panaroo FASTA headers often look like ">group_1234 gene_name"
                # We take just the first part (ID) to match the CSV cluster names
                seq_id = line[1:].split()[0]
                seq = []
            else:
                seq.append(line)
        if seq_id:
            seqs[seq_id] = "".join(seq)
    return seqs

def reformat_string_for_filepath(s):
    replacements = {
        ' ': '_', '\\': '', '/': '', ':': '', '*': '',
        '?': '', '"': '', '<': '', '>': '', '|': '', '.':'_', '~':''
    }
    for k, v in replacements.items():
        s = s.replace(k, v)
    return re.sub(r'[^a-zA-Z0-9_.-]', '', s)

# --- Load data ---
print(f"📖 Loading data...")
ibd_ecoli = pd.read_csv(ibd_ecoli_csv)
simphyni = pd.read_csv(simphyni_csv)
ibd_sequences = read_fasta(panaroo_faa)

# --- Build Mapping ---
# We create a lookup: {SanitizedName: OriginalName}
# Note: If two different genes sanitize to the same name, this picks the first one.
name_lookup = {reformat_string_for_filepath(c): c for c in ibd_ecoli.columns}

# --- Apply Mapping ---

simphyni["T1_original_name"] = simphyni["T1"].apply(lambda x: name_lookup.get(x, None))
simphyni["mapping_status_T1"] = simphyni["T1_original_name"].apply(lambda x: "missing" if pd.isna(x) else "found")

simphyni["T2_original_name"] = simphyni["T2"].apply(lambda x: name_lookup.get(x, None))
simphyni["mapping_status_T2"] = simphyni["T2_original_name"].apply(lambda x: "missing" if pd.isna(x) else "found")

# Add amino acid sequences (mark missing if not found)
def get_sequence(name):
    if pd.isna(name):
        return "no_sequence_found"
    return ibd_sequences.get(name, "no_sequence_found")

simphyni["T1_aa_sequence"] = simphyni["T1_original_name"].apply(get_sequence)
simphyni["T2_aa_sequence"] = simphyni["T2_original_name"].apply(get_sequence)


# --- Save output ---
simphyni.to_csv(outfile, index=False)
print(f"💾 Saved to: {outfile}")