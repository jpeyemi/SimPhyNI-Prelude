import pandas as pd
import argparse
import sys
import os

def main():
    parser = argparse.ArgumentParser(description="Format Panaroo output for Simphyni")
    parser.add_argument("--panaroo", required=True, help="Panaroo gene_presence_absence.csv")
    parser.add_argument("--pheno", required=False, default=None, help="Phenotype CSV file (Optional)")
    parser.add_argument("--output", required=True, help="Output path for Simphyni")
    args = parser.parse_args()

    # 1. Load Panaroo Data
    print(f"Loading Panaroo data: {args.panaroo}")
    pan = pd.read_csv(args.panaroo, low_memory=False)
    
    metadata = [
        'Gene', 'Non-unique Gene name', 'Annotation', 'No. isolates', 
        'No. sequences', 'Avg sequences per isolate', 'Genome Fragment', 
        'Order within Fragment', 'Accessory Fragment', 'Accessory Order with Fragment', 
        'QC', 'Min group size nuc', 'Max group size nuc', 'Avg group size nuc'
    ]
    
    gene_ids = pan['Gene']
    samples_only = pan.drop(columns=[c for c in metadata if c in pan.columns])
    binary_matrix = samples_only.notnull().astype(int)
    binary_matrix.index = gene_ids

    # 2. Transpose
    gene_data = binary_matrix.T
    gene_data.index.name = 'sample_id'
    gene_data.index = gene_data.index.astype(str)

    # 3. Handle Optional Phenotype
    if args.pheno and os.path.exists(args.pheno):
        print(f"Matching phenotypes from: {args.pheno}")
        pheno_df = pd.read_csv(args.pheno, index_col=0)
        pheno_df.index.name = 'sample_id'
        pheno_df.index = pheno_df.index.astype(str)
        
        # Inner Join
        final_df = pheno_df.join(gene_data, how='inner')
        
        if final_df.empty:
            print("ERROR: Phenotype provided but no matching IDs found!")
            sys.exit(1)
    else:
        print("No phenotype file provided or file missing. Exporting gene matrix only.")
        final_df = gene_data

    # 4. Save
    os.makedirs(os.path.dirname(args.output), exist_ok=True)
    final_df.to_csv(args.output)
    print(f"Saved output to: {args.output} ({final_df.shape})")

if __name__ == "__main__":
    main()