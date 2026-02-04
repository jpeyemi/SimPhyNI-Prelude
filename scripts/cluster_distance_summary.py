import os
import csv
import glob
import argparse
from collections import defaultdict
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor, as_completed
# from tqdm import tqdm

import re
id_pattern = re.compile(r'ID=([^;]+)')

def parse_gff(filepath):
    gene_positions = {}
    contig_size = {}

    with open(filepath) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            contig = parts[0]

            if ftype == "CDS":
                start, end = int(parts[3]), int(parts[4])
                match = id_pattern.search(parts[8])
                if match:
                    gene_id = match.group(1)
                    gene_positions[gene_id] = (contig, (start + end) // 2)

            elif ftype == "region":
                length = int(parts[4])
                attrs = parts[8]
                if "is_circular=true" not in attrs.lower():
                    length = 0
                contig_size[contig] = (contig, length)

    return gene_positions, contig_size

# def load_cluster_pairs(pair_file):
#     pairs = set()
#     with open(pair_file) as f:
#         reader = csv.reader(f)
#         next(reader)  # Skip header
#         for row in reader:
#             if len(row) >= 2:
#                 geneA, geneB = row[0].strip(), row[1].strip()
#                 pairs.add((geneA, geneB))
#     return pairs

def load_cluster_pairs(pair_file, col1, col2, pval_col='pval_bh',alpha = 0.05):
    df = pd.read_csv(pair_file)
    significant_df = df[df[pval_col] <= alpha]
    pairs = set(significant_df[[col1, col2]].itertuples(index=False, name=None))
    return pairs


def prepare_presence_maps(presence_df):
    """Return a dict of genome -> {cluster: gene_id}"""
    presence_maps = {}
    for genome in presence_df.columns[3:]:
        presence_maps[genome] = presence_df[genome].dropna().to_dict()
    return presence_maps

def process_genome(gff_path, cluster_pairs, presence_map):
    genome_name = os.path.basename(gff_path)[:-4]
    gene_positions, contig_size = parse_gff(gff_path)
    genome_distances = defaultdict(list)

    # Fast membership filter
    present_clusters = set(presence_map.keys())
    filtered_pairs = [(A, B) for A, B in cluster_pairs if A in present_clusters and B in present_clusters]

    for clusterA, clusterB in filtered_pairs:
        geneA = presence_map[clusterA]
        geneB = presence_map[clusterB]

        if geneA not in gene_positions or geneB not in gene_positions:
            genome_distances[(clusterA, clusterB)].append(np.nan)
            continue

        contigA, posA = gene_positions[geneA]
        contigB, posB = gene_positions[geneB]

        if contigA != contigB:
            genome_distances[(clusterA, clusterB)].append(np.nan)
            continue

        contig_name, length = contig_size.get(contigA, (None, 0))
        dist = abs(posA - posB)
        if length:
            dist = min(dist, abs(length - dist))
        genome_distances[(clusterA, clusterB)].append(dist)

    return genome_distances

def merge_results(all_results): 
    combined = defaultdict(list)
    for result in all_results:
        for key, dists in result.items():
            if dists:  # Only add non-empty results
                combined[key].extend(dists)
            else:
                print(f"Skipping empty pair: {key}")  # Debugging line
    print(f"Combined entries: {len(combined)}")  # Check the number of unique gene pairs
    return combined

def summarize_and_write(distances, output_csv):
    with open(output_csv, "w", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["geneA", "geneB", "average_distance", "min_distance", "max_distance", "same_contig_count", "total_count"])
        
        for (geneA, geneB), dists in distances.items():
            # Filter out np.nan values and calculate stats only if there are valid distances
            valid_dists = [d for d in dists if not np.isnan(d)]
            
            if valid_dists:  # If there are valid distances (i.e., not all are np.nan)
                avg = round(np.nanmean(valid_dists), 2)
                writer.writerow([geneA, geneB, avg, np.nanmin(valid_dists), np.nanmax(valid_dists), len(valid_dists), len(dists)])
            else:
                # If no valid distances, write NaN for all columns
                writer.writerow([geneA, geneB, np.nan, np.nan, np.nan, 0, len(dists)])


def main(simphyni_file, presence_absence_file, gff_dir, output_csv, num_workers):
    cluster_pairs = load_cluster_pairs(simphyni_file, "T1_original_name","T2_original_name",'pval_bh',0.05)
    presence_df = pd.read_csv(presence_absence_file, dtype=str)
    if presence_df.index.name != 'Gene':
        presence_df.set_index('Gene', inplace=True)

    presence_maps = prepare_presence_maps(presence_df)
    genome_names = list(presence_maps.keys())
    gff_files = [os.path.join(gff_dir, genome, f"{genome}.gff") for genome in genome_names]

    results = []
    with ProcessPoolExecutor(max_workers=num_workers) as pool:
        futures = {pool.submit(process_genome, gff, cluster_pairs, presence_maps[genome]): genome
                   for genome, gff in zip(genome_names, gff_files)}

        for future in as_completed(futures):
            try:
                results.append(future.result())
                print(f"✅ Finished {futures[future]}")
            except Exception as e:
                print(f"Error processing {futures[future]}: {e}")

    combined = merge_results(results)
    summarize_and_write(combined, output_csv)
    print(f"Done! Results written to {output_csv}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate average gene distances for gene cluster pairs.")
    parser.add_argument("simphyni_file", help="CSV output from simphyni")
    parser.add_argument("presence_absence_file", help="Panaroo gene_presence_absence.csv file")
    parser.add_argument("gff_dir", help="Directory containing GFF files")
    parser.add_argument("output_csv", help="Output CSV file path")
    parser.add_argument("--threads", type=int, default=16, help="Number of parallel workers (default: 16)")
    args = parser.parse_args()

    main(args.simphyni_file, args.presence_absence_file, args.gff_dir, args.output_csv, args.threads)
    
