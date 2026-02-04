import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def main():
    parser = argparse.ArgumentParser(description="Convert nucleotide FASTA (.fa) to amino acid FASTA (.faa)")
    parser.add_argument("input_file", help="Input nucleotide FASTA file (e.g., pan_genome_reference.fa)")
    parser.add_argument("output_file", help="Output protein FASTA file (e.g., pan_genome_reference.faa)")
    args = parser.parse_args()

    input_file = args.input_file
    output_file = args.output_file

    protein_records = []
    for record in SeqIO.parse(input_file, "fasta"):
        seq = record.seq
        # Ensure length is multiple of 3 to avoid translation errors
        if len(seq) % 3 != 0:
            seq = seq[: len(seq) - (len(seq) % 3)]
        try:
            translated_seq = seq.translate(to_stop=True)
        except Exception as e:
            sys.stderr.write(f"Warning: failed to translate {record.id}: {e}\n")
            continue
        protein_records.append(SeqRecord(translated_seq, id=record.id, description="translated"))

    SeqIO.write(protein_records, output_file, "fasta")
    print(f"✅ Translation complete. Output saved to {output_file} ({len(protein_records)} records).")

if __name__ == "__main__":
    main()
