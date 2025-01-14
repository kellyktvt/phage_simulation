import csv
import time
import pandas as pd
import urllib.error
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from codon_tools import FopScorer

def fetch_genbank_record():
    for attempt in range(5):
        try:
            handle = Entrez.efetch(db="nuccore", id=["NC_001604"], rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            return record
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(2 ** attempt)  # Exponential backoff
                print(f"Rate limit exceeded. Retrying...")
            else:
                print(f"Failed to fetch GenBank record: {e}")
                break
    return None


def extract_gene_sequence(record, gene_name):
    for feature in record.features:
        if feature.type == "gene" and "note" in feature.qualifiers and gene_name in feature.qualifiers["note"]:
            return feature.extract(record.seq)
    return None


def process_id_map(id_map_file):
    genes = []
    with open(id_map_file, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene_name = row["pinetree"]
            if gene_name:
                genes.append(gene_name)
    return genes


def main():
    Entrez.email = "kelly.to@utexas.edu"
    
    # Process the CSV file to get gene names from the 'pinetree' column
    genes = process_id_map("output/id_map.csv")

    # Fetch the GenBank record for T7 wild-type
    record = fetch_genbank_record()
    if not record:
        print("Failed to fetch GenBank record.")
        return

    data = []
    for gene_name in genes:
        # Extract the gene sequence
        gene_sequence = extract_gene_sequence(record, gene_name)

        if gene_sequence is None:
            print(f"Warning: No sequence found for gene {gene_name}")
            continue

        print(f"{gene_name}: sequence length {len(gene_sequence)}")

        # Calculate Fop score
        fop_score = None
        if len(gene_sequence) % 3 == 0:
            scorer = FopScorer()
            fop_score = scorer.score(gene_sequence)
        else:
            print(f"Skipping {gene_name} due to invalid sequence length.")

        # Append results to the data list
        data.append({
            "gene_name": gene_name,
            "fop_score": fop_score,
            "length": len(gene_sequence)
        })

    # Create a DataFrame and save to a CSV file
    df = pd.DataFrame(data)
    output_filename = "src/python/models/trna_phage_model/fop_scores_and_lengths.csv"
    df.to_csv(output_filename, index=False)
    print(f"FOP scores and lengths saved to {output_filename}.")


if __name__ == "__main__":
    main()
