import csv
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
import time
import urllib.error

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
        if feature.type == "gene" and gene_name in feature.qualifiers["note"]:
            gene_seq = feature.extract(record.seq)
            return gene_seq
    return None

def process_id_map(id_map_file):
    genes = []
    with open(id_map_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene_name = row['pinetree']
            if gene_name:  # Only include rows with gene_name
                genes.append(gene_name)
    return genes

def main():
    Entrez.email = "kelly.to@utexas.edu"
    
    # Process the CSV file to get gene names from the 'pinetree' column
    genes = process_id_map("id_map.csv")

    # Fetch the GenBank record for T7 wild-type 
    record = fetch_genbank_record()

    all_sequences = []
    
    for gene_name in genes:
        # Extract the gene sequence 
        gene_sequence = extract_gene_sequence(record, gene_name)

        # Append gene sequence to the FASTA file
        all_sequences.append(gene_sequence)

    # Save all sequences to a single FASTA file
    with open("all_id_map_sequences.fasta", "w") as fasta_file:
        for i, sequence in enumerate(all_sequences):
            gene_name = genes[i]
            gene_record = SeqRecord(sequence, id=gene_name, description=f"{gene_name} sequence")
            SeqIO.write(gene_record, fasta_file, "fasta")

if __name__ == "__main__":
    main()
