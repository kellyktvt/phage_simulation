#!/usr/bin/env python3
import argparse, sys

from Bio import SeqIO 
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord    

from codon_tools import FopScorer

def scorer(seq, gene_description):
    print(f"Sequence length: {len(seq)}")  # Add this line to debug
    print(f"Sequence: {seq}")  # Debug: print the sequence
    
    if len(seq) % 3 != 0:
        raise ValueError(f"Sequence length is not a multiple of 3: {len(seq)}")
    
    scorer = FopScorer()
    score = scorer.score(seq)
    records = []

    seq_record = SeqRecord(seq, id = '', description = gene_description + " -- Fop = %f" % score)
    records.append(seq_record)
    
    return records
        
# when run as its own script
if __name__ == "__main__":
    # set up command line arguments
    parser = argparse.ArgumentParser(description='Calculate Fop scores for gene sequences.')
    parser.add_argument('filename', metavar='filename',
                        help='fasta-formatted file holding the sequence')
    args = parser.parse_args()

    # read in the sequence     
    record = SeqIO.parse(open(args.filename, "r"), "fasta").__next__()
    # run analysis    
    records = scorer(record.seq, record.description)
    # save the scores to a new FASTA file
    output_filename = "src/python/models/trna_phage_model/id_map_fop_scores.fasta"
    with open(output_filename, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
