#!/usr/bin/env python3
import argparse, sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord    
    
from codon_tools import CodonOptimizer, FopScorer, opt_codons_E_coli


def deoptimize(seq, gene_description, Fop_targets, start_window = 0, end_window = 0, max_wait_count = 5000, opt_codons = opt_codons_E_coli):
    """Parameters:
* ``seq``: coding sequence to de-optimize
* ``start_window``: number of codons to omit at beginning of sequence
* ``end_window``: number of codons to omit at end of sequence
* ``opt_codons``: set of optimal codons

"""
    assert len(seq) % 3 == 0
    
    scorer = FopScorer()
    o = CodonOptimizer(scorer)

    score = scorer.score(seq)
    seq_orig = seq
    Fop_orig = score
    tolerance = 0.01 # how closely do we want to match the final number?
    records = [SeqRecord(seq, id='',
                 description = gene_description + " -- Fop = %f" % Fop_orig)]
    
    for Fop_target in Fop_targets:
        if Fop_target < score:
            maximize = False
        else:
            maximize = True
        seq, score = o.hillclimb(seq, start_window, end_window,
                        Fop_target, tolerance, max_wait_count,
                        maximize, verbosity = 0)
        assert seq_orig.translate() == seq.translate()
        
        description = gene_description + " -- recoded to Fop = %f (keeping first %i and last %i codons unchanged)" % (score, start_window, end_window) 
        seq_record = SeqRecord(seq, id = '', description = description)
        records.append(seq_record)
    
    return records
    
        
# when run as its own script
if __name__ == "__main__":
    # set up command line arguments
    parser = argparse.ArgumentParser(description='Codon-deoptimize gene sequence.')
    parser.add_argument('filename', metavar='filename',
                        help='fasta-formatted file holding the sequence')
    parser.add_argument('-x', '--exclude-front', default=14, type=int,
                        metavar='n',
                        help='number of codons to exclude from the front of the sequence, default = 14')
    parser.add_argument('-X', '--exclude-back', default=14, type=int,
                        metavar='n',
                        help='number of codons to exclude from the back of the sequence, default = 14')
    parser.add_argument('-M', '--max-wait', default=5000, type=int,
                        metavar='n',
                        help='maximum number of attempted codon substitutions before we give up the optimization procedure, default = 5000')


    args = parser.parse_args()

    Fop_targets = [0, 1]

    # read in the sequence     
    record = SeqIO.parse(open(args.filename, "r"), "fasta").__next__()
    # run analysis    
    records = deoptimize(record.seq, record.description,
                         Fop_targets,
                         args.exclude_front, args.exclude_back,
                         args.max_wait)
    # save the deoptimized sequences to a new FASTA file
    output_filename = "gene_10A_deoptimized_extend.fasta"
    with open(output_filename, "w") as output_handle:
        SeqIO.write(records, output_handle, "fasta")
