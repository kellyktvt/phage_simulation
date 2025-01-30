from Bio import Entrez, SeqIO
import pinetree as pt
import sys
import os

CELL_VOLUME = 1.1e-15
PHI10_BIND = 1.82e7  # Binding constant for phi10

IGNORE_REGULATORY = ["E. coli promoter E[6]",
                     "T7 promoter phiOR",
                     "T7 promoter phiOL",
                     "E. coli promoter A0 (leftward)"]

IGNORE_GENES = ["gene 10B",
                "possible gene 5.5-5.7",
                "gene 1.5",
                "gene 1.6",
                "gene 4.1",
                "gene 4B",
                "gene 0.6A",
                "gene 0.6B",
                "possible gene 0.6B",
                "gene 0.5",
                "gene 0.4"]

RELABEL_GENES = {"gene 2": "gp-2",
                 "gene 1": "rnapol-1",
                 "gene 3.5": "lysozyme-3.5",
                 "gene 0.7": "protein_kinase-0.7"}

# Optimal E. coli codons
OPT_CODONS_E_COLI = {'A': ['GCT'],
                     'R': ['CGT', 'CGC'],
                     'N': ['AAC'],
                     'D': ['GAC'],
                     'C': ['TGC'],
                     'Q': ['CAG'],
                     'E': ['GAA'],
                     'G': ['GGT', 'GGC'],
                     'H': ['CAC'],
                     'I': ['ATC'],
                     'L': ['CTG'],
                     'F': ['TTC'],
                     'P': ['CCG'],
                     'S': ['TCT', 'TCC'],
                     'T': ['ACT', 'ACC'],
                     'Y': ['TAC'],
                     'V': ['GTT', 'GTA']}


def get_promoter_interactions(name):
    '''
    Calculate promoter binding strengths. The relative strengths defined here
    come from 2012 Covert, et al paper.
    '''
    ecoli_strong = ["E. coli promoter A1",
                    "E. coli promoter A2",
                    "E. coli promoter A3"]
    ecoli_weak = ["E. coli B promoter",
                  "E. coli C promoter"]
    phi1_3 = ["T7 promoter phi1.1A",
              "T7 promoter phi1.1B",
              "T7 promoter phi1.3",
              "T7 promoter phi1.5",
              "T7 promoter phi1.6"]
    phi3_8 = ["T7 promoter phi2.5",
              "T7 promoter phi3.8",
              "T7 promoter phi4c",
              "T7 promoter phi4.3",
              "T7 promoter phi4.7"]
    phi6_5 = ["T7 promoter phi6.5"]
    phi9 = ["T7 promoter phi9"]
    phi10 = ["T7 promoter phi10"]
    phi13 = ["T7 promoter phi13",
             "T7 promoter phi17"]

    if name in ecoli_strong:
        return {'ecolipol': 10e4,
                'ecolipol-p': 3e4}
    elif name in ecoli_weak:
        return {'ecolipol': 1e4,
                'ecolipol-p': 0.3e4}
    elif name in phi1_3:
        return {'rnapol-1': PHI10_BIND * 0.01,
                'rnapol-3.5': PHI10_BIND * 0.01 * 0.5}
    elif name in phi3_8:
        return {'rnapol-1': PHI10_BIND * 0.01,
                'rnapol-3.5': PHI10_BIND * 0.01 * 0.5}
    elif name in phi6_5:
        return {'rnapol-1': PHI10_BIND * 0.05,
                'rnapol-3.5': PHI10_BIND * 0.05}
    elif name in phi9:
        return {'rnapol-1': PHI10_BIND * 0.2,
                'rnapol-3.5': PHI10_BIND * 0.2}
    elif name in phi10:
        return {'rnapol-1': PHI10_BIND,
                'rnapol-3.5': PHI10_BIND}
    elif name in phi13:
        return {'rnapol-1': PHI10_BIND * 0.1,
                'rnapol-3.5': PHI10_BIND * 0.1}
    else:
        raise ValueError(
            "Promoter strength for {0} not assigned.".format(name))


def get_terminator_interactions(name):
    '''
    Get terminator efficiencies.
    '''
    if name == "E. coli transcription terminator TE":
        return {'ecolipol': 1.0,
                'ecolipol-p': 1.0,
                'rnapol-1': 0.0,
                'rnapol-3.5': 0.0}
    elif name == "T7 transcription terminator Tphi":
        return {'rnapol-1': 0.85,
                'rnapol-3.5': 0.85}
    else:
        return {'name': 0.0}


def compute_cds_weights(record, feature, opt_weight, nonopt_weight, weights):
    # Extract nucleotide sequence
    nuc_seq = feature.location.extract(record).seq
    # Extract translated amino acid sequence
    aa_seq = feature.qualifiers["translation"][0]
    
    # Iterate over nucleotide sequence
    for index, nuc in enumerate(nuc_seq):
        # Calculate amino acid index (Divide nucleotide index by 3 since an amino acid is 3 nucleotides aka 1 codon)
        aa_index = int(index / 3)
        # Determine start index of codon
        codon_start = aa_index * 3
        # Extract codon from nucleotide sequence
        codon = nuc_seq[codon_start:codon_start + 3]
        # Calculate genome index of nucleotide, accounting feature's location
        genome_index = feature.location.start + index

        # Check if amino acid seq is long enough to have an amino acid at the current index
        if aa_index < len(aa_seq):
            # Check if amino acid at current index has optimal codons
            if aa_seq[aa_index] in OPT_CODONS_E_COLI and codon in OPT_CODONS_E_COLI[aa_seq[aa_index]]:
                weights[genome_index] = opt_weight
            else:
                weights[genome_index] = nonopt_weight

    return weights


def normalize_weights(weights):
    # Average over all CDSs, which will have non-zero weights
    non_zero = sum(1 if i != 0 else 0.0 for i in weights)
    mean_weight = sum(weights) / non_zero
    norm_weights = [i / mean_weight for i in weights]
    # Replace non-CDS weights with 1
    norm_weights = [1 if i == 0 else i for i in norm_weights]
    return norm_weights


def main(opt_weight, nonopt_weight):
    sim = pt.Model(cell_volume=CELL_VOLUME)

    # Download T7 wild-type genbank records
    Entrez.email = "kelly.to@utexas.edu"
    handle = Entrez.efetch(db="nuccore",
                           id=["NC_001604"],
                           rettype="gb",
                           retmode="text")

    record = SeqIO.read(handle, "genbank")
    genome_length = len(record.seq)
    phage = pt.Genome(name="phage", length=genome_length)

    for feature in record.features:
        weights = [nonopt_weight] * len(record.seq)
        # Convert to inclusive genomic coordinates
        start = feature.location.start + 1
        stop = feature.location.end
        name = ''
        if "note" in feature.qualifiers:
            name = feature.qualifiers["note"][0]
        # Grab promoters and terminators
        if feature.type == "regulatory":
            if name in IGNORE_REGULATORY:
                continue
            # Construct promoter
            if "promoter" in feature.qualifiers["regulatory_class"]:
                length = stop - start
                if length < 35:
                    start = start - 35
                interactions = get_promoter_interactions(name)
                phage.add_promoter(name, start, stop, interactions)
            # Construct terminator params
            if "terminator" in feature.qualifiers["regulatory_class"]:
                interactions = get_terminator_interactions(name)
                phage.add_terminator(name, start, stop, interactions)
        # Grab genes/CDSes
        if feature.type == "gene":
            if name in IGNORE_GENES:
                continue
            if name in RELABEL_GENES:
                name = RELABEL_GENES[name]
            # Construct CDS parameters for this gene
            phage.add_gene(name=name, start=start, stop=stop,
                           rbs_start=start - 30, rbs_stop=start, rbs_strength=1e7)
        if feature.type == "CDS":
            weights = compute_cds_weights(record, feature, opt_weight, nonopt_weight, weights)

    mask_interactions = ["rnapol-1", "rnapol-3.5",
                         "ecolipol", "ecolipol-p", "ecolipol-2", "ecolipol-2-p"]
    phage.add_mask(500, mask_interactions)

    norm_weights = normalize_weights(weights)
#-------------------------------------------------------------------------------

    # commenting out the use of weights allows you to run this sim on the dev version of pinetree
    phage.add_weights(norm_weights)

#-------------------------------------------------------------------------------

    sim.register_genome(phage)

    sim.add_polymerase("rnapol-1", 35, 230, 0)
    sim.add_polymerase("rnapol-3.5", 35, 230, 0)
    sim.add_polymerase("ecolipol", 35, 45, 0)
    sim.add_polymerase("ecolipol-p", 35, 45, 0)
    sim.add_polymerase("ecolipol-2", 35, 45, 0)
    sim.add_polymerase("ecolipol-2-p", 35, 45, 0)

    sim.add_ribosome(30, 30, 0)

    sim.add_species("bound_ribosome", 10000)

    sim.add_species("bound_ecolipol", 1800)
    sim.add_species("bound_ecolipol_p", 0)
    sim.add_species("ecoli_genome", 0)
    sim.add_species("ecoli_transcript", 0)

    sim.add_reaction(1e6, ["ecoli_transcript", "__ribosome"], [
                     "bound_ribosome"])

    sim.add_reaction(0.04, ["bound_ribosome"], [
                     "__ribosome", "ecoli_transcript"])

    sim.add_reaction(0.001925, ["ecoli_transcript"], ["degraded_transcript"])

    sim.add_reaction(1e7, ["ecolipol", "ecoli_genome"], ["bound_ecolipol"])

    sim.add_reaction(
        0.3e7, ["ecolipol-p", "ecoli_genome"], ["bound_ecolipol_p"])

    sim.add_reaction(0.04, ["bound_ecolipol"], [
                     "ecolipol", "ecoli_genome", "ecoli_transcript"])

    sim.add_reaction(0.04, ["bound_ecolipol_p"], [
                     "ecolipol-p", "ecoli_genome", "ecoli_transcript"])

    sim.add_reaction(3.8e7, ["protein_kinase-0.7", "ecolipol"],
                     ["ecolipol-p", "protein_kinase-0.7"])

    sim.add_reaction(3.8e7, ["protein_kinase-0.7", "ecolipol-2"],
                     ["ecolipol-2-p", "protein_kinase-0.7"])

    sim.add_reaction(3.8e7, ["gp-2", "ecolipol"], ["ecolipol-2"])

    sim.add_reaction(3.8e7, ["gp-2", "ecolipol-p"], ["ecolipol-2-p"])

    sim.add_reaction(1.1, ["ecolipol-2-p"], ["gp-2", "ecolipol-p"])

    sim.add_reaction(1.1, ["ecolipol-2"], ["gp-2", "ecolipol"])

    sim.add_reaction(3.8e9, ["lysozyme-3.5", "rnapol-1"], ["rnapol-3.5"])

    sim.add_reaction(3.5, ["rnapol-3.5"], ["lysozyme-3.5", "rnapol-1"])

    sim.seed(seed_val)

    # generate a unique filename based on the charge
    output_dir = f"data/simulation/phage/revised_weighted_opt1.0_nonopt0.5"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    output_filename = os.path.join(output_dir, f"revised_weighted_opt{opt_weight}_nonopt{nonopt_weight}_seed{seed_val}.tsv")

    print("Starting simulation...")
    try:
        sim.simulate(time_limit=1200, time_step=5, 
                    output=output_filename)
        print("Simulation completed.")
    except Exception as e:
        print(f"Simulation failed: {e}")


if __name__ == "__main__":
    # Take in arguments from the jobs file
    # optimal codon weight 
    opt_weight = float(sys.argv[1])
    # nonoptimal codon weight
    nonopt_weight = float(sys.argv[2])
    # seed value
    seed_val = int(sys.argv[3])

    main(opt_weight, nonopt_weight)