from Bio import Entrez, SeqIO
import pinetree as pt
import time
import urllib.error
import random
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

#----------------------------------------------------------------------------

def tRNA_map_maker():
    # https://github.com/clauswilke/codon_tools/blob/56491dc51bad4d1bc7bd748b0bf69948843214db/codon_tools/lookup_tables.py
    # Optimal codons for E. coli
    # T. Zhou, M. Weems, C. O. Wilke (2009). Translationally optimal codons associate with structurally sensitive sites in proteins. Mol. Biol. Evol. 26:1571–1580. 
    opt_codons_E_coli = { 'A':['GCT'], 'R':['CGT', 'CGC'], 'N':['AAC'], 'D':['GAC'], 'C':['TGC'], 'Q':['CAG'], 
                         'E':['GAA'], 'G':['GGT','GGC'], 'H':['CAC'], 'I':['ATC'], 'L':['CTG'], 'F':['TTC'], 
                         'P':['CCG'], 'S':['TCT','TCC'], 'T':['ACT','ACC'], 'Y':['TAC'], 'V':['GTT','GTA'] }
    # Reverse genetic code
    # For each amino acid, provides a list of all the codons that translate into that amino acid
    reverse_genetic_code = {  'A':['GCA', 'GCC', 'GCG', 'GCT'], 'R':['AGA', 'AGG', 'CGA', 'CGT', 'CGC', 'CGG'], 
                            'N':['AAC', 'AAT'], 'D':['GAC', 'GAT'], 'C':['TGC', 'TGT' ], 'Q':['CAA', 'CAG'], 
                            'E':['GAA', 'GAG'], 'G':['GGA', 'GGC', 'GGG', 'GGT'], 'H':['CAC', 'CAT'], 
                            'I':['ATA', 'ATC', 'ATT'], 'L':['CTA', 'CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 
                            'F':['TTT', 'TTC'], 'P':['CCA', 'CCC', 'CCG', 'CCT'], 
                            'S':['AGC', 'AGT', 'TCA', 'TCC','TCG','TCT'], 'T':[ 'ACA', 'ACT', 'ACC', 'ACG' ], 
                            'Y':['TAC', 'TAT'], 'V':['GTA', 'GTC', 'GTG', 'GTT'], 'W':['TGG'], 'M':['ATG'], 
                            'K':['AAA', 'AAG'], '*':['TAA', 'TAG', 'TGA']}
    
    # Initialize tRNA_map
    tRNA_map = {}

    # Create tRNA_map
    for amino_acid, codons in reverse_genetic_code.items():
        for codon in codons:
            if amino_acid in opt_codons_E_coli and codon in opt_codons_E_coli[amino_acid]:
                tRNA_map[codon] = ["pref"]
            else:
                tRNA_map[codon] = ["nonpref"]

    return tRNA_map

def randomize_codons(record, feature, fop):
    # Extract nucleotide sequence
    nuc_seq = feature.location.extract(record).seq

    # Split the sequence into codons
    codons = [str(nuc_seq[i:i + 3]) for i in range(0, len(nuc_seq), 3)]

    # Categorize codons into optimal and non-optimal categories
    optimal_codons = []
    nonoptimal_codons = []

    for codon in codons:
        found_optimal = False
        for aa, opt_list in OPT_CODONS_E_COLI.items():
            if codon in opt_list:
                optimal_codons.append(codon)
                found_optimal = True
                break
        if not found_optimal:
            nonoptimal_codons.append(codon)

    # Calculate the number of codons to sample from each category based on the given percentages
    num_optimal = int(len(codons) * fop)
    num_nonoptimal = len(codons) - num_optimal  # The rest will be non-optimal

    # Ensure the number of sampled codons does not exceed the available codons and repeat if necessary
    randomized_optimal = (random.sample(optimal_codons, len(optimal_codons)) * (num_optimal // len(optimal_codons))) + random.sample(optimal_codons, num_optimal % len(optimal_codons))
    randomized_nonoptimal = (random.sample(nonoptimal_codons, len(nonoptimal_codons)) * (num_nonoptimal // len(nonoptimal_codons))) + random.sample(nonoptimal_codons, num_nonoptimal % len(nonoptimal_codons))

    # Combine and shuffle the selected codons
    randomized_codons = randomized_optimal + randomized_nonoptimal
    random.shuffle(randomized_codons)

    # Reconstruct the nucleotide sequence from the randomized codons
    randomized_seq = "".join(randomized_codons)

    return randomized_seq
    
def main(fop, charge_rate, pref_proportion, seed_val, ribo_speed, trna_count):
    sim = pt.Model(cell_volume=CELL_VOLUME)
    Entrez.email = "kelly.to@utexas.edu"

    # Download T7 wild-type genbank records
    for attempt in range(20):
        try:
            handle = Entrez.efetch(db="nuccore", id=["NC_001604"], rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            break
        except urllib.error.HTTPError as e:
            if e.code == 429:
                time.sleep(2 ** attempt)  # Exponential backoff
                print(f"Rate limit exceeded. Retrying...")
            else:
                print(f"Failed to fetch GenBank record: {e}")
                break

    # Extract the sequence of gene 10A
    # gene_10A_sequence = None
    gene10_start = None
    gene10_stop = None
    gene10_feature = None
    for feature in record.features:
        if feature.type == "gene" and "gene 10A" in feature.qualifiers["note"]:
            # print(f"Found gene 10A: {feature}")
            # gene_10A_sequence = feature.extract(record.seq)
            gene10_start = feature.location.start
            gene10_stop = feature.location.end
            gene10_feature = feature
            break

    # Use randomize_codons to generate a randomized gene 10A sequence and replace the original gene 10A sequence
    randomized_sequence = randomize_codons(record, gene10_feature, fop)

    # Replace the original gene 10A sequence with the randomized sequence
    record.seq = record.seq[:gene10_start] + randomized_sequence + record.seq[gene10_stop:]
    # print(f"Randomized gene 10A sequence with {fop*100}% optimal codons")
    # print(record.seq[gene10_start:gene10_stop])

    phage = pt.Genome(name="phage", length=len(record.seq))
    phage.add_sequence(str(record.seq))

#----------------------------------------------------------------------------

    for feature in record.features:
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

    mask_interactions = ["rnapol-1", "rnapol-3.5",
                         "ecolipol", "ecolipol-p", "ecolipol-2", "ecolipol-2-p"]
    phage.add_mask(500, mask_interactions)

    sim.register_genome(phage)

    sim.add_polymerase("rnapol-1", 35, 230, 0)
    sim.add_polymerase("rnapol-3.5", 35, 230, 0)
    sim.add_polymerase("ecolipol", 35, 45, 0)
    sim.add_polymerase("ecolipol-p", 35, 45, 0)
    sim.add_polymerase("ecolipol-2", 35, 45, 0)
    sim.add_polymerase("ecolipol-2-p", 35, 45, 0)

    # footprint, speed, copy number
    sim.add_ribosome(30, ribo_speed, 0)   # originally speed of 30

    sim.add_species("bound_ribosome", 10000)   # originally 10000

    sim.add_species("bound_ecolipol", 1800)   # originally 1800
    sim.add_species("bound_ecolipol_p", 0)
    sim.add_species("ecoli_genome", 0)
    sim.add_species("ecoli_transcript", 0)

    # translation reactions
    sim.add_reaction(1e6, ["ecoli_transcript", "__ribosome"], [
                     "bound_ribosome"])

    sim.add_reaction(0.04, ["bound_ribosome"], [
                     "__ribosome", "ecoli_transcript"])

    sim.add_reaction(0.001925, ["ecoli_transcript"], ["degraded_transcript"])

    # transcription reactions
    sim.add_reaction(1e7, ["ecolipol", "ecoli_genome"], ["bound_ecolipol"])

    sim.add_reaction(
        0.3e7, ["ecolipol-p", "ecoli_genome"], ["bound_ecolipol_p"])

    sim.add_reaction(0.04, ["bound_ecolipol"], [
                     "ecolipol", "ecoli_genome", "ecoli_transcript"])

    sim.add_reaction(0.04, ["bound_ecolipol_p"], [
                     "ecolipol-p", "ecoli_genome", "ecoli_transcript"])

    # # protein modification reactions
    sim.add_reaction(3.8e7, ["protein_kinase-0.7", "ecolipol"],
                     ["ecolipol-p", "protein_kinase-0.7"])

    sim.add_reaction(3.8e7, ["protein_kinase-0.7", "ecolipol-2"],
                     ["ecolipol-2-p", "protein_kinase-0.7"])

    sim.add_reaction(3.8e7, ["gp-2", "ecolipol"], ["ecolipol-2"])

    sim.add_reaction(3.8e7, ["gp-2", "ecolipol-p"], ["ecolipol-2-p"])

    sim.add_reaction(1.1, ["ecolipol-2-p"], ["gp-2", "ecolipol-p"])

    sim.add_reaction(1.1, ["ecolipol-2"], ["gp-2", "ecolipol"])

    # # RNA polymerase and lysozyme reactions
    sim.add_reaction(3.8e9, ["lysozyme-3.5", "rnapol-1"], ["rnapol-3.5"])

    sim.add_reaction(3.5, ["rnapol-3.5"], ["lysozyme-3.5", "rnapol-1"])

#--------------------------------------------------------------------------

    # strength of tRNA re-charging reaction [tRNA_a, tRNA_b]   
    TRNA_CHRG_RATES = [charge_rate, charge_rate]   # originally [100.0, 100.0]

    # tRNA proportions, i.e. 10% total tRNA is pref, other 90% is nonpref 
    nonpref_proportion = 1 - pref_proportion
    TRNA_PROPORTIONS = (pref_proportion, nonpref_proportion)   # originally (0.1, 0.9)

    TOTAL_TRNA = trna_count

    # tRNA/codon mapping: 
    tRNA_map = tRNA_map_maker()

    # initial tRNA species counts:
    # pref: 250 (total tRNA * 0.1) charged, 0 uncharged
    # nonpref: 2250 (total tRNA * 0.9) charged, 0 uncharged
    tRNA_counts = {"pref": [int(TOTAL_TRNA*TRNA_PROPORTIONS[0]), 0], "nonpref": [int(TOTAL_TRNA*TRNA_PROPORTIONS[1]), 0]}
    # tRNA charging rates: pref: 100.0, nonpref: 100.0
    tRNA_rates = {"pref": TRNA_CHRG_RATES[0], "nonpref": TRNA_CHRG_RATES[1]}
    sim.add_trna(tRNA_map, tRNA_counts, tRNA_rates)

    sim.seed(seed_val)

    # generate a unique filename based on the charge
    # output_dir = f"/scratch/10081/kellyktvt/trna_parallel_output/revised_dynamic_pref{pref_proportion}_charge{charge_rate}_ribospeed{ribo_speed}"
    output_dir = f"data/simulation/phage/revised_dynamic_pref{pref_proportion}_charge{charge_rate}_ribospeed{ribo_speed}"
    output_filename = os.path.join(output_dir, f"trna_phage_pref{pref_proportion}_{seed_val}_fop{fop}_ribospeed{ribo_speed}_trna{trna_count}.tsv")

    sim.simulate(time_limit=1200, time_step=100, 
                 output=output_filename)

#--------------------------------------------------------------------------

if __name__ == "__main__":
    # Take in arguments from the jobs file
    # fraction of optimal codons
    fop = float(sys.argv[1])
    # tRNA charging rate
    charge_rate = float(sys.argv[2])
    # proportion of preferred tRNAs
    pref_proportion = float(sys.argv[3])
    # seed value
    seed_val = int(sys.argv[4])
    # ribosome speed
    ribo_speed = float(sys.argv[5])
    # trna_count
    trna_count = int(sys.argv[6])

    main(fop, charge_rate, pref_proportion, seed_val, ribo_speed, trna_count)