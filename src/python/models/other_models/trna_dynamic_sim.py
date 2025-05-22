import pinetree as pt


## parameters
RB_COPY = 500 # number of ribosomes in simulation
TS_COPY = 100 # number of transcripts (mRNAs) in simulation
RBS_STRENGTH = 100000.0 # strength of ribosome binding to an mRNA
TRNA_CHRG_RATES = [100.0, 100.0] # strength of tRNA re-charging reaction [tRNA_a, tRNA_b]
TRNA_PROPORTIONS = (0.1, 0.9) # tRNA proportions, i.e. 90% total tRNA is type A, other 10% is type B
TOTAL_TRNA = 2500 # total tRNA
SPEED = 0.5 # rate constant describing how efficient ribosome movement is
TIME_LIMIT = 100
TIME_STEP = 5
SEED = 1


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


def trna_simulation(dir):

    # 300 nt coding region + 50 nt buffer (30 on the left, 20 on the right)
    # transcripts coding region is 10% codon "AAA" and 90% codon "TAT"
    seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAATATTATTATTATAAATATTATTATTATTATTATTATAAATATTATTATTATTATTATTATTATTATAAATATAAATATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATAAATATTATTATTATTATAAATATTATTATTATTATTATTATTATTATTATTATTATTATTATTATTATAAATATTATTATTATTATTATTATTATTATAAATATTATTATAAAAAATATTATTATTATTATTATTATTATAAAAAAAAAAAAAAAAAAAA"

    sim = pt.Model(cell_volume=8e-16)
    sim.seed(SEED)
    sim.add_ribosome(copy_number=RB_COPY, speed=SPEED, footprint=15)

    i = 0
    while i < TS_COPY: # add 100 copies of transcript to the simulation (each copy is identical)
        transcript = pt.Transcript("transcript", 350)
        transcript.add_gene(name="proteinX", start=31, stop=330,
                     rbs_start=(31 - 15), rbs_stop=31, rbs_strength=RBS_STRENGTH)
        transcript.add_seq(seq=seq)
        sim.register_transcript(transcript)
        i += 1
    
    #--------------------------------------------------------------------------

    TRNA_CHRG_RATES = [100.0, 100.0] # strength of tRNA re-charging reaction [tRNA_a, tRNA_b]
    TRNA_PROPORTIONS = (0.1, 0.9) # tRNA proportions, i.e. 90% total tRNA is type A, other 10% is type B
    TOTAL_TRNA = 2500 # total tRNA

    ## tRNA/codon mapping: 
    tRNA_map = tRNA_map_maker()

    # initial tRNA species counts:
    tRNA_counts = {"pref": [int(TOTAL_TRNA*TRNA_PROPORTIONS[0]), 0], "nonpref": [int(TOTAL_TRNA*TRNA_PROPORTIONS[1]), 0]}
    # tRNA charging rates: pref: 100.0, nonpref: 100.0
    tRNA_rates = {"pref": TRNA_CHRG_RATES[0], "nonpref": TRNA_CHRG_RATES[1]}
    sim.add_trna(tRNA_map, tRNA_counts, tRNA_rates)

    #--------------------------------------------------------------------------
    
    output = "data/simulation/phage/trna_dynamic_sim.tsv"
    sim.simulate(time_limit=TIME_LIMIT, time_step=TIME_STEP, output=f"{dir}/{output}")

    
if __name__ == "__main__":
    trna_simulation(".")