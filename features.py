import pandas as pd
import RNA

def features(gene_sequence, start_codon_index, stop_codon_index):
    start_codon = "AUG"
    stop_codon = ["UAA", "UAG", "UGA"]

    if gene_sequence[start_codon_index:start_codon_index + 3] == start_codon and \
            gene_sequence[stop_codon_index:stop_codon_index + 3] in stop_codon:
        cds_sequence = gene_sequence[start_codon_index :stop_codon_index + 3]
        cds_sequence = cds_sequence[3:-3]
    else:
        return None
    gene_length = len(cds_sequence)
    folding_energy_80 = calculate_folding_energy_80(gene_sequence,start_codon_index)
    folding_energy_70 = calculate_folding_energy_70(gene_sequence)
    length_of_5prime_utr = len(gene_sequence[:start_codon_index])
    koz_score = kozak_score(gene_sequence)
    in_frame_AUG = calculate_in_frame_AUG(cds_sequence)

    encoding = {"A": 1, "U": 2, "G": 3, "C": 4}

    N1 = encoding.get(gene_sequence[start_codon_index - 5])
    N4 = encoding.get(gene_sequence[start_codon_index - 3])
    
    gene_features = {
    "gene_length": gene_length,
    "folding_energy_80": folding_energy_80,
    "folding_energy_70": folding_energy_70,
    "length of 5prime utr": length_of_5prime_utr,
    "kozak score": koz_score,
    "N1": N1,
    "N4": N4,
    "in_frame AUG" : in_frame_AUG
    }

    return gene_features


def calculate_folding_energy_70(gene_sequence):
    sequence_70 = gene_sequence[:70]
    (ss, mfe) = RNA.fold(sequence_70)
    return "{:.2f}".format(mfe)


def calculate_folding_energy_80(gene_sequence, start_codon_index):
    sequence_80 = gene_sequence[start_codon_index - 40:start_codon_index + 40]
    (ss, mfe) = RNA.fold(sequence_80)
    return "{:.2f}".format(mfe)


def kozak_score(gene_sequence):
    koz = gene_sequence[50-6:50] + gene_sequence[50+3:50+6]
    score = 0

    if len(koz) < 9:
          return 0

    if koz[0] == "A" or koz[0] == "U":
          score += 1
    if koz[1] == "A":
          score += 1
    if koz[2] == "A" or koz[2] == "C":
          score += 1
    if koz[3] == "A":
          score += 1
    if koz[4] == "A" or koz[2] == "C":
          score += 1
    if koz[5] == "A":
          score += 1
    if koz[6] == "U":
          score += 1
    if koz[7] == "C":
          score += 1
    if koz[8] == "U" or koz[8] == "C":
          score += 1

    return score


def calculate_in_frame_AUG(cds_sequence):
    num_in_frame_aug = 0
    for i in range(0, len(cds_sequence), 3):
        if cds_sequence[i:i + 3] == "AUG":
            num_in_frame_aug += 1
    return num_in_frame_aug

if __name__ == "__main__":
    gene_sequence = "CACCAGGUUUUUGGCUUUUUAGAUUUUAUCCCCUUCCAGCAUGAGGAUUGGCACGGAUGCUAACGUGAUAAUCUUGGCUGUAG"
    start_codon_index = 40
    stop_codon_index = 80
    gene_features = features(gene_sequence, start_codon_index, stop_codon_index)
    print(gene_features)
