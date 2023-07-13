import pandas as pd
import RNA


def features(gene_sequence, start_codon, stop_codons):

    start_codon = gene_sequence.find("AUG")
    if start_codon == -1:
        raise ValueError("No start codon found in gene sequence")

    stop_codons = ("UAA", "UAG", "UGA")
    stop_codon = -1

    for stop_codon in stop_codons:
        index = gene_sequence.find(stop_codon, start_codon + 3)
        if index != -1:
            stop_codon = index
            break

    if stop_codon == -1:
        raise ValueError("No stop codon found in gene sequence")

    cds_sequence = gene_sequence[start_codon:stop_codon + 3]

    folding_energy_80 = calculate_folding_energy_80(gene_sequence,start_codon)
    folding_energy_70 = calculate_folding_energy_70(gene_sequence)
    length_of_5prime_utr = len(gene_sequence[:start_codon])
    koz_score = kozak_score(gene_sequence)
    #in_frame_AUG = in_frame_AUG(cds_sequence)
    in_frame_AUG = 0

    encoding = {"A": 1, "U": 2, "G": 3, "C": 4}

   # N1 = encoding.get(cds_sequence[start_codon - 6])
   # N4 = encoding.get(cds_sequence[start_codon - 3])
   
    N1 = 2
    N4 = 4
    #koz_score = 4
    
   

    gene_features = {
        "gene_length": len(cds_sequence),
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


def calculate_folding_energy_80(gene_sequence, start_codon):
    sequence_80 = gene_sequence[start_codon - 40:start_codon + 40]
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


def in_frame_AUG(cds_sequence):
    num_in_frame_aug = 0
    for i in range(0, len(cds_sequence), 3):
        if cds_sequence[i:i + 3] == "AUG":
            num_in_frame_aug += 1
    return num_in_frame_aug

def main():

    gene_sequence = "CACCAGGUUUUUGGCUUUUUAGAUUUUAUCCCCUUCCAGCAUGAGGAUUGGCACGGAUGCUAACGUGAUAAUCUUGGCUGUAG"
    start_codon = 41
    stop_codons = 81

    gene_features = features(gene_sequence, start_codon, stop_codons)
    print(gene_features)



if __name__ == "__main__":
    main()

