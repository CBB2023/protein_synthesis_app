import streamlit as st
import pickle
import pandas as pd
import RNA
import base64

# Page title
st.markdown("""
# Translation Initiation Rate Prediction App

This app allows you to predict Translation Initiation Rate in Saccharomyces cerevisiae using Machine Learning methods.

**Credits**
- App built in `Python` + `Streamlit` by Sulagno Chakraborty, Inayat Ullah Irshad, Mahima, and Ajeet K. Sharma
[[Read the Paper]]().
---
""")

# Function to calculate CDS length
def calculate_cds(sequence):
    start_codon = 'AUG'
    stop_codons = ['UAA', 'UAG', 'UGA']

    cds_start = sequence.find(start_codon)
    cds_end = -1

    for i in range(cds_start + len(start_codon), len(sequence), 3):
        codon = sequence[i:i + 3]
        if codon in stop_codons:
            cds_end = i + len(codon)
            break

    if cds_start != -1 and cds_end != -1:
        cds_length = cds_end - cds_start
    else:
        cds_length = 0

    return cds_length

# Function to calculate Length of 5' UTR
def calculate_five_prime_utr(sequence):
    start_codon = 'AUG'

    if start_codon in sequence:
        return sequence.index(start_codon)
    else:
        return 0

# Function to calculate 1st position of each Kozak sequence
def calculate_kozak_pos_1(sequence):
    kozak_start = 50 - 6
    return sequence[kozak_start]

# Function to calculate 4th position of each Kozak sequence
def calculate_kozak_pos_4(sequence):
    kozak_start = 50 + 3
    return sequence[kozak_start]

# Function to calculate features
def calculate_features(sequence):
    # Create DataFrame
    df = pd.DataFrame({'Sequence': [sequence]})
    
    # Exclude empty sequences
    df = df[df['Sequence'] != '']

    if not df.empty:
        # Calculate Gene Length
        df['Gene Length'] = df['Sequence'].str.len()

        # Calculate Length of 5' UTR
        df['Length of 5\' UTR'] = df['Sequence'].apply(calculate_five_prime_utr)

        # Calculate CDS length
        df['CDS Length'] = df['Sequence'].apply(calculate_cds)

        # Calculate Kozak pos. 1
        df['Kozak pos. 1'] = df['Sequence'].apply(calculate_kozak_pos_1)

        # Calculate Kozak pos. 4
        df['Kozak pos. 4'] = df['Sequence'].apply(calculate_kozak_pos_4)

        # Calculate folding energy of first 70 base pairs
        df['Folding Energy 70'] = df['Sequence'].apply(calculate_folding_energy_70)

        # Calculate folding energy of 40 base pairs left of "AUG" plus 40 base pairs of "AUG"
        df['Folding Energy 80'] = df['Sequence'].apply(calculate_folding_energy_80)

        X = df[['Gene Length', 'Length of 5\' UTR', 'CDS Length', 'Kozak pos. 1', 'Kozak pos. 4', 'Folding Energy 70', 'Folding Energy 80']]
        return X
    else:
        return None

# Rest of the code...

