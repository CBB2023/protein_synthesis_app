import streamlit as st
import pandas as pd
import pickle

# Page title
st.markdown("""
# Translation Initiation Rate Prediction App

This app allows you to predict Translation Initiation Rate in Saccharomyces cerevisiae using Machine Learning methods.

**Credits**
- App built in `Python` + `Streamlit` by Sulagno Chakraborty, Inayat Ullah Irshad, and Dr. Ajeet K. Sharma
[[Read the Paper]]().
---
""")

# Function to calculate Kozak Score
def calculate_kozak_score(sequence):
    koz = sequence[50-6:50] + sequence[50+3:50+6]
    score = 0
    
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
    score += 1
    
    return score

# Load Models
rf_model_path = "tir_rf_model.pkl"

with open(rf_model_path, 'rb') as f:
    rf_model = pickle.load(f)

# Streamlit app
def main():
    
    # User input
    option = st.radio("Select option:", ("Enter sequence", "Upload .txt file"))
    
    if option == "Enter sequence":
        sequence = st.text_area("Enter the sequence:")
        sequences = [sequence]
    else:
        uploaded_file = st.file_uploader("Upload .txt file", type="txt")
        if uploaded_file is not None:
            content = uploaded_file.read().decode("utf-8")
            sequences = content.split("\n")
    
    if sequences:
        # Create DataFrame
        df = pd.DataFrame({'Sequence': sequences})

        # Calculate Gene Length
        df['Gene Length'] = df['Sequence'].str.len()

        # Calculate Length of 5' UTR
        start_codon = 'AUG'
        df['Length of 5\' UTR'] = df['Sequence'].apply(lambda seq: seq.index(start_codon) if start_codon in seq else 0)

        # Calculate Kozak Score
        kozak_scores = []
        for sequence in df['Sequence']:
            kozak_score = calculate_kozak_score(sequence)
            kozak_scores.append(kozak_score)
        df['Kozak Score'] = kozak_scores

        # Define a dictionary that maps each letter to its corresponding value
        encoding = {"A": 1, "U": 2, "G": 3, "C": 4}

        # Use the map() function to apply the encoding to the first and fourth letters of each string
        df["First Letter"] = df["Sequence"].str[50-6].map(encoding)
        df["Fourth Letter"] = df["Sequence"].str[50+3].map(encoding)

        # Perform predictions
        X = df[['Gene Length', 'Length of 5\' UTR', 'Kozak Score', 'First Letter', 'Fourth Letter']]
        df['Translation Initiation Rate (RF)'] = rf_model.predict(X)

        # Save predictions to CSV
        df.to_csv('predictions.csv', index=False)

        # Display predictions
        st.write(df)
    else:
        st.warning("No sequences found.")

if __name__ == "__main__":
    main()
