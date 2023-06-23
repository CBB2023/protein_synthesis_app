import streamlit as st
import pandas as pd

# Page title
st.markdown("""
# Translation Initation Rate Prediction App 

This app allows you to predict Translation Initation Rate in Saccharomyces cerevisiae using Machine Learning methods

**Credits**
- App built in `Python` + `Streamlit` by Sulagno Chakraborty, Inayat Ullah Irshad and Dr. Ajeet K. Sharma
[[Read the Paper]]().
---
""")

import streamlit as st
import pandas as pd
import joblib

# Load the TIR prediction model
model = joblib.load("tir_rf_model.pkl")

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

# Streamlit app
def main():
    st.title("Gene Sequence Analysis")
    
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
        df['Kozak Score'] = df['Sequence'].apply(calculate_kozak_score)

        # Define a dictionary that maps each letter to its corresponding value
        encoding = {"A": 1, "U": 2, "G": 3, "C": 4}

        # Use the map() function to apply the encoding to the first and fourth letters of each string
        df["First Letter"] = df["Sequence"].str[50-6].map(encoding)
        df["Fourth Letter"] = df["Sequence"].str[50+3].map(encoding)

        # Predict TIR using the model
        df['TIR Prediction'] = model.predict(df[['Gene Length', "Length of 5' UTR", 'Kozak Score', 'First Letter', 'Fourth Letter']])

        # Save the predicted TIRs to a CSV file
        df.to_csv("predicted_tirs.csv", index=False)

        # Print the dataset
        st.write(df)
        
        # Provide a download link for the CSV file
        st.markdown("[Download Predicted TIRs CSV](predicted_tirs.csv)")
    else:
        st.warning("No sequences found.")

if __name__ == "__main__":
    main()
