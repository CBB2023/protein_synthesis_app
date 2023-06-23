import streamlit as st
import pandas as pd
import pickle

# Page title
st.markdown("""
# Translation Initiation Rate Prediction App

This app allows you to predict Translation Initiation Rate in Saccharomyces cerevisiae using Machine Learning methods.

**Credits**
- App built in `Python` + `Streamlit` by Sulagno Chakraborty, Inayat Ullah Irshad, Mahima and Dr. Ajeet K. Sharma
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
    # Title of the dialogue box
    st.subheader("Enter a gene sequence")

    # User input - Text area for entering the sequence
    sequence = st.text_area("Sequence")

    # Title for uploading file section
    st.subheader("Or, upload a file")

    # File upload
    uploaded_file = st.file_uploader("Upload .txt file", type="txt")

    # Display selected file name or "No file selected"
    if uploaded_file is not None:
        file_name = uploaded_file.name
        st.write("Selected file:", file_name)
    else:
        file_name = "No file selected"

    # Calculate features and generate dataset button
    if st.button("Calculate Features and Generate Dataset"):
        if sequence or uploaded_file:
            # Create DataFrame
            if sequence:
                df = pd.DataFrame({'Sequence': [sequence]})
            else:
                content = uploaded_file.read().decode("utf-8")
                sequences = content.split("\n")
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
            df['Initiation Rate'] = rf_model.predict(X)

            # Download dataset
            csv = df.to_csv(index=False)
            b64 = base64.b64encode(csv.encode()).decode()
            href = f'<a href="data:file/csv;base64,{b64}" download="dataset.csv">Download Dataset</a>'
            st.markdown(href, unsafe_allow_html=True)

    # Prediction button
    if st.button("PREDICT"):
        if sequence or uploaded_file:
            # Create DataFrame
            if sequence:
                df = pd.DataFrame({'Sequence': [sequence]})
            else:
                content = uploaded_file.read().decode("utf-8")
                sequences = content.split("\n")
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
            df['Initiation Rate'] = rf_model.predict(X)

            # Download predictions
            csv = df.to_csv(index=False)
            b64 = base64.b64encode(csv.encode()).decode()
            href = f'<a href="data:file/csv;base64,{b64}" download="predictions.csv">Download Predictions</a>'
            st.markdown(href, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
