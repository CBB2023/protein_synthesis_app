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

# Function to calculate Kozak Score
def calculate_kozak_score(sequence):
    if len(sequence) < 9:
        return 0

    koz = sequence[50-6:50] + sequence[50+3:50+6]
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
    score += 1

    return score

# Function to calculate folding energy of first 70 base pairs
def calculate_folding_energy_70(sequence):
    sequence_70 = sequence[:70]
    (ss, mfe) = RNA.fold(sequence_70)
    return "{:.2f}".format(mfe)

# Function to calculate folding energy of 40 base pairs left of "AUG" plus 40 base pairs of "AUG"
def calculate_folding_energy_80(sequence):
    aug_index = sequence.find("AUG")
    sequence_80 = sequence[aug_index - 40:aug_index + 43]
    (ss, mfe) = RNA.fold(sequence_80)
    return "{:.2f}".format(mfe)

def evaluate_model(model, X):
    y_pred = model.predict(X)
    return y_pred    

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
                sequences = [line for line in content.split("\n") if not line.startswith('>')]
                df = pd.DataFrame({'Sequence': sequences})

            # Exclude empty sequences
            df = df[df['Sequence'] != '']

            if not df.empty:
                # Calculate Gene Length
                df['Gene Length'] = df['Sequence'].str.len()

                # Calculate Length of 5' UTR
                start_codon = 'AUG'
                df['Length of 5\' UTR'] = df['Sequence'].apply(lambda seq: seq.index(start_codon) if start_codon in seq else 0)

                # Calculate Kozak Score
                df['Kozak Score'] = df['Sequence'].apply(calculate_kozak_score)

                # Calculate folding energy of first 70 base pairs
                df['Folding Energy 70'] = df['Sequence'].apply(calculate_folding_energy_70)

                # Calculate folding energy of 40 base pairs left of "AUG" plus 40 base pairs of "AUG"
                df['Folding Energy 80'] = df['Sequence'].apply(calculate_folding_energy_80)

                # Define a dictionary that maps each letter to its corresponding value
                encoding = {"A": 1, "U": 2, "G": 3, "C": 4}

                # Use the map() function to apply the encoding to the first and fourth letters of each string
                df["Kozak pos. 1"] = df["Sequence"].str[50-6].map(encoding)
                df["Kozak pos. 4"] = df["Sequence"].str[50+3].map(encoding)

                X = df[['Gene Length', 'Length of 5\' UTR', 'Kozak Score', 'Kozak pos. 1', 'Kozak pos. 4', 'Folding Energy 70', 'Folding Energy 80']]
                # Download dataset
                csv = X.to_csv(index=False)
                b64 = base64.b64encode(csv.encode()).decode()
                href = f'<a href="data:file/csv;base64,{b64}" download="dataset.csv">Download Dataset</a>'
                st.markdown(href, unsafe_allow_html=True)

                # Start Prediction Button
                if st.button("Start Prediction"):
                    # Load Model
                    rf_model_path = "tir_rf_model.pkl"

                    with open(rf_model_path, 'rb') as f:
                        rf_model = pickle.load(f)

                    # Evaluate Random Forest Model
                    rf_y_pred = evaluate_model(rf_model, X)
                    
                    # Create a DataFrame with predictions
                    X_predictions = pd.DataFrame({
                        'Gene Sequence': sequence,
                        'Random Forest Predictions': rf_y_pred 
                    })

                    # Provide a download link for predictions
                    csv = X_predictions.to_csv(index=False)
                    b64 = base64.b64encode(csv.encode()).decode()  
                    # Convert DataFrame to base64 encoding
                    href = f'<a href="data:file/csv;base64,{b64}" download="predictions.csv">Download Predictions</a>'
                    st.markdown("Download Predictions:")
                    st.markdown(href, unsafe_allow_html=True)


if __name__ == "__main__":
    main()
