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


# Function to calculate Length of 5' UTR
def calculate_five_prime_utr(sequence):
    start_codon = 'AUG'

    if start_codon in sequence:
        return sequence.index(start_codon)
    else:
        return 0

def calculate_kozak_pos_1(sequence):
    kozak_start = 50 - 6
    encoding = {"A": 1, "U": 2, "G": 3, "C": 4}
    return encoding.get(sequence[kozak_start], 0)

def calculate_kozak_pos_4(sequence):
    kozak_start = 50 + 3
    encoding = {"A": 1, "U": 2, "G": 3, "C": 4}
    return encoding.get(sequence[kozak_start], 0)

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

def evaluate_model(model, X_test):
    y_pred = model.predict(X_test)
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

    # Calculate features button
    if st.button("Calculate Features"):
        if sequence or uploaded_file:
            if sequence:
                X = calculate_features(sequence)
            else:
                content = uploaded_file.read().decode("utf-8")
                sequences = [line for line in content.split("\n") if not line.startswith('>')]
                dfs = []
                for seq in sequences:
                    df = calculate_features(seq)
                    if df is not None:
                        dfs.append(df)
                if dfs:
                    X = pd.concat(dfs, ignore_index=True)
                    st.write("Calculated Features:")
                    st.write(X)  # Print the calculated features dataset X
                else:
                    st.write("No valid sequences found in the uploaded file.")
                    return

    # Start Prediction Button
    if st.button("Start Prediction"):
        if sequence:
            X = calculate_features(sequence)
            if X is None:
                st.write("Invalid sequence.")
                return
        else:
            st.write("No sequence available")
            return

        # Load Models
        rf_model_path = "tir_rf_model.pkl"

        with open(rf_model_path, 'rb') as f:
            rf_model = pickle.load(f)

        # Evaluate Random Forest Model
        rf_y_pred = evaluate_model(rf_model, X)

        # Create a DataFrame with predictions
        df_predictions = pd.DataFrame({
            'Gene Sequence': sequence,
            'Random Forest Predictions': rf_y_pred
        })

        # Display the predictions
        st.write("Prediction Results:")
        st.write(df_predictions)

if __name__ == "__main__":
    main()
