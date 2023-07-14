import streamlit as st
import numpy as np
import pandas as pd
import RNA
import pickle
import features
import base64
import warnings
warnings.filterwarnings("ignore")

def load_data(file):
    df = pd.read_csv(file)
    return df

def evaluate_model(model, X_test):
    y_pred = model.predict(X_test)
    return y_pred

def main():
    gene_sequence = st.text_input("Enter the mRNA sequence:")
    start_codon = st.text_input("Enter the start codon position:")
    stop_codons = st.text_input("Enter the stop codons:")

    if st.button("Calculate features"):
        start_codon_index = int(start_codon)

        gene_features = features.features(gene_sequence, int(start_codon), int(stop_codons))
        data = pd.DataFrame(gene_features, index=["Value"])

        st.subheader("Calculated Features:")
        st.write(data.to_html(index=False), unsafe_allow_html=True)

        # Download link for the DataFrame
        csv = data.to_csv(index=False)
        b64 = base64.b64encode(csv.encode()).decode()
        href = f'<a href="data:file/csv;base64,{b64}" download="gene_features.csv">Download here</a>'
        st.markdown(href, unsafe_allow_html=True)

    # File Upload
    file = st.file_uploader("Upload the file", type=["csv","xlsx"], key='file_uploader')
    if not file:
        st.warning("Please upload a file in excel or csv format")
        return


    # Load Data
    df = load_data(file)
    st.write("Data:")
    st.write(df)
    
    

    # Start Prediction Button
    if st.button("Predict translation initiation rate"):  
        # Load Models
        rf_model_path = "tir_rf_model.pkl"

        with open(rf_model_path, 'rb') as f:
            rf_model = pickle.load(f)

        # Evaluate Random Forest Model
        rf_y_pred = evaluate_model(rf_model, df)
        #st.write(rf_y_pred)

        # Create a DataFrame with predictions
        df_predictions = pd.DataFrame({
            'Random Forest Predictions': rf_y_pred,
        })

        # Provide a download link for predictions
        csv = df_predictions.to_csv(index=False)
        st.write(csv)

if __name__ == "__main__":
    main()
