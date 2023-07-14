import streamlit as st
import numpy as np
import pandas as pd
import RNA
import pickle
import features
import base64

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
        df = pd.DataFrame(gene_features, index=["Value"])

        st.subheader("Calculated Features:")
        st.write(df.to_html(index=False, justify="center"), unsafe_allow_html=True)

        # Download link for the DataFrame
        csv = df.to_csv(index=False)
        b64 = base64.b64encode(csv.encode()).decode()
        href = f'<a href="data:file/csv;base64,{b64}" download="gene_features.csv">Download here</a>'
        st.markdown(href, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
