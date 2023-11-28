import numpy as np
import pandas as pd
import streamlit as st
import random
import features,InitiationRate
import OptimizationUTRcodon, OptimizationUTR
import warnings
warnings.filterwarnings("ignore")

# Page title
st.markdown("""
# TIR Predictor 

This allows user to predict Translation Initation Rate in Saccharomyces cerevisiae using mRNA features using Machine Learning methods.

### Introduction: 

Translation initiation, which is the rate-limiting step in protein synthesis, can vary significantly and have a profound impact on cellular protein levels. 
Multiple molecular factors, such as mRNA structure stability, coding sequence length, and specific motifs in mRNA, influence the translation initiation rate, 
allowing precise control of protein synthesis. Despite the crucial role of translation initiation rate, accurately predicting its absolute values based on mRNA 
sequence features remains challenging. To address this issue, we developed a machine learning model specifically trained to predict the in vivo 
initiation rate in S. cerevisiae transcripts. 

Further using this app the user can optimize the gene and achieve their desired target initiation rate using 2 methods:
1. Optimization with UTR
2. Optimization with UTR and codon

This has been developed on python 3.9 

### How to use:

1. Enter the CDS sequence, start codon and stop codon sites. 
2. Click on the "Calculate Features & Predict Initiation rate" to initiate the analysis.
3. If you want to optimize enter the target inititation rate.
4. Enter how many time you want to iterate the process.
5. Wait for the message "Optimization completed!"


Note: The output file will contain the perdcited translation initiation rate of the input given for specific given genes.It works properly with one or more genes.

**Credits**
- Built in `Python` + `Streamlit` by Sulagno Chakraborty, Inayat Ullah Irshad, Mahima and Ajeet K. Sharma
[[Read the Paper]]().
---
""")

# Streamlit web app
def main():
    st.title("Initiation rate Prediction")

    # Input fields for gene sequence, start codon index, stop codon index, and target initiation rate
    gene_sequence = st.text_input("Enter the mRNA sequence:")
    start_codon_index = st.number_input("Start Codon Index:", min_value=1, value=18, step=1)
    stop_codon_index = st.number_input("Stop Codon Index:", min_value=1, value=732, step=1)

    if st.button("Calculate Features & Predict Initiation rate"):
        gene_features = features.features(gene_sequence, int(start_codon_index), int(stop_codon_index))
        df = pd.DataFrame([gene_features], columns=["gene_length", "folding_energy_70", "folding_energy_80", "length_of_5prime_utr", "kozak_score", "N1", "N4", "in_frame AUG"])
        st.write("The dataset of calculated features are as follows:")
        st.dataframe(df)
        rate = InitiationRate.InitiationRate(gene_sequence, int(start_codon_index), int(stop_codon_index))
        st.write("The Translation Initiation rate predicted is:", rate)
    
    st.title("Gene Optimization")
    
    target_initiation_rate = st.number_input("Target Initiation Rate:", min_value=0.0, value=0.8, step=0.01)
    iterations = st.number_input("Enter the number of iterations:", min_value=1, value=100, step=1)

    method_options = ["Optimization with UTR", "Optimization with UTR and codon"]
    selected_method = st.selectbox("Choose a method", method_options)
    

    if st.button("Optimize"):
        if selected_method == "Optimization with UTR":
            opt = OptimizationUTRcodon.OptimizationUTRcodon(gene_sequence, start_codon_index, stop_codon_index, target_initiation_rate, iterations)
        elif selected_method == "Optimization with UTR and codon":
            opt = OptimizationUTR.OptimizationUTR(gene_sequence, start_codon_index, stop_codon_index, target_initiation_rate, iterations)
            
            
        print(opt)
        st.success("Optimization completed!")
    
    

if __name__ == "__main__":
    main()
