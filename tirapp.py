import streamlit as st
import pickle
import pandas as pd
import RNA


def calculate_five_prime_utr(sequence, start_codon):
    start_codon = str(start_codon)
    if start_codon in sequence:
        return sequence.index(start_codon)
    else:
        return 0


def calculate_kozak_pos_1(sequence, kozak_start):
    encoding = {"A": 1, "U": 2, "G": 3, "C": 4}
    return encoding.get(sequence[kozak_start], 0)


def calculate_kozak_pos_4(sequence, kozak_start):
    encoding = {"A": 1, "U": 2, "G": 3, "C": 4}
    return encoding.get(sequence[kozak_start], 0)


def calculate_folding_energy_70(sequence):
    sequence_70 = sequence[:70]
    (ss, mfe) = RNA.fold(sequence_70)
    return "{:.2f}".format(mfe)


def calculate_folding_energy_80(sequence, aug_index):
    sequence_80 = sequence[aug_index - 40:aug_index + 43]
    (ss, mfe) = RNA.fold(sequence_80)
    return "{:.2f}".format(mfe)


def calculate_in_frame_aug(sequence, start_codon, stop_codon):
    sequence_cds = sequence[start_codon:stop_codon]
    num_in_frame_aug = sequence_cds.count("AUG")
    return num_in_frame_aug


def calculate_features(sequence, start_codon, stop_codon):
    if sequence:
        df = pd.DataFrame({'Sequence': [sequence]})
        df = df[df['Sequence'] != '']

        if not df.empty:
            df['CDS Length'] = stop_codon - start_codon
            df['Length of 5\' UTR'] = df['Sequence'].apply(lambda x: calculate_five_prime_utr(x, start_codon))
            df['Kozak pos. 1'] = df['Sequence'].apply(lambda x: calculate_kozak_pos_1(x, start_codon - 6))
            df['Kozak pos. 4'] = df['Sequence'].apply(lambda x: calculate_kozak_pos_4(x, start_codon - 3))
            df['Folding Energy 70'] = df['Sequence'].apply(calculate_folding_energy_70)
            df['Folding Energy 80'] = df['Sequence'].apply(lambda x: calculate_folding_energy_80(x, start_codon))
            df['in_frame AUG'] = df['Sequence'].apply(lambda x: calculate_in_frame_aug(x, start_codon, stop_codon))
            X = df[['CDS Length', 'Length of 5\' UTR', 'Kozak pos. 1', 'Kozak pos. 4', 'Folding Energy 70',
                    'Folding Energy 80', 'in_frame AUG']]
            return X
        else:
            return None


def evaluate_model(model, X_test):
    y_pred = model.predict(X_test)
    return y_pred


def main():
    st.markdown(
        """
        # Translation Initiation Rate Prediction App
        
        This app allows you to predict Translation Initiation Rate in Saccharomyces cerevisiae using Machine Learning methods.
        
        **Credits**
        - App built in `Python` + `Streamlit` by Sulagno Chakraborty, Inayat Ullah Irshad, Mahima, and Ajeet K. Sharma
        [[Read the Paper]]().
        ---
        """
    )

    start_codon = 'AUG'
    stop_codon = 'UAA'or 'UAG' or 'UGA'
    sequences = []
    start_codons = []
    stop_codons = []

    num_sequences = st.number_input("Number of Sequences", min_value=1, step=1, value=1)

    for i in range(num_sequences):
        col1, col2 = st.columns(2)
        sequence = col1.text_area(f"Sequence {i+1}")
        start_codon = col2.text_input(f"Start Codon Position {i+1}")
        stop_codon = col2.text_input(f"Stop Codon Position {i+1}")

        sequences.append(sequence)
        start_codons.append(start_codon)
        stop_codons.append(stop_codon)

    calculate_features_button = st.button("Calculate Features")
    predict_button = st.button("Predict")

    if calculate_features_button:
        df = pd.DataFrame()
        for sequence, start_codon, stop_codon in zip(sequences, start_codons, stop_codons):
            if sequence and start_codon and stop_codon:
                try:
                    start_codon = int(start_codon)
                    stop_codon = int(stop_codon)

                    X = calculate_features(sequence, start_codon, stop_codon)

                    if X is not None:
                        df = pd.concat([df, X])
                    else:
                        st.write(f"Invalid sequence or codon positions for Sequence {i+1}. "
                                 f"Please enter valid values.")

                except ValueError:
                    st.write(f"Invalid input for codon positions for Sequence {i+1}. "
                             f"Please enter numeric values.")
            else:
                st.write(f"Incomplete input for Sequence {i+1}. Please enter sequence and codon positions.")

        if not df.empty:
            st.subheader("Calculated Features:")
            st.write(df)
        else:
            st.write("No valid sequences found. Please enter valid values for all sequences.")

    if predict_button:
        df = pd.DataFrame()
        for sequence, start_codon, stop_codon in zip(sequences, start_codons, stop_codons):
            if sequence and start_codon and stop_codon:
                try:
                    start_codon = int(start_codon)
                    stop_codon = int(stop_codon)

                    X = calculate_features(sequence, start_codon, stop_codon)

                    if X is not None:
                        df = pd.concat([df, X])
                    else:
                        st.write(f"Invalid sequence or codon positions for Sequence {i+1}. "
                                 f"Please enter valid values.")

                except ValueError:
                    st.write(f"Invalid input for codon positions for Sequence {i+1}. "
                             f"Please enter numeric values.")
            else:
                st.write(f"Incomplete input for Sequence {i+1}. Please enter sequence and codon positions.")

        if not df.empty:
            model = pickle.load(open("tir_rf_model.pkl", "rb"))
            predictions = evaluate_model(model, df)

            st.subheader("Predictions:")
            st.write(predictions)
        else:
            st.write("No valid sequences found. Please enter valid values for all sequences.")


if __name__ == "__main__":
    main()
