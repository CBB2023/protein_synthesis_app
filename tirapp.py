import streamlit as st
import pickle
import pandas as pd
import RNA


def calculate_five_prime_utr(sequence, start_codon):
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
            df['Kozak pos. 4'] = df['Sequence'].apply(lambda x: calculate_kozak_pos_4(x, start_codon + 3))
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
        [[Readthe Paper]]().
        ---
        """
    )

    sequence = st.text_area("Enter the sequence:")
    start_codon = st.text_input("Enter the position of the Start codon:")
    stop_codon = st.text_input("Enter the position of the Stop codon:")

    calculate_features_button = st.button("Calculate Features")
    predict_button = st.button("Predict")

    if calculate_features_button:
        if sequence and start_codon and stop_codon:
            try:
                start_codon = int(start_codon)
                stop_codon = int(stop_codon)

                X = calculate_features(sequence, start_codon, stop_codon)

                if X is not None:
                    st.subheader("Calculated Features:")
                    st.write(X)
                else:
                    st.write("Invalid sequence or codon positions. Please enter valid values.")
            except ValueError:
                st.write("Invalid input. Please enter numeric values for codon positions.")

    if predict_button:
        if sequence and start_codon and stop_codon:
            try:
                start_codon = int(start_codon)
                stop_codon = int(stop_codon)

                X = calculate_features(sequence, start_codon, stop_codon)

                if X is not None:
                    model = pickle.load(open("tir_rf_model.pkl", "rb"))
                    predictions = evaluate_model(model, X)

                    st.subheader("Predictions:")
                    st.write(predictions)
                else:
                    st.write("Invalid sequence or codon positions. Please enter valid values.")
            except ValueError:
                st.write("Invalid input. Please enter numeric values for codon positions.")


if __name__ == "__main__":
    main()
