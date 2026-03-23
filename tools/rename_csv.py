import os
import numpy as np
import pandas as pd

column_names = ['0', '11', '12', '15', '16', '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '32', '33',
                '34', '35', '37', '38', '39', '40', '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51',
                '52', '53', '54', '55', '56', '57', '59', '60', '61', '62', '63', '64', '65', '66',
                '67', '68', '69', '70', '71', '72', '73', '74', '75', '76', '77', '78', '79', '80', '81', '82', '83',
                '84', '85', '86', '87', '88', '89', '90', '91', '92', '93', '94']


def rename_csv_columns(input_file, output_file, column_names):

    # Load the CSV file
    df = pd.read_csv(input_file, sep=",", header=None)

    df.columns = column_names

    # Save the updated DataFrame to a new CSV file
    df.to_csv(output_file, sep=",", na_rep="NULL", index=False)


if __name__ == "__main__":

    dataset = "lha_data_2026-03-23"

    input_csv = f"./data/Germany/pydata/5315/{dataset}/lha_synthetic_data.csv"
    output_csv = f"./data/Germany/pydata/5315/{dataset}/lha_synthetic_data_renamed.csv"

    input_csv_altered_vacc = f"./data/Germany/pydata/5315/{dataset}/lha_synthetic_data_altered_vaccinations.csv"
    output_csv_altered_vacc = f"./data/Germany/pydata/5315/{dataset}/lha_synthetic_data_altered_vaccinations_renamed.csv"

    rename_csv_columns(input_csv, output_csv, column_names)

    rename_csv_columns(input_csv_altered_vacc,
                       output_csv_altered_vacc, column_names)

    print()
