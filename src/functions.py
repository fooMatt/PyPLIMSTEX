import pandas as pd
import numpy as np
from pathlib import Path

def import_dynamx_csv(csv_path):
    """
    Import Cluster CSV exported from DynamX data analysis
    Returns pandas dataframe
    """
    csv_path_abs = Path(csv_path).resolve()
    data = pd.read_csv(csv_path_abs)
    return data

def make_output_dir(output_dir):
    """
    Makes output directory and returns its absolute path
    """
    output_dir_abs = Path(output_dir).resolve()
    output_dir_abs.mkdir(parents=True, exist_ok=True)
    return output_dir_abs

def remove_extra_charge_states(data):
    """
    Create a new series counting the number of occurrences of each charge state for each peptide
    in case the user forgets to deselect the other charge states in the DynamX
    """
    most_frequent_z = data.groupby('Sequence')['z'].apply(lambda x: x.value_counts().idxmax())

    # create a boolean mask: for each row, check if the 'z' value in that row is equal to the most frequent 'z' for that peptide
    clean_data = data[data.apply(lambda row: row['z'] == most_frequent_z[row['Sequence']], axis=1)]

    return clean_data

def rename_exposure_entry(data, maxd_exists=False):
    """
    Rename rows of t0, maxD and 0 ligand equivalents in the dataframe
    """
    data['Exposure'] = data['Exposure'].astype(str)

    if maxd_exists:
        data.loc[data['Exposure'] == list(data['Exposure'].unique())[0], 'Exposure'] = 't0'
        data.loc[data['Exposure'] == list(data['Exposure'].unique())[1], 'Exposure'] = 'maxD'
        data.loc[data['Exposure'] == list(data['Exposure'].unique())[2], 'Exposure'] = '0'

    else:
        data.loc[data['Exposure'] == list(data['Exposure'].unique())[0], 'Exposure'] = 't0'
        data.loc[data['Exposure'] == list(data['Exposure'].unique())[1], 'Exposure'] = '0'

    return data

def setup_peptide_data(data, normalise=False):
    """
    Reads pandas dataframe and calculates average t0 deuteration
    For each 'equivalent' of ligand, calculates the delta deuteration (i.e. shift in deuteration compared to t0)
    Returns nested dictionary: for each peptide, lists delta deuteration for every ligand equivalent tested (or at least present in the data)
    Each peptide will eventually produce 1 titration curve plotting delta deuteration against ligand equivalent  
    """
    # calculate average t0 for each sequence
    t0_all = []
    for i in data['Sequence'].unique():
        t0_cum = []
        for j in data[(data['Sequence'] == i) & (data['Exposure'] == 't0')]['Center']:
            t0_cum.append(j)
            t0_avg = np.mean(t0_cum)
        t0_all.append(t0_avg)

    # create a list of unique peptides and ligands
    peptide_column = list(data['Sequence'].unique())
    ligand_eq = list(data['Exposure'].unique())[2:]

    peptide_data = {}

    # create framework for peptide data
    for i in peptide_column:
        ligand_dict = {j:[] for j in ligand_eq}
        peptide_data[i] = ligand_dict

    # fill in the peptide data
    for i in peptide_column:
        current_seq = data[data['Sequence'] == i]
        for j in ligand_eq:
            current_eq = current_seq[current_seq['Exposure'] == j]
            for k in current_eq['Center']:
                peptide_data[i][j].append(k - t0_all[peptide_column.index(i)])

    # if normalisation is true...
    if normalise:
        # loop to calculate average maxD for each sequence
        maxD_all = []
        for i in data['Sequence'].unique():
            maxD_cum = []
            for j in data[(data['Sequence'] == i) & (data['Exposure'] == 'maxD')]['Center']:
                maxD_cum.append(j)
                maxD_avg = np.mean(maxD_cum)
            maxD_all.append(maxD_avg)

        # normalise the data by maxD
        peptide_data_maxD = {}

        # create framework for normalised peptide data
        for i in peptide_column:
            ligand_dict = {j:[] for j in ligand_eq}
            peptide_data_maxD[i] = ligand_dict

        # fill in the normalised peptide data
        for i in peptide_column:
            current_seq = data[data['Sequence'] == i]
            for j in ligand_eq:
                current_eq = current_seq[current_seq['Exposure'] == j]
                for k in current_eq['Center']:
                    norm_value = 100 * (k - t0_all[peptide_column.index(i)]) / (maxD_all[peptide_column.index(i)] - t0_all[peptide_column.index(i)])
                    peptide_data_maxD[i][j].append(norm_value)
        
        peptide_data = peptide_data_maxD

    return peptide_data

def ode_model():
    pass

def analytic_model():
    pass

def fit_model():
    pass

def remove_outliers():
    pass

def plot_plimstex_curve():
    pass

def pseudo_bootstrap():
    pass

