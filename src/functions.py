import pandas as pd
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

def setup_peptide_data():
    pass

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

