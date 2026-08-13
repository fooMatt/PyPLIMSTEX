import pandas as pd
import numpy as np
from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.optimize import minimize

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

def define_model(do_ODE):
    """
    Set up the model to calculate total deuteration as a function of bound ligand concentration (i.e. conc of protein-ligand complex) 
    Parameters for the model are KD (dissociation constant), D0 (deuteration with no ligand), deltaD1 (), Ltot (total ligand concentration), and Ptot (total protein concentration)
    Following Zhu et al. (2004), bound ligand concentration can be calculated by solving an ODE
    In our 1:1 stoichiometry, we can also solve it analytically... this is computationally faster
    ODE version is included for compatibility with the original PLIMSTEX method and also in case we expand it to 1:N stoichiometries in the future
    """
    if do_ODE:
        def freeligand(Ltot, Lfree, Ptot, KD):
            dLfreedLtot = (KD + Lfree) / (KD + 2*Lfree + Ptot - Ltot)
            return dLfreedLtot

        def totaldeut(Ltot, Ptot, D0, deltaD1, KD):
            # solve ODE to determine Lfree for each value of Ltot
            Lfree_init = 0 # if Ltot is 0, then obviously Lfree is 0
            Ltot_span = [Ltot[0], Ltot[-1]]
            Ltot_eval = np.linspace(Ltot_span[0], Ltot_span[1], 1000)

            sol = solve_ivp(freeligand, t_span=Ltot_span, y0=[Lfree_init],
                            t_eval=Ltot_eval, args=(Ptot, KD), method='Radau', dense_output=True)

            # create lookup lists with Lfree values for each Ltot value
            Lfree_lookup = sol.y[0]
            Ltot_lookup = sol.t
            Deut = []
            # check which Ltot value in lookup list is closest to actual Ltot and return corresponding Lfree
            Lfree_est = []
            for L in Ltot:
                idx = (np.abs(Ltot_lookup - L)).argmin()
                Lfree_est.append(Lfree_lookup[idx])
            Lfree_est = np.array(Lfree_est)
            Deut = D0 - deltaD1 * ((Ltot - Lfree_est) / Ptot)
            return Deut

    else:
        def totaldeut(Ltot, Ptot, D0, deltaD1, KD):
            Ptot = float(Ptot)
            D0 = float(D0)
            deltaD1 = float(deltaD1)
            KD = float(KD)
            Lfree = []
            for L in Ltot:
                L = float(L)
                Lfree_calc = ((L-KD-Ptot)+np.sqrt((Ptot-L+KD)**2 + 4*KD*L))/2
                Lfree.append(float(Lfree_calc))
            Lfree = np.array(Lfree)
            Ltot = np.array(Ltot)
            Deut = D0 - deltaD1 * ((Ltot - Lfree) / Ptot)
            return Deut
        
    return totaldeut

def fit_model(params, myargs):
    """
    Fits the model to the experimental data by varying the parameters D0, deltaD1 and KD in order to minimise the mean squared error from the observed deuteration
    Requires initial estimates of these parameters as a starting point
    Uses the L-BFGS-B algorithm for optimisation from scipy.optimize
    No bounds on the parameters except for KD which must be positive by definition (since it is a ratio of concentrations)
    """
    # params should be (D0, deltaD1, KD)
    # model is defined in totaldeut()
    def mse(params, model, data, Ltot_exp, Ptot):
        predict = model(Ltot_exp, Ptot, params[0], params[1], params[2])
        sq_err = [(predict[i] - data[i])**2 for i in range(len(predict))]
        return np.mean(sq_err)

    # params is a list of the initial estimates
    # myargs is a tuple of the model, data, Ltot_exp and Ptot
    init_est = np.array(params)
    res = minimize(mse, init_est, args=myargs, method='L-BFGS-B', bounds=[(None,None),(None,None),(0,None)])

    # res.x should contain 3 values - D0, deltaD1 and KD
    return res.x, res.fun
    
def remove_outliers(data, threshold):
  """
  Removes outliers for each ligand equivalent and for each peptide
  Roughly, outliers are data points (an observed deuterium uptake) that lie far from the rest
  This is done by comparing the ratios of distances between the two extreme points and their neighbours;
  if the difference between the ratios is larger in magnitude than the threshold, the furthest point is removed.
  Repeats until stable (no more points removed) or fewer than 3 points remain.
  """
  threshold = float(threshold)
  clean_data = []

  for i in data.values():
    i = tuple(sorted(i))  # sorts in ascending order

    while len(i) >= 3:
      gap_low = i[1] - i[0]      # gap between smallest two points
      gap_high = i[-1] - i[-2]   # gap between largest two points
      span = i[-1] - i[0]        # total range

      gl, gh = (g or 1e-8 for g in (gap_low, gap_high))  # avoid division by 0

      if (span/gl) - (span/gh) > threshold:
        i = i[:-1]   # drop the max, re-check
      elif (span/gl) - (span/gh) < -threshold:
        i = i[1:]    # drop the min, re-check
      else:
        break        # stable, stop removing

    clean_data.append(i)

  return clean_data

def plot_plimstex_curve():
    pass

def pseudo_bootstrap():
    pass

def loading_message():
    """
    Prints regular updates while PyPLIMSTEX analysis and plotting is taking place
    """
    pass