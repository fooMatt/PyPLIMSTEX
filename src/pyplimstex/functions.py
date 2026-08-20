import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random, time, csv
from pathlib import Path
from scipy.integrate import solve_ivp
from scipy.optimize import minimize
from sklearn.metrics import r2_score

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
    unique_exposures = sorted(data['Exposure'].unique(), key=float)

    data['Exposure'] = data['Exposure'].astype(str)

    if maxd_exists:
        data.loc[data['Exposure'] == str(unique_exposures[0]), 'Exposure'] = 't0'
        data.loc[data['Exposure'] == str(unique_exposures[1]), 'Exposure'] = 'maxD'
        data.loc[data['Exposure'] == str(unique_exposures[2]), 'Exposure'] = '0'

    else:
        data.loc[data['Exposure'] == str(unique_exposures[0]), 'Exposure'] = 't0'
        data.loc[data['Exposure'] == str(unique_exposures[1]), 'Exposure'] = '0'

    return data

def setup_peptide_data(data, normalise=False):
    """
    Reads pandas dataframe and calculates average t0 deuteration
    For each 'equivalent' of ligand, calculates the delta deuteration (i.e. shift in deuteration compared to t0)
    Returns nested dictionary: for each peptide, lists delta deuteration for every ligand equivalent tested (or at least present in the data)
    Each peptide will eventually produce 1 titration curve plotting delta deuteration against ligand equivalent
    Peptides missing t0 (or maxD, if normalising) data are skipped with a warning rather than aborting the run
    """
    t0_all = {}
    skipped_peptides = set()
    for i in data['Sequence'].unique():
        t0_cum = list(data[(data['Sequence'] == i) & (data['Exposure'] == 't0')]['Center'])
        if len(t0_cum) == 0:
            print(f"[WARNING] No t0 data found for peptide {i} - skipping this peptide")
            skipped_peptides.add(i)
            continue
        t0_all[i] = np.mean(t0_cum)

    if normalise:
        maxD_all = {}
        for i in data['Sequence'].unique():
            if i in skipped_peptides:
                continue
            maxD_cum = list(data[(data['Sequence'] == i) & (data['Exposure'] == 'maxD')]['Center'])
            if len(maxD_cum) == 0:
                print(f"[WARNING] No maxD data found for peptide {i} - skipping this peptide")
                skipped_peptides.add(i)
                continue
            maxD_all[i] = np.mean(maxD_cum)

    # exclude skipped peptides from further processing
    peptide_column = [i for i in data['Sequence'].unique() if i not in skipped_peptides]
    ligand_eq = list(data['Exposure'].unique())[2:]

    peptide_data = {}
    for i in peptide_column:
        peptide_data[i] = {j: [] for j in ligand_eq}

    for i in peptide_column:
        current_seq = data[data['Sequence'] == i]
        for j in ligand_eq:
            current_eq = current_seq[current_seq['Exposure'] == j]
            for k in current_eq['Center']:
                peptide_data[i][j].append(k - t0_all[i])

    if normalise:
        peptide_data_maxD = {}
        for i in peptide_column:
            peptide_data_maxD[i] = {j: [] for j in ligand_eq}

        for i in peptide_column:
            current_seq = data[data['Sequence'] == i]
            for j in ligand_eq:
                current_eq = current_seq[current_seq['Exposure'] == j]
                for k in current_eq['Center']:
                    norm_value = 100 * (k - t0_all[i]) / (maxD_all[i] - t0_all[i])
                    peptide_data_maxD[i][j].append(norm_value)

        peptide_data = peptide_data_maxD

    return peptide_data

def _totaldeut_ode(Ltot, Ptot, D0, deltaD1, KD):
    Lfree_init = 0
    Ltot_span = [Ltot[0], Ltot[-1]]
    Ltot_eval = np.linspace(Ltot_span[0], Ltot_span[1], 1000)

    sol = solve_ivp(_freeligand, t_span=Ltot_span, y0=[Lfree_init],
                    t_eval=Ltot_eval, args=(Ptot, KD), method='Radau', dense_output=True)

    Lfree_lookup = sol.y[0]
    Ltot_lookup = sol.t
    Lfree_est = []
    for L in Ltot:
        idx = (np.abs(Ltot_lookup - L)).argmin()
        Lfree_est.append(Lfree_lookup[idx])
    Lfree_est = np.array(Lfree_est)
    Deut = D0 - deltaD1 * ((Ltot - Lfree_est) / Ptot)
    return Deut

def _freeligand(Ltot, Lfree, Ptot, KD):
    dLfreedLtot = (KD + Lfree) / (KD + 2*Lfree + Ptot - Ltot)
    return dLfreedLtot

def _totaldeut_analytical(Ltot, Ptot, D0, deltaD1, KD):
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

def define_model(do_ODE):
    """
    Returns the appropriate module-level deuteration model function (ODE-based or analytical),
    both defined outside this function so they remain picklable for ProcessPoolExecutor workers
    """
    return _totaldeut_ode if do_ODE else _totaldeut_analytical

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

def pseudo_bootstrap(deut_data, peptide, Ptot, KD_init, D0_init, dD1_init, 
                     outlier_threshold, deuteration_model, bootstrap,
                     output_dir, orig_df, normalise=False, remove_badfits=False):
    """
    Runs pseudo-bootstrap for one peptide to estimate KD 
    by fitting model deuteration to experimental deuteration
    Calculates average KD and standard deviation over all runs
    Plots all runs on one plot and saves to output directory 
    """ 
    # data here refers to deuteration data for a given peptide
    ligand_conc_eqs = [float(eq) for eq in list(deut_data.keys())]
    ligand_conc_eqs = np.array(ligand_conc_eqs)

    # convert to absolute concentrations in µM
    ligand_conc_abs = [float(eq) * Ptot for eq in list(deut_data.keys())] 
    ligand_conc_abs = np.array(ligand_conc_abs)

    init_estimates = [D0_init, dD1_init, KD_init] # initial estimates for D0, deltaD1 and KD

    cleaned_data = remove_outliers(deut_data, outlier_threshold)

    # define function for modelling deuteration
    totaldeut = deuteration_model

    fig, ax = plt.subplots()

    KD_list = []
    r2_list = []
    completed_runs = 0
    attempts = 0
    max_attempts = bootstrap * 10  # tune as needed

    while completed_runs < bootstrap:
        attempts += 1
        if attempts > max_attempts:
            raise RuntimeError(
                f"[ERROR] {peptide}: only {completed_runs}/{bootstrap} fits passed "
                f"R² >= 0.25 after {max_attempts} attempts. Check initial estimates or data quality."
            )

        deuteration_exp = np.array([np.random.choice(t) for t in cleaned_data])
        estimate = fit_model(init_estimates, (totaldeut, deuteration_exp, ligand_conc_abs, Ptot))
        
        # extract parameter estimates
        D0_est = estimate[0][0]
        deltaD1_est = estimate[0][1]
        KD_est = estimate[0][2]

        # calculate model-estimated deuteration for given ligand concentrations 
        deut_est = totaldeut(ligand_conc_abs, Ptot, D0_est, deltaD1_est, KD_est)

        # calculate r-squared (goodness of fit of model)
        r_squared = r2_score(deuteration_exp, deut_est)

        if remove_badfits and r_squared < 0.25:
            continue

        r2_list.append(r_squared)
        KD_list.append(KD_est)
        completed_runs += 1

        # plotting fitted curve
        deut_est_plotted = totaldeut(np.linspace(ligand_conc_abs[0],ligand_conc_abs[-1]), Ptot, D0_est, deltaD1_est, KD_est)

        ax.plot(ligand_conc_eqs, deuteration_exp, 'o', label='Data', color='black', alpha=0.5)
        ax.plot(np.linspace(ligand_conc_eqs[0],ligand_conc_eqs[-1]), deut_est_plotted, '-', label='Fit', alpha=0.75)
        ax.set_xlabel('[Ligand total]/[Protein total]')
        if normalise:
            ax.set_ylabel('Relative deuteration (% of max D)')
        else:
            ax.set_ylabel('Absolute deuteration (Da)')
        ax.set_title(f'Fitted PLIMSTEX curve for {peptide}')

    KD_avg = np.average(KD_list)
    KD_sd = np.std(KD_list)
    r2_avg = np.average(r2_list)

    ax.text(0.5, -0.15, f"Average R-squared: {r2_avg:.3f}\nAverage KD: {KD_avg:.3f} ± {KD_sd:.3f} µM",
            horizontalalignment='center',
            verticalalignment='top',
            transform=ax.transAxes) # Use transAxes to position text relative to the axes

    # Save the plot to a PNG file
    pep_start = str(list(set(orig_df[orig_df['Sequence']==peptide]['Start'].unique()))[0])
    pep_end = str(list(set(orig_df[orig_df['Sequence']==peptide]['End'].unique()))[0])
    peptide_name = pep_start+"-"+pep_end+"_"+peptide
    figure_path = Path(output_dir).resolve() / f"PLIMSTEX fit for {peptide_name}.png"
    fig.savefig(figure_path, bbox_inches='tight') # Use bbox_inches='tight' to include the text
    plt.close(fig)

    # also save KD_list and r2_list for further analysis
    rows = zip(KD_list, r2_list)
    save_csv_path = Path(output_dir).resolve() / f"Estimated KD and R-squared of fits for {peptide_name}.csv"
    with open(save_csv_path, "w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(["Estimated KD", "R-squared of fit"])
        writer.writerows(rows)
    
def monitor_progress(loading_text, peptide_dict, task_futures):
    """
    Prints regular updates (every 15 seconds) while PyPLIMSTEX analysis and plotting is taking place
    """
    loading_text = Path(loading_text).resolve()
    lines = Path(loading_text).read_text(encoding="utf-8").splitlines()

    total_peptides = len(peptide_dict.keys())

    first_iteration = True

    while True:
        succeeded = 0
        failed = 0
        cancelled = 0

        # Inspect state of each future
        for f in task_futures:
            if f.done():
                if f.cancelled():
                    cancelled += 1
                elif f.exception() is not None:
                    failed += 1
                else:
                    succeeded += 1

        finished_peptides = succeeded + failed + cancelled
        selected = random.choice(lines)   

        if not first_iteration:
            # \033[F moves cursor UP 1 line. \033[4F moves cursor UP 4 lines.
            print("\033[6F", end="")
        else:
            first_iteration = False

        # Print your 4 lines as normal (use \033[K to clear any leftover characters per line)
        print(f"\033[KFinished processing {finished_peptides} / {total_peptides} peptides...")
        print(f"\033[K  • Succeeded: {succeeded}")
        print(f"\033[K  • Failed:    {failed}")
        print(f"\033[K  • Cancelled: {cancelled}")
        print("\033[K=====")
        print(f"\033[KDid you know? {selected}")
        
        if finished_peptides == total_peptides:
            break

        time.sleep(15)

    print("Finished processing all peptides!")