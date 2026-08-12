import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.integrate import solve_ivp
from sklearn.metrics import r2_score
import argparse
from concurrent.futures import ProcessPoolExecutor

from functions import import_dynamx_csv, make_output_dir, remove_extra_charge_states, rename_exposure_entry

# === parsing arguments ===
p = argparse.ArgumentParser()
p.add_argument('-i', '--input',         type=str,  required=True, help='Path to Cluster CSV from DynamX')
p.add_argument('-o', '--output',        type=str,  required=True, help='Output folder where results will be saved.')
p.add_argument('-r', '--renumber',      type=int,  required=False, default=0, help='Renumber amino acid residues with this offset')
p.add_argument('-pc', '--protein-conc', type=float,  required=True, help='Total protein concentration in deuterated sample (in μM)')
p.add_argument('-kd', '--kd-init',      type=float,  required=True, help='Best initial estimate of the dissociation constant (in μM)')
p.add_argument('-d0', '--d0-init',      type=float,  required=True, help='Best initial estimate of the deuteration amount when no ligand is present (in Da if not normalising, else in % uptake)')
p.add_argument('-dd1', '--dd1-init',    type=float,  required=True, help='Best initial estimate of the difference in deuteration between D0 and the fully saturated system (in Da if not normalising, else in % uptake)')
p.add_argument('--outlier-threshold',   type=float,  required=False, default=None, help='How tolerant you are of outliers in the data, the higher the number, the more tolerant you are. (default None i.e. will not remove outliers)')
p.add_argument('--bootstrap',           type=int,    required=False, default=100, help='Number of iterations for the pseudo-bootstrapping of each peptide (default 100)')
p.add_argument('--maxd-exists',      action='store_true', help='Maximum deuteration data is included in DynamX analysis')
p.add_argument('--normalise-maxd',   action='store_true', help='Normalise deuteration as a % of experimentally obtained maximum deuteration')
p.add_argument('--do-ode',           action='store_true', help="Solve for free ligand concentration using scipy's solve_ivp Radau method instead of an analytical solution")

def main(input_data, output_dir, renumber, prot_conc, kd_init, d0_init, dd1_init, threshold, bootstrap, maxd_exists, normalise, do_ode):
    # import data and make output folder
    df = import_dynamx_csv(input_data)
    output_dir_path = make_output_dir(output_dir)

    # sanity check: cannot normalise if maxD doesn't exist
    if normalise and not maxd_exists:
        normalise = False
        print("[WARNING] Cannot normalise if maximum deuteration data is not provided!")
        print("[WARNING] Skipping normalisation...")

    # renumber residues
    if not isinstance(renumber, int):
        raise ValueError("[ERROR] Residue renumber offset must be an integer")
    
    df['Start'] = df['Start'] + renumber
    df['End'] = df['End'] + renumber 

    # clean data - remove extra charge states that user might have forgotten
    df = remove_extra_charge_states(df)

    # rename rows of t0, maxD (if it exists) and 0 ligand equivalents
    df = rename_exposure_entry(df, maxd_exists)



if __name__ == '__main__':
    args = p.parse_args()
    main(input_data=args.input,
         output_dir=args.output,
         renumber=args.renumber,
         prot_conc=args.protein_conc,
         kd_init=args.kd_init,
         d0_init=args.d0_init,
         dd1_init=args.dd1_init,
         threshold=args.outlier_threshold,
         bootstrap=args.bootstrap,
         maxd_exists=args.maxd_exists,
         normalise=args.normalise_maxd,
         do_ode=args.do_ode
        )