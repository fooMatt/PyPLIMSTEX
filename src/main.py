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
p.add_argument('-i', '--input',    type=str,  required=True, help='Path to Cluster CSV from DynamX')
p.add_argument('-o', '--output',   type=str,  required=True, help='Output folder where results will be saved.')
p.add_argument('-r', '--renumber',   type=int,  required=False, default=0, help='Renumber amino acid residues with this offset')
p.add_argument('--maxd-exists',    action='store_true', help='Maximum deuteration data is included in DynamX analysis')
p.add_argument('--normalise-maxd', action='store_true', help='Normalise deuteration as a % of experimentally obtained maximum deuteration')

def main(input_data, output_dir, renumber, maxd_exists, normalise):
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
         maxd_exists=args.maxd_exists,
         normalise=args.normalise_maxd
        )