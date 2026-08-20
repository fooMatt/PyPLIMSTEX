import argparse, tomllib
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor

from pyplimstex.functions import import_dynamx_csv, make_output_dir, remove_extra_charge_states, rename_exposure_entry
from pyplimstex.functions import setup_peptide_data, define_model, pseudo_bootstrap, monitor_progress

def plimstex(input, output, num_workers,
             renumber=0, protein_conc=1, 
             kd_init=0.1, d0_init=1, dd1_init=1, 
             outlier_threshold=0, bootstrap=100, 
             maxd_exists=False, normalise_maxd=False, do_ode=False, remove_badfits=False):
    # import data and make output folder
    df = import_dynamx_csv(csv_path=input)
    output_dir_path = make_output_dir(output_dir=output)

    # sanity check: cannot normalise if maxD doesn't exist
    if normalise_maxd and not maxd_exists:
        normalise_maxd = False
        print("[WARNING] Cannot normalise if maximum deuteration data is not provided!")
        print("[WARNING] Skipping normalisation...")

    # renumber residues
    if not isinstance(renumber, int):
        raise ValueError("[ERROR] Residue renumber offset must be an integer")
    
    df['Start'] = df['Start'] + renumber
    df['End'] = df['End'] + renumber 

    # clean data - remove extra charge states that user might have forgotten
    df = remove_extra_charge_states(data=df)

    # rename rows of t0, maxD (if it exists) and 0 ligand equivalents
    df = rename_exposure_entry(data=df, maxd_exists=maxd_exists)

    # set up dictionary to store peptide deuteration data from dataframe
    peptide_dict = setup_peptide_data(data=df, normalise=normalise_maxd)

    # define the deuteration model function (either ODE or analytical)
    deuteration_model = define_model(do_ODE=do_ode)

    # run fitting, parameter estimation, and plotting in parallel
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [
            executor.submit(pseudo_bootstrap,
                            deut_data=peptide_data,
                            peptide=peptide_name,
                            Ptot = protein_conc,
                            KD_init = kd_init, D0_init = d0_init, dD1_init =  dd1_init,
                            outlier_threshold = outlier_threshold,
                            deuteration_model = deuteration_model,
                            bootstrap = bootstrap,
                            output_dir = output_dir_path,
                            orig_df = df,
                            normalise = normalise_maxd,
                            remove_badfits = remove_badfits
                            ) 
                            for peptide_name, peptide_data in peptide_dict.items()
            ]

        # monitor progress of jobs
        loading_text = Path(__file__).parent / "loading.txt"
        monitor_progress(loading_text=loading_text, peptide_dict=peptide_dict, task_futures=futures)

        for f in futures:
            f.result()   

def launch_plimstex(config_path, num_workers):
    config_abs_path = Path(config_path).resolve()
    with open(config_abs_path, "rb") as file:
        config_file = tomllib.load(file)

    settings_dict = config_file["settings"]

    plimstex(**settings_dict, num_workers=num_workers)

def cli():
    p = argparse.ArgumentParser()
    p.add_argument('-c', '--config', type=str, required=True, help='Path to config.toml file')
    p.add_argument('-w', '--workers', type=int, required=False, default=4, help='Number of workers to spawn (default 4)')

    args = p.parse_args()
    launch_plimstex(config_path=args.config,
                    num_workers=args.workers
                    )

if __name__ == '__main__':
    cli()