## PLIMSTEX
### Python implementation of PLIMSTEX

This is a Python-based implementation of the approach to quantifying protein-ligand interaction by mass spectrometry, titration and H/D exchange (PLIMSTEX) proposed by [Zhu *et al.*, 2004](https://doi.org/10.1016/j.jasms.2003.11.007).

For now, only 1:1 protein-ligand stoichiometry can be analysed with this code, which was developed to analyse interactions of GPCRs with G proteins, a known 1:1 stoichiometry. The theoretical basis of 1:N stoichiometry is outlined in the 2004 paper by Zhu *et al.*. 

---
### Installation
#### Pre-requisites
**Python 3.11** or higher

---
#### Installing in a clean Conda environment 
   ```bash
   conda create -y -n pyplimstex python=3.11
   conda activate pyplimstex
   pip install git+https://github.com/fooMatt/PyPLIMSTEX.git
   ```

### Quick start
   You can modify to the template `config.toml` file in `assets/`.
   
   Note that the user still needs to 'pre-process' the HDX-MS data in DynamX and export this as a cluster CSV file as this code is unable to read the raw HDX-MS data.

   Then run:
   ```bash
   pyplimstex --config path/to/config.toml (optionally: --workers NUM_WORKERS)
   ```

---
Made at the *Institut de Génomique Fonctionnelle*, Montpellier

Granier-Mouillac Team

Matthew Chee, 2026