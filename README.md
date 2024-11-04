# TERNIFY: Efficient sampling of PROTAC-induced ternary complexes

This repo accompanies the preprint ["Efficient sampling of PROTAC-induced ternary complexes"](https://www.biorxiv.org/content/10.1101/2024.10.30.619573v1).

### Setting up the environment 

Create a new conda environment and install rdkit:
`conda create -c conda-forge -n <your-env-name> rdkit`

Activate your environment and install the rest of the requirements:
`conda activate <your-env-name>`
`pip install -r requirements.txt`

### How to run

`chmod +x src/ternify.py`

`cd data/6hr2`

`../../src/ternify.py`

Analyze the results

`python ../../scripts/rmsd.py`

or 

`cd data/wdr5/7q2j`

`../../../src/ternify.py -p ../input/tcs_wdr5.inp`

Analyze the results

`python ../../../scripts/rmsd_m.py  poi_7q2j_aligned_on_5nvx_l_by_VHL.pdb`

### how to setup protein-protein interface
The box shall cover the entire PPIs interface of interest and extend at least 5 angstrom away from the surface of the anchroing protein. 

Check the parameter file (default: tcs.inp) in a data directory for preparation of input files and parameters

