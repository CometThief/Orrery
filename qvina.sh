#!/bin/bash
set -e

# Filenames from Python
ligand_mol2=$1
ligand_pdbqt=$2
ligand_out_pdbqt=$3

# Prepare both files
source activate mgltools
prepare_ligand4.py -l $ligand_mol2
prepare_receptor4.py -r receptor_65.pdb

# Run QVina2
source activate Orrery
qvina2 --config config.vd --ligand $ligand_pdbqt
obabel $ligand_out_pdbqt -O output.mol2 -m

conda deactivate