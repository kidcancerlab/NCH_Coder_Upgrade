#!/usr/bin/bash

set -e

if [ ! -d input/nichenetr ]
then
    mkdir -p input/nichenetr
fi

cd input/nichenetr/

# Download data from https://zenodo.org/record/3260758

wget \
    -O ligand_target_matrix_mouse.rds \
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds"
wget \
    -O lr_network_mouse.rds \
    "https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds"
wget \
    -O weighted_networks_mouse.rds \
    "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds"
