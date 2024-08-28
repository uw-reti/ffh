#!/bin/bash

# Extract the files from the zipped file (maybe change to endf8 at some point)
tar -xf endfb71_hdf5.xz

# If necessary, configure environment variables
export OPENMC_CROSS_SECTIONS=endfb71_hdf5/cross_sections.xml

ls

# Run OpenMC (via python) with 16 cores
python3 mgxs_MOSART_chtc.py

ls

# zip the mgxs directory
tar -czf mgxs_out_MOSART.tar.gz mgxs

# move files to staging directory (they aren't automatically transferred back to home)
cp mgxs_out_MOSART.tar.gz /staging/jeickman

# Print a list of environment variables
env

ls

# Sleep
sleep 3m
