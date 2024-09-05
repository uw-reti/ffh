import numpy as np
import matplotlib.pyplot as plt

import openmc
import sys
import os
import openmc.mgxs as mgxs
from datetime import datetime
import pandas as pd
import argparse
import shutil
import csv

# Script to convert csv created from this script to a GeN-Foam usable external source: https://github.com/pshriwise/openmc-src-mesh-mapping

# Create argparser for some non-required variables such as filepath
parser = argparse.ArgumentParser(prog="srcMap_CSV2GF", 
                                 description="Script to convert csv created from this script to a GeN-Foam usable external source: https://github.com/pshriwise/openmc-src-mesh-mapping", 
                                 epilog="V0.1")
parser.add_argument('-f', '--filepath', help='filepath where csv file lives and output should be written', default="./")
# parser.add_argument('-n', '--filename', help='filename where the OpenMC produced source csv exists', default="mesh_src_strengths.csv")
parser.add_argument('-df', '--defaultFluxFile', help='Filepath of default flux file to use for last section')
args=parser.parse_args()

if args.filepath[-1] != '/':
    filePath=args.filepath+'/'
else:
    filePath=args.filepath

csv_data = np.genfromtxt(f"{filePath}mesh_src_strengths.csv", delimiter=',',dtype=None)
numElements=np.size(csv_data,0)-1


f = open(f"{filePath}defaultExternalSourceFlux","w")

# Use existing defaultFlux file as a template if argmuent exists
if args.defaultFluxFile is not None:
    defaultFile= open(f"{args.defaultFluxFile}", "r")
    endSection=False
    for lines in defaultFile.readlines():
        if lines[0:13] == "internalField":
            f.write("internalField   nonuniform List<scalar>\n")
            f.write(f"{numElements}\n")
            f.write("(\n")
            for i in range(0,numElements):
                # get volumetric strength in 1/m^3-s
                volStrength = (10**6)*(float(csv_data[i+1,6])/float(csv_data[i+1,4]))
                f.write(f"{volStrength}\n")

            f.write(")\n")
            f.write(";\n\n")
        else:
            f.write(lines)

else:
    f.write("boundaryField \n\{\nTODO\n\}")

    f.write("/*--------------------------------*- C++ -*----------------------------------*\ \n")
    f.write("|       ______          _   __           ______                               |\n")
    f.write("|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     |\n")
    f.write("|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    |\n")
    f.write("|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    |\n")
    f.write("|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     |\n")
    f.write("|    Copyright (C) 2015 - 2022 EPFL                                           |\n")
    f.write("\*---------------------------------------------------------------------------*/\n")
    f.write("FoamFile\n")
    f.write("{\n")
    f.write("    version     2.0;\n")
    f.write("    format      ascii;\n")
    f.write("    class       volScalarField;\n")
    f.write("    object      defaultFlux;\n")
    f.write("}\n")
    f.write("// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n")
    f.write("dimensions      [ 0 -3 -1 0 0 0 0 ];\n\n")

    f.write("internalField   nonuniform List<scalar>\n")
    f.write(f"{numElements}\n")
    f.write("(\n")
    for i in range(0,numElements):
        # get volumetric strength in 1/m^3-s
        volStrength = (10**6)*(float(csv_data[i+1,6])/float(csv_data[i+1,4]))
        f.write(f"{volStrength}\n")

    f.write(")\n")
    f.write(";\n\n")


f.close()