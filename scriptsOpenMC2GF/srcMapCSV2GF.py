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
parser.add_argument('-df', '--defaultFluxFile', help='Filepath of default flux file to use for formatting and boundary field section')
args=parser.parse_args()

if args.filepath[-1] != '/':
    filePath=args.filepath+'/'
else:
    filePath=args.filepath

csv_data = np.genfromtxt(f"{filePath}mesh_src_strengths.csv", delimiter=',',dtype=None, encoding=None)
numElements=np.size(csv_data,0)-1


fileString=""

# Use existing defaultFlux file as a template if argmuent exists
if args.defaultFluxFile is not None:
    defaultFile= open(f"{args.defaultFluxFile}", "r")
    endSection=False
    for lines in defaultFile.readlines():
        if lines[0:13] == "internalField":
            fileString+="internalField   nonuniform List<scalar>\n"
            fileString+=f"{numElements}\n"
            fileString+="(\n"
            for i in range(0,numElements):
                # get volumetric strength in 1/m^3-s
                volStrength = (10**6)*(float(csv_data[i+1,6])/float(csv_data[i+1,4]))
                fileString+=f"{volStrength}\n"

            fileString+=")\n"
            fileString+=";\n\n"
        else:
            fileString+=lines
# If default flux template is not used, create from scratch
# There are "TODO"s in places where more information needs to be added
else:
    fileString+="/*--------------------------------*- C++ -*----------------------------------*\ \n"
    fileString+="|       ______          _   __           ______                               | \n"
    fileString+="|      / ____/  ___    / | / /          / ____/  ____   ____ _   ____ ___     | \n"
    fileString+="|     / / __   / _ \  /  |/ /  ______  / /_     / __ \ / __ `/  / __ `__ \    | \n"
    fileString+="|    / /_/ /  /  __/ / /|  /  /_____/ / __/    / /_/ // /_/ /  / / / / / /    | \n"
    fileString+="|    \____/   \___/ /_/ |_/          /_/       \____/ \__,_/  /_/ /_/ /_/     | \n"
    fileString+="|    Copyright (C) 2015 - 2022 EPFL                                           | \n"
    fileString+="\*---------------------------------------------------------------------------*/ \n"
    fileString+="FoamFile\n"
    fileString+="{\n"
    fileString+="    version     2.0;\n"
    fileString+="    format      ascii;\n"
    fileString+="    class       volScalarField;\n"
    fileString+="    object      defaultFlux;\n"
    fileString+="}\n"
    fileString+="// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n"
    fileString+="dimensions      [ 0 -3 -1 0 0 0 0 ];\n\n"

    fileString+="internalField   nonuniform List<scalar>\n"
    fileString+=f"{numElements}\n"
    fileString+="(\n"
    for i in range(0,numElements):
        # get volumetric strength in 1/m^3-s
        volStrength = (10**6)*(float(csv_data[i+1,6])/float(csv_data[i+1,4]))
        fileString+=f"{volStrength}\n"

    fileString+=")\n"
    fileString+=";\n\n"
    fileString+="boundaryField\n{\nTODO\n}\n\n"
    fileString+="// ************************************************************************* //"

f = open(f"{filePath}defaultExternalSourceFlux","w")
f.write(fileString)
f.close()
