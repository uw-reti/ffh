# 0-D Transport Solver for OpenMC Validation
import scipy
import numpy as np
import pandas as pd
from scipy.optimize import fsolve

# Define Known Variables from OpenMC Output

# Unit Conversion
cmm = 100

# Phi_1
p1 = 1

# Absorption
absp = [0, 0, 0, 0, 0, 0]
# df = pd.read_csv(r"\\wsl$\Ubuntu\home\jeickman\testbin\OpenMC\for_joey_8_22\cross_section_gen\mgxs\absorption-xs.csv")
df = pd.read_csv("~/Desktop/ffh/scriptsOpenMC2GF/testing/MOSART_fuel/mgxs/absorption-xs.csv")
i=0
for row in df['mean']:
    absp[i] = row*cmm
    i += 1 
print('absorption: ', absp)
# Scattering
sct = [[0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0]]
# df = pd.read_csv(r"\\wsl$\Ubuntu\home\jeickman\testbin\OpenMC\for_joey_8_22\cross_section_gen\mgxs\scatteringmatrix-xs.csv")
df = pd.read_csv("~/Desktop/ffh/scriptsOpenMC2GF/testing/MOSART_fuel/mgxs/scatteringmatrix-xs.csv")
i=0
j=0
for row in df['mean']:
    sct[i][j] = row*cmm
    j += 1
    if j == 6:
        i += 1 
        j = 0
print('scattering: ',sct)
# NuFission
nf = [0, 0, 0, 0, 0, 0]
# df = pd.read_csv(r"\\wsl$\Ubuntu\home\jeickman\testbin\OpenMC\for_joey_8_22\cross_section_gen\mgxs\nu-SigmaEff.csv")
df = pd.read_csv("~/Desktop/ffh/scriptsOpenMC2GF/testing/MOSART_fuel/mgxs/nu-SigmaEff.csv")
i=0
for row in df['mean']:
    nf[i] = row*cmm
    i += 1 
print('nufission: ', nf)
# Chi
chi = [0, 0, 0, 0, 0, 0]
# df = pd.read_csv(r"\\wsl$\Ubuntu\home\jeickman\testbin\OpenMC\for_joey_8_22\cross_section_gen\mgxs\chi.csv")
df = pd.read_csv("~/Desktop/ffh/scriptsOpenMC2GF/testing/MOSART_fuel/mgxs/chi.csv")
i=0
for row in df['mean']:
    chi[i] = row
    i += 1 
print('chi: ', chi)
# Define Function and Equations for 6 groups
def zeroD(unks):
    p2, p3, p4, p5, p6, k = unks
    return [(absp[0]+sct[0][1]+sct[0][2]+sct[0][3]+sct[0][4]+sct[0][5])*p1 - (1/k)*(nf[0]*p1+nf[1]*p2+nf[2]*p3+nf[3]*p4+nf[4]*p5+nf[5]*p6)*chi[0],
            (absp[1]+sct[1][2]+sct[1][3]+sct[1][4]+sct[1][5])*p2 - (1/k)*(nf[0]*p1+nf[1]*p2+nf[2]*p3+nf[3]*p4+nf[4]*p5+nf[5]*p6)*chi[1] - sct[0][1]*p1,
            (absp[2]+sct[2][3]+sct[2][4]+sct[2][5])*p3 - (1/k)*(nf[0]*p1+nf[1]*p2+nf[2]*p3+nf[3]*p4+nf[4]*p5+nf[5]*p6)*chi[2] - sct[0][2]*p1-sct[1][2]*p2,
            (absp[3]+sct[3][4]+sct[3][5])*p4 - (1/k)*(nf[0]*p1+nf[1]*p2+nf[2]*p3+nf[3]*p4+nf[4]*p5+nf[5]*p6)*chi[3] - sct[0][3]*p1-sct[1][3]*p2-sct[2][3]*p3,
            (absp[4]+sct[4][5])*p5 - (1/k)*(nf[0]*p1+nf[1]*p2+nf[2]*p3+nf[3]*p4+nf[4]*p5+nf[5]*p6)*chi[4] - sct[0][4]*p1-sct[1][4]*p2-sct[2][4]*p3-sct[3][4]*p4,
            (absp[5])*p6 - (1/k)*(nf[0]*p1+nf[1]*p2+nf[2]*p3+nf[3]*p4+nf[4]*p5+nf[5]*p6)*chi[5] - sct[0][5]*p1-sct[1][5]*p2-sct[2][5]*p3-sct[3][5]*p4-sct[4][5]*p5]
root = fsolve(zeroD, [1, 1, 1, 1, 1, 1])
print('roots: ', root)

print('zeroes: ', zeroD(root))