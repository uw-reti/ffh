from datetime import datetime
import numpy as np
import pandas as pd
import sys
import argparse

parser = argparse.ArgumentParser(prog="HO_OpenMC2GF", 
                                 description="Conversion of CSVs from OpenMC into Cross-sections usable in a GeN-Foam simulation", 
                                 epilog="convert script only")
parser.add_argument('-f', '--csv_filepath', help='filepath where OpenMC csv files are stored (mgxs directory)')
args=parser.parse_args()


fileName = "./nuclearData.txt"
if args.filepath:
    baseFilePath = args.filepath
    if baseFilePath[-1] != "/":
        baseFilePath=baseFilePath+"/"
base_filepath_csvs = baseFilePath+"mgxs/"
resultsFile = "results"

# Begin writing nuclearData file	
f = open(fileName,"w")

# Opening information
now = datetime.now()

dt_string = now.strftime("%d/%m/%Y %H:%M:%S")

f.write("/*\n crossSection dictionary\n Generated by OpenMC2FoamXS \n %s\n OpenMC results file: %s\n" % (dt_string, resultsFile))

f.write("\neffective delayed neutron fraction, prompt neutron spectrum\n*/\n")

f.write("\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    location    constant;\n    object      %s;\n}\n" % resultsFile);

# Initial simulation parameters
f.write("\nfastNeutrons %s;\n" % '         true' )
f.write("\nadjustDiscFactors %s;\n" % '    false' )
f.write("\nuseGivenDiscFactors %s;\n" % '  false' )

# prompt generation and delayed parameters
# Kinetics properties have to be defined at the beginning of the file
f.write("\npromptGenerationTime %.2e;\n" % 6.70e-07) #Calculate a new value for this

energyGroups = 6
precGroups = 6

beta_eff      = [0, 0, 0, 0, 0, 0]  
decay_const   = [0, 0, 0, 0, 0, 0] 


# Beta (note this is not beta_eff yet, it is the default beta from OpenMC)
df = pd.read_csv(f"{base_filepath_csvs}betaone.csv")
i=0
for row in df['mean']:
    beta_eff[i] = row
    i += 1 

# Lambda
df = pd.read_csv(f"{base_filepath_csvs}lambdaone.csv")
i=0
for row in df['mean']:
    decay_const[i] = row
    i += 1 

# Initial print beta 
f.write("\nBeta (%s  )" % str(' '.join(['{:.6e}'.format(x) for x in beta_eff])) + ";\n")

# Initial print decay constants for precursors
f.write("\nlambda (%s  )" % str(' '.join(['{:.6e}'.format(x) for x in decay_const])) + ";\n")

# Feedback coefficients
f.write("\nfeedbackCoeffFastDoppler %.9f;\n" % 2.2808e-05) #Calculate a new value for this
f.write("\nfeedbackCoeffTFuel %i;\n" % 0)
f.write("\nfeedbackCoeffTClad %i;\n" % 0)
f.write("\nfeedbackCoeffTCool %i;\n" % 0)
# f.write("\nfeedbackCoeffRhoCool %.5e;\n" % 4.55597e-05) #Calculate a new value for this
f.write("\nfeedbackCoeffRhoCool %.5e;\n" % 14.2138e-05) #Calculate a new value for this
f.write("\nfeedbackCoeffTStruct %i;\n" % 0)
f.write("\nabsoluteDrivelineExpansionCoeff %i;\n" % 0)
f.write("\ncontrolRodReactivityMap %s; \n" % "( ( 0.1 -0.01 ) ( 0 0 ) ( -0.1 0.01 ) )")

# Write Number of energy & precursor groups
f.write("\ninitPrecursorsLiquidFuel %s;\n\n" % 'true' )
f.write("\n energyGroups %i ;\n" % energyGroups)
f.write("\n precGroups %i ;\n" % precGroups)

# Write the main nuclear data

# Add array of region names here, may need another for loop to cycle through the different regions
OF_NAME = ["hx","intermed","main_fd","pump"]  # The name of the GeN-Foam mesh region this cross section is for

# Data for GeN-Foam
fuel_fraction = 1.000000e+00        # The volumetric fuel fraction
inv_velocity  = [0, 0, 0, 0, 0, 0]  # m/s
diffCoeff     = [0, 0, 0, 0, 0, 0]  # 
nusigmaf      = [0, 0, 0, 0, 0, 0]  # 
sigmaFissPow  = [0, 0, 0, 0, 0, 0]  # 
scattering_P0 = [[0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0],[0, 0, 0, 0, 0, 0]]   # 
scattering    = [0, 0, 0, 0, 0, 0]  # used to calculate sigmaDisapp
scattering_ODS= [0, 0, 0, 0, 0, 0]  # used to calculate sigmaDisapp
sigmaDisapp   = [0, 0, 0, 0, 0, 0]  # 
chi_prompt    = [0, 0, 0, 0, 0, 0]  # 
chi_delayed   = [0, 0, 0, 0, 0, 0]  # 
beta_eff      = [0, 0, 0, 0, 0, 0]  #
decay_const   = [0, 0, 0, 0, 0, 0]  #
generic_6     = [1, 1, 1, 1, 1, 1]  # Used for any arrays that are all 1 (discFactor, integralFlux)
absorption    = [0, 0, 0, 0, 0, 0]  # 
fission       = [0, 0, 0, 0, 0, 0]  # 
total         = [0, 0, 0, 0, 0, 0]  # used to calculate sigmaDisapp
cmm           = 100                 # used to convert from cm^-1 to m^-1
ev2j         = 1.602176487e-19      # used to convert eV to J

## Example in case data is different for each zone
# fuel_fraction_zone1 = 1.000000e-01
# fuel_fraction_zone2 = 1.000000e-01
# fuel_fraction_zone3 = 1.000000e-01
# fuel_fraction_all = [fuel_fraction_zone1, fuel_fraction_zone2, fuel_fraction_zone3]
# for zone in (OF_NAME):
#   for i in range(OF_NAME): 
#       f.write("  fuelFraction %s ; \n" % "{:.6e}".format(fuel_fraction[i]))

# Fill variable lists with data from OpenMC

# inverse velocity
df = pd.read_csv(f"{base_filepath_csvs}inverse-velocity.csv")
i=0
for row in df['mean']:
    inv_velocity[i] = row*cmm
    i += 1 

# Diffusion Coefficient
df = pd.read_csv(f"{base_filepath_csvs}diffusion-coefficient.csv")
i=0
for row in df['mean']:
    diffCoeff[i] = row/cmm
    i += 1

# nusigmaf
df = pd.read_csv(f"{base_filepath_csvs}nu-SigmaEff.csv")
i=0
for row in df['mean']:
    nusigmaf[i] = row*cmm
    i += 1 
    
# sigmaPow
df = pd.read_csv(f"{base_filepath_csvs}sigmaPow.csv")
i=0
for row in df['mean']:
    sigmaFissPow[i] = row*cmm*ev2j
    i += 1

# sigmaDisapp (TOTAL - SCATTERING) with transport correction accounted for
df = pd.read_csv(f"{base_filepath_csvs}total-xs.csv")
i=0
for row in df['mean']:
    total[i] = row*cmm
    i += 1 

# Scattering Matrix
df = pd.read_csv(f"{base_filepath_csvs}scatteringmatrix-xs.csv")
i=0
j=0
for row in df['mean']:
    scattering_P0[i][j] = row*cmm
    j += 1
    if j == 6:
        i += 1 
        j = 0

# Scattering XS (array)
df = pd.read_csv(f"{base_filepath_csvs}scattering-xs.csv")
i=0
for row in df['mean']:
    scattering[i] = row*cmm
    i += 1 

# Define off-diagonal sum (ODS) of scattering matrix terms for each group
scattering_ODS[0] = scattering_P0[0][1] + scattering_P0[0][2] + scattering_P0[0][3] + scattering_P0[0][4] + scattering_P0[0][5]
scattering_ODS[1] = scattering_P0[1][2] + scattering_P0[1][3] + scattering_P0[1][4] + scattering_P0[1][5]
scattering_ODS[2] = scattering_P0[2][3] + scattering_P0[2][4] + scattering_P0[2][5]
scattering_ODS[3] = scattering_P0[3][4] + scattering_P0[3][5]
scattering_ODS[4] = scattering_P0[4][5]
scattering_ODS[5] = 0

# Subtract ODS from scattering for Non-Transport-Corrected Scattering (NTC)? Still not sure exactly how this works
scattering_NTC = np.subtract(scattering, scattering_ODS)

# Finally write sigmaDisapp
for i in range(len(sigmaDisapp)):
    sigmaDisapp[i] = total[i] - scattering_NTC[i]

# chi_prompt
df = pd.read_csv(f"{base_filepath_csvs}chi-prompt.csv")
i=0
for row in df['mean']:
    chi_prompt[i] = row
    i += 1 
print('chi_prompt: ', chi_prompt)

# chi_delayed
df = pd.read_csv(f"{base_filepath_csvs}chi-delayed.csv")
i=0
for row in df['mean']:
    chi_delayed[i] = row
    i += 1 

# Scattering Matrix
df = pd.read_csv(f"{base_filepath_csvs}scatteringmatrix-xs.csv")
i=0
j=0
for row in df['mean']:
    scattering_P0[i][j] = row*cmm
    j += 1
    if j == 6:
        i += 1 
        j = 0

# Beta (note this is not beta_eff yet, it is the default beta from OpenMC)
df = pd.read_csv(f"{base_filepath_csvs}betaone.csv")
i=0
for row in df['mean']:
    beta_eff[i] = row
    i += 1 

# Lambda
df = pd.read_csv(f"{base_filepath_csvs}lambdaone.csv")
i=0
for row in df['mean']:
    decay_const[i] = row
    i += 1 

# Start writing nuclearData
f.write("\nzones \n ( \n")

for zone in (OF_NAME):
    f.write("\n%s \n{ \n" % zone)

    f.write(" fuelFraction %s ; \n" % "{:.6e}".format(fuel_fraction) + "\n")

    f.write(" IV nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in inv_velocity])) + "  );\n")

    # Print groupwise diffusion coefficients
    f.write("\n D nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in diffCoeff])) + "  );\n")

    # Print groupwise nu sigma f
    f.write("\n nuSigmaEff nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in nusigmaf])) + "  );\n")

    # Print groupwise fission powers
    f.write("\n sigmaPow nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in sigmaFissPow])) + "  );\n")

    # Print scattering production matrix - P0 only
    f.write("\n scatteringMatrixP0 %i %i (\n" % (energyGroups, energyGroups))

    for i in range (energyGroups):
        f.write(" (")
        f.write(str(' '.join(['{:.6e}'.format(x) for x in scattering_P0[i][:]])) + " )\n")
    f.write(" );\n")

    # Print sigma absorption

    f.write("\n sigmaDisapp nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in sigmaDisapp])) + "  );\n") # sigmaAbs 0.85 keff

    f.write("\n chiPrompt nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in chi_prompt])) + "  );\n")

    f.write("\n chiDelayed nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in chi_delayed])) + "  );\n")

    # Print beta 
    f.write("\n Beta nonuniform List<scalar> %i (" % precGroups + str(' '.join(['{:.6e}'.format(x) for x in beta_eff])) + "  );\n")

    # Print decay constants for precursors
    f.write("\n lambda nonuniform List<scalar> %i (" % precGroups + str(' '.join(['{:.6e}'.format(x) for x in decay_const])) + "  );\n")

    # Disc factors - use generic unity matrix
    f.write("\n discFactor nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in generic_6])) + "  );\n")

    # Print integral fluxes 
    f.write("\n integralFlux nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in generic_6])) + "  );\n")

    f.write("\n }\n")

f.write("\n); \n")
f.close()

print("DONE")