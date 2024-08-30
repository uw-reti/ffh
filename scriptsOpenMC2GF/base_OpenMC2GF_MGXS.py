# %matplotlib inline
from fileinput import filename
from operator import inv
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



def run_MGXS():
    isModelXML=False
    MC_model=None
    args=parser.parse_args()

    # Allow for base filepath other than CWD
    baseFilePath=""

    if args.filepath:
        baseFilePath = args.filepath

    baseFilePath = os.getcwd() + "/" + baseFilePath
    if baseFilePath[-1] != "/":
        baseFilePath=baseFilePath+"/"
    print("filepath: "+baseFilePath)

    # Set variables as OpenMC objects for later use, default to model.xml
    if os.path.isfile(baseFilePath+"model.xml"):
        isModel=True
        MC_model=openmc.model.from_model_xml(f"{baseFilePath}model.xml")
        MC_mats=MC_model.materials
        MC_geom=MC_model.geometry
        MC_settings=MC_model.settings
    else:
        print((f"{baseFilePath}materials.xml"))
        MC_mats=openmc.Materials.from_xml(f"{baseFilePath}materials.xml")
        MC_geom=openmc.Geometry.from_xml(f"{baseFilePath}geometry.xml", materials=MC_mats)
        MC_settings=openmc.Settings.from_xml(f"{baseFilePath}settings.xml")

    # Set particles if option argument exists
    if args.particles:
        batches = 100
        inactive = 10
        particles = int(args.particles/batches)
        MC_settings.batches = batches
        MC_settings.inactive = inactive
        MC_settings.particles = particles

    MC_settings.output={'tallies': True}

    # Create an initial uniform spatial source distribution over fissionable zones
    """     if openmc.Geometry.bounding_box:
        boundingBox=openmc.Geometry.bounding_box
        print(boundingBox)
        uniform_dist=openmc.stats.Box(boundingBox.lower_left, boundingBox.upper_right, only_fissionable=True)
    else:  """
    box_bound=args.boundingBoxSize/2
    bounds = [-box_bound, -box_bound, -box_bound, box_bound, box_bound, box_bound]
    uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
    MC_settings.source = openmc.Source(space=uniform_dist)

    # Instantiate a 6-group EnergyGroups object
    groups = mgxs.EnergyGroups()
    groups.group_edges = np.array([0., 748.5, 5531.0, 24790.0, 497900.0, 2.231e6, 20.0e6])

    # Instantiate a 100-group EnergyGroups object - not used
    # supergroup = mgxs.EnergyGroups()
    # supergroup.group_edges = np.logspace(-3, 7.3, 101)

    # Instantiate a 1-group EnergyGroups object - used for delayed properties
    one_group = mgxs.EnergyGroups()
    one_group.group_edges = np.array([groups.group_edges[0], groups.group_edges[-1]])

    delayed_groups = list(range(1,7)) # 6 delayed groups

    # Currently homogenizing over the entire universe
    # Could change this to directly convert models with multiple regions/materials
    full_universe=MC_geom.root_universe
    # Instantiate more than a few different sections
    total = mgxs.TotalXS(domain=full_universe, energy_groups=groups)
    absorption = mgxs.AbsorptionXS(domain=full_universe, energy_groups=groups)
    scattering = mgxs.ScatterXS(domain=full_universe, energy_groups=groups)
    fission = mgxs.FissionXS(domain=full_universe, energy_groups=groups, nu=False)
    nu_fission_matrix = mgxs.NuFissionMatrixXS(domain=full_universe, energy_groups=groups)
    scattering_matrix = mgxs.ScatterMatrixXS(domain=full_universe, energy_groups=groups)
    inverse_velocity = mgxs.InverseVelocity(domain=full_universe, energy_groups=groups)
    chi_prompt=mgxs.Chi(domain=full_universe, energy_groups=groups, prompt=True)
    chi=mgxs.Chi(domain=full_universe, energy_groups=groups)
    chi_delayed=mgxs.ChiDelayed(domain=full_universe, energy_groups=groups)
    nuSigmaEff=mgxs.FissionXS(domain=full_universe, energy_groups=groups, nu=True)
    diffusion_coefficient=mgxs.DiffusionCoefficient(domain=full_universe, energy_groups=groups, nu=True) # should nu be true or false? - doesn't seem to make a difference
    sigmaPow=mgxs.KappaFissionXS(domain=full_universe, energy_groups=groups)
    beta=mgxs.Beta(domain=full_universe, energy_groups=one_group, delayed_groups=delayed_groups)
    lambda1=mgxs.DecayRate(domain=full_universe, energy_groups=one_group, delayed_groups=delayed_groups)


    # Instantiate an empty Tallies object
    MC_MGXS_tallies = openmc.Tallies()

    # Add total tallies to the tallies file
    MC_MGXS_tallies += total.tallies.values()

    # Add absorption tallies to the tallies file
    MC_MGXS_tallies += absorption.tallies.values()

    # Add scattering tallies to the tallies file
    MC_MGXS_tallies += scattering.tallies.values()

    # Add fission tallies to the tallies file
    MC_MGXS_tallies += fission.tallies.values()

    # Add nu_fission_matrix tallies to the tallies file
    MC_MGXS_tallies += nu_fission_matrix.tallies.values()

    # Add scattering_matrix tallies to the tallies file
    MC_MGXS_tallies += scattering_matrix.tallies.values()

    # Add inverse velocity tallies to the tallies file
    MC_MGXS_tallies += inverse_velocity.tallies.values()

    # Add nuSigmaEff tallies to the tallies file
    MC_MGXS_tallies += nuSigmaEff.tallies.values()

    # Add chi prompt tallies to the tallies file
    MC_MGXS_tallies += chi_prompt.tallies.values()

    # Add chi tallies to the tallies file
    MC_MGXS_tallies += chi.tallies.values()

    # Add chi delayed tallies to the tallies file
    MC_MGXS_tallies += chi_delayed.tallies.values()

    # Add diffusion coefficient tallies to the tallies file
    MC_MGXS_tallies += diffusion_coefficient.tallies.values()

    # Add sigmaPow tallies to the tallies file
    MC_MGXS_tallies += sigmaPow.tallies.values()

    # Add beta tallies to the tallies file
    MC_MGXS_tallies += beta.tallies.values()

    # Add lambda tallies to the tallies file
    MC_MGXS_tallies += lambda1.tallies.values()

    # Export to "tallies.xml"
    # MC_MGXS_tallies.export_to_xml()

    finalModel=openmc.Model(geometry=MC_geom, materials=MC_mats, tallies=MC_MGXS_tallies, settings=MC_settings)
    finalModel.export_to_model_xml(path=f"{baseFilePath}MGXS_model.xml")
    # Run OpenMC
    openmc.run(path_input=f"{baseFilePath}MGXS_model.xml", cwd=baseFilePath)

    # Load the last statepoint file
    # sp = openmc.StatePoint('statepoint.60.h5') #change this with the simulation parameters 'statepoint.n.h5' n=number of batches
    sp = openmc.StatePoint(f"{baseFilePath}statepoint.{batches}.h5")

    # Load the tallies from the statepoint into each MGXS object
    total.load_from_statepoint(sp)
    absorption.load_from_statepoint(sp)
    scattering.load_from_statepoint(sp)
    fission.load_from_statepoint(sp)
    nu_fission_matrix.load_from_statepoint(sp)
    scattering_matrix.load_from_statepoint(sp)
    inverse_velocity.load_from_statepoint(sp)
    chi_prompt.load_from_statepoint(sp)
    chi.load_from_statepoint(sp)
    chi_delayed.load_from_statepoint(sp)
    nuSigmaEff.load_from_statepoint(sp)
    diffusion_coefficient.load_from_statepoint(sp)
    sigmaPow.load_from_statepoint(sp)
    beta.load_from_statepoint(sp)
    lambda1.load_from_statepoint(sp)
    # Close the statepoint file now that we're done getting info
    sp.close()

    total.print_xs()
    df = scattering.get_pandas_dataframe()
    df.head(10)

    # Originally 'excel' output type, changed to csv
    absorption.export_xs_data(filename='absorption-xs', directory=f"{baseFilePath}mgxs", format='csv')
    scattering.export_xs_data(filename='scattering-xs', directory=f"{baseFilePath}mgxs", format='csv')
    fission.export_xs_data(filename='fission-xs', directory=f"{baseFilePath}mgxs", format='csv')
    nu_fission_matrix.export_xs_data(filename='nufissionmatrix-xs', directory=f"{baseFilePath}mgxs", format='csv')
    scattering_matrix.export_xs_data(filename='scatteringmatrix-xs', directory=f"{baseFilePath}mgxs", format='csv')
    total.export_xs_data(filename='total-xs', directory=f"{baseFilePath}mgxs", format='csv')
    inverse_velocity.export_xs_data(filename='inverse-velocity', directory=f"{baseFilePath}mgxs", format='csv')
    chi_prompt.export_xs_data(filename='chi-prompt', directory=f"{baseFilePath}mgxs", format='csv')
    chi.export_xs_data(filename='chi', directory=f"{baseFilePath}mgxs", format='csv')
    chi_delayed.export_xs_data(filename='chi-delayed', directory=f"{baseFilePath}mgxs", format='csv')
    nuSigmaEff.export_xs_data(filename='nu-SigmaEff', directory=f"{baseFilePath}mgxs", format='csv')
    diffusion_coefficient.export_xs_data(filename='diffusion-coefficient', directory=f"{baseFilePath}mgxs", format='csv')
    sigmaPow.export_xs_data(filename='sigmaPow', directory=f"{baseFilePath}mgxs", format='csv')
    beta.export_xs_data(filename='betaone', directory=f"{baseFilePath}mgxs", format='csv')
    lambda1.export_xs_data(filename='lambdaone', directory=f"{baseFilePath}mgxs", format='csv')

    total.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    absorption.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    scattering.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    fission.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    nu_fission_matrix.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    scattering_matrix.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    inverse_velocity.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    chi_prompt.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    chi.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    chi_delayed.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    nuSigmaEff.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    diffusion_coefficient.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    sigmaPow.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    beta.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)
    lambda1.build_hdf5_store(filename='mgxs', directory=f"{baseFilePath}mgxs", append=True)


def CSV_to_GF():
    args=parser.parse_args()
    fileName = "nuclearData.txt"
    if args.csv_filepath:
        baseFilePath = args.csv_filepath
        if baseFilePath[-1] != "/":
            baseFilePath=baseFilePath+"/"

    elif args.filepath:
        baseFilePath = args.filepath
        if baseFilePath[-1] != "/":
            baseFilePath=baseFilePath+"/"
    base_filepath_csvs = baseFilePath+"mgxs/"
    resultsFile = "results"

    # Begin writing nuclearData file	
    f = open(f"{baseFilePath}{fileName}","w")

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

    if args.copyToAllRequiredFiles:
        dataFileList=["Ksalt_pertdens", "Ksalt_perttemp", 
                      "nuclearData", "nuclearDataFuelTemp",
                      "nuclearDataRhoCool", "nuclearDataTCool"]
        os.mkdir(f"{baseFilePath}/GF_neutroRegion_data")
        originalFile = f"{baseFilePath}{fileName}"
        for dataFileName in dataFileList:
            shutil.copyfile(originalFile, f"{baseFilePath}/GF_neutroRegion_data/{dataFileName}")


    print("DONE")





# Create argparser for some non-required variables such as filepath
parser = argparse.ArgumentParser(prog="OpenMC2GF_MGXS", 
                                 description="Conversion of OpenMC Model into Cross-sections usable in a GeN-Foam simulation", 
                                 epilog="V0.1_base")
parser.add_argument('-f', '--filepath', help='filepath where OpenMC xml files are stored (relative from run directory)', default="./")
parser.add_argument('-p', '--particles', help='A separate option for setting the total number of active particles for XS generation', default=1000000)
parser.add_argument('-b', '--boundingBoxSize', help="Side length of bounding box if it doesn't exist in the geometry (centered at origin)", default=1000)
parser.add_argument('--csv_filepath', help='filepath where OpenMC csv files are stored (mgxs directory)')
parser.add_argument('-m', '--mgxs_run', action="store_true")
parser.add_argument('-c', '--CSV2GF_run', action="store_true")
parser.add_argument('-mc', '--combinedRun', action="store_true", help="Run OpenMC MGXS based on xml files, then run CSV2GF Conversion Script")
parser.add_argument('-copyToAllRequiredFiles', help="Copy the produced GF data file to all files required for a transient run", default=True)


args=parser.parse_args()
if args.combinedRun:
    run_MGXS()
    CSV_to_GF()
elif args.mgxs_run:
    run_MGXS()
elif args.CSV2GF_run:
    CSV_to_GF()
else:
    print("Choose either -m (mgxs run), -c (CSV2GF run), or -mc (combined run) option")
