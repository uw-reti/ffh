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



def run_MGXS(perturbation=""):
    MC_model=None
    args=parser.parse_args()

    # Allow for base filepath other than CWD, set base file path
    if args.filepath:
        baseFilePath = args.filepath
    else:
        baseFilePath=""

    baseFilePath = os.getcwd() + "/" + baseFilePath
    if baseFilePath[-1] != "/":
        baseFilePath=baseFilePath+"/"

    # Set variables as OpenMC objects for later use, default to model.xml
    if os.path.isfile(baseFilePath+"model.xml"):
        MC_model=openmc.model.from_model_xml(f"{baseFilePath}model.xml")
        MC_mats=MC_model.materials
        MC_geom=MC_model.geometry
        MC_settings=MC_model.settings
    else:
        MC_mats=openmc.Materials.from_xml(f"{baseFilePath}materials.xml")
        MC_geom=openmc.Geometry.from_xml(f"{baseFilePath}geometry.xml", materials=MC_mats)
        MC_settings=openmc.Settings.from_xml(f"{baseFilePath}settings.xml")
        

    # Set up perturbations if it is run
    if perturbation:
        global original_value
        global perturbed_value
        for i in range(0,np.size(MC_mats)):
            if int(MC_mats[i].id) == int(args.mat_id):
                perturbed_mat=MC_mats[i]
                mat_list_position=i
        if perturbation == "rho":
            rho_mat=perturbed_mat.density
            original_value = rho_mat
            MC_mats[mat_list_position].set_density('g/cm3',rho_mat*1.1)
            perturbed_value = rho_mat*1.1
            # TODO create a way to output the perturbed state files

        elif perturbation == "temp":
            t_mat=perturbed_mat.temperature
            original_value = t_mat
            MC_mats[mat_list_position].temperature = t_mat*1.2
            perturbed_value = t_mat*1.2
            # TODO create a way to output the perturbed state files
    else:
        perturbation=""

    # Set particles if option argument exists
    if args.particles:
        batches = 100
        inactive = 10
        particles = int(args.particles/batches)
        MC_settings.batches = batches
        MC_settings.inactive = inactive
        MC_settings.particles = particles

    MC_settings.output={'tallies': True}
    MC_settings.temperature = {'method' : 'interpolation'}

    # Create an initial uniform spatial source distribution over fissionable zones
    """     if openmc.Geometry.bounding_box:
        boundingBox=openmc.Geometry.bounding_box
        print(boundingBox)
        uniform_dist=openmc.stats.Box(boundingBox.lower_left, boundingBox.upper_right, only_fissionable=True)
    else:  """
    box_bound=float(args.boundingBoxSize)/2
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
    diffusion_coefficient=mgxs.DiffusionCoefficient(domain=full_universe, energy_groups=groups, nu=True) # TODO Check effect of nu true/false
    sigmaPow=mgxs.KappaFissionXS(domain=full_universe, energy_groups=groups)
    beta=mgxs.Beta(domain=full_universe, energy_groups=one_group, delayed_groups=delayed_groups)
    lambda1=mgxs.DecayRate(domain=full_universe, energy_groups=one_group, delayed_groups=delayed_groups)

    XS_list=[total,absorption,scattering,fission,nu_fission_matrix,scattering_matrix,inverse_velocity,
             chi_prompt,chi,chi_delayed,nuSigmaEff,diffusion_coefficient,sigmaPow,beta,lambda1]

    # Instantiate an empty Tallies object
    MC_MGXS_tallies = openmc.Tallies()

    # Add all tallies to the tallies file
    for XS_objs in XS_list:
        MC_MGXS_tallies += XS_objs.tallies.values()

    finalModel=openmc.Model(geometry=MC_geom, materials=MC_mats, tallies=MC_MGXS_tallies, settings=MC_settings)
    finalModel.export_to_model_xml(path=f"{baseFilePath}MGXS_model.xml")
    # Run OpenMC
    openmc.run(path_input=f"{baseFilePath}MGXS_model.xml", cwd=baseFilePath)

    # Load the last statepoint file
    sp = openmc.StatePoint(f"{baseFilePath}statepoint.{batches}.h5")

    # Load the tallies from the statepoint into each MGXS object
    for XS_objs in XS_list:
        XS_objs.load_from_statepoint(sp)
    sp.close()

    total.print_xs()
    df = scattering.get_pandas_dataframe()
    df.head(10)

    # Export XS data as csv files for use in the CSV_TO_GF function
    csv_filenames=['total-xs', 'absorption-xs', 'scattering-xs', 'fission-xs', 'nufissionmatrix-xs', 'scatteringmatrix-xs', 'inverse-velocity', 
               'chi-prompt', 'chi', 'chi-delayed', 'nu-SigmaEff', 'diffusion-coefficient', 'sigmaPow', 'betaone', 'lambdaone']
    for i in range(0,np.size(csv_filenames)):
        XS_list[i].export_xs_data(filename=f"{perturbation}{csv_filenames[i]}", directory=f"{baseFilePath}mgxs", format='csv')
    for XS_objs in XS_list:
        XS_objs.build_hdf5_store(filename=f'{perturbation}mgxs', directory=f"{baseFilePath}mgxs", append=True)
    if perturbation:
        os.system(f"mkdir {baseFilePath}{perturbation}")
        os.system(f'cp {baseFilePath}MGXS_model.xml {baseFilePath}{perturbation}/MGXS_model.xml') 
        os.system(f'cp -r {baseFilePath}mgxs {baseFilePath}{perturbation}/')

# Write the standard header for the files
def writeHeader():
    fnWriteString = ""

    # make header:
    resultsFile = "results"
    now = datetime.now()
    dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
    fnWriteString += "/*\n crossSection dictionary\n Generated by OpenMC2FoamXS \n %s\n OpenMC results file: %s\n" % (dt_string, resultsFile)
    fnWriteString+="\neffective delayed neutron fraction, prompt neutron spectrum\n*/\n"
    fnWriteString+="\nFoamFile\n{\n    version     2.0;\n    format      ascii;\n    class       dictionary;\n    location    constant;\n    object      %s;\n}\n" % resultsFile
    
    return fnWriteString


def CSV_to_GF(perturbation=""):
    args=parser.parse_args()
    fileName = "nuclearData.txt"
    baseFilePath=""
    # Deal with potential filepath arguments
    if args.csv_filepath:
        baseFilePath = args.csv_filepath
    elif args.filepath:
        baseFilePath = args.filepath
    if baseFilePath[-1] != "/":
        baseFilePath=baseFilePath+"/"
    base_filepath_csvs = baseFilePath+"mgxs/"
    resultsFile = "results"

    # Begin writing nuclearData file
    if perturbation:
        f = open(f"{baseFilePath}{perturbation}/{fileName}","w")
    else: 
        f = open(f"{baseFilePath}{fileName}","w")

    fileString=""

    # File header
    fileString+=writeHeader()

    # Initial simulation parameters
    fileString+="\nfastNeutrons %s;\n" % '         true' 
    fileString+="\nadjustDiscFactors %s;\n" % '    false' 
    fileString+="\nuseGivenDiscFactors %s;\n" % '  false'

    if perturbation:
        fileString += export_perturbed(perturbation)

    # prompt generation and delayed parameters
    # Kinetics properties have to be defined at the beginning of the file
    fileString+="\npromptGenerationTime TODO;\n" # Default: 6.70e-07 # TODO Calculate a new value for this for each fuel

    energyGroups = 6
    precGroups = 6

    beta_eff      = [0, 0, 0, 0, 0, 0]  
    decay_const   = [0, 0, 0, 0, 0, 0] 


    # Beta (note this is not beta_eff yet, it is the default beta from OpenMC)
    beta_eff = pd.read_csv(f"{base_filepath_csvs}betaone.csv")['mean']

    # Lambda
    decay_const = pd.read_csv(f"{base_filepath_csvs}lambdaone.csv")['mean']

    # Initial print beta 
    fileString+="\nBeta (%s  )" % str(' '.join(['{:.6e}'.format(x) for x in beta_eff])) + ";\n"

    # Initial print decay constants for precursors
    fileString+="\nlambda (%s  )" % str(' '.join(['{:.6e}'.format(x) for x in decay_const])) + ";\n"

    # Feedback coefficients
    # Many of these require new values for each fuel
    fileString+="\nfeedbackCoeffFastDoppler TODO;\n" # default: 2.2808e-05
    fileString+="\nfeedbackCoeffTFuel TODO;\n" # default: 0
    fileString+="\nfeedbackCoeffTClad TODO;\n" # default: 0
    fileString+="\nfeedbackCoeffTCool TODO;\n" # default: 0
    fileString+="\nfeedbackCoeffRhoCool TODO;\n" # default: 4.55597e-05
    fileString+="\nfeedbackCoeffTStruct TODO;\n" # default: 0
    fileString+="\nabsoluteDrivelineExpansionCoeff TODO;\n" # default: 0
    fileString+="\ncontrolRodReactivityMap ( ( 0.1 -0.01 ) ( 0 0 ) ( -0.1 0.01 ) ); \n" 

    # Write Number of energy & precursor groups
    fileString+="\ninitPrecursorsLiquidFuel %s;\n\n" % 'true' 
    fileString+="\n energyGroups %i ;\n" % energyGroups
    fileString+="\n precGroups %i ;\n" % precGroups

    # Write the main nuclear data
    f.write(fileString)
    fileString=""

    # Add array of region names here, may need another for loop to cycle through the different regions
    OF_NAME = ["hx","intermed","main_fd","pump"]  # The name of the GeN-Foam mesh region this cross section is for

    # Data required for GeN-Foam
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
    #       fileString+="  fuelFraction %s ; \n" % "{:.6e}".format(fuel_fraction[i])

    # Fill variable lists with data from OpenMC

    # inverse velocity
    inv_velocity = pd.read_csv(f"{base_filepath_csvs}inverse-velocity.csv")['mean']*cmm

    # Diffusion Coefficient
    diffCoeff = pd.read_csv(f"{base_filepath_csvs}diffusion-coefficient.csv")['mean']/cmm

    # nusigmaf
    nusigmaf = pd.read_csv(f"{base_filepath_csvs}nu-SigmaEff.csv")['mean']*cmm
        
    # sigmaPow
    sigmaFissPow = pd.read_csv(f"{base_filepath_csvs}sigmaPow.csv")['mean']*cmm*ev2j

    # sigmaDisapp (TOTAL - SCATTERING) with transport correction accounted for
    total = pd.read_csv(f"{base_filepath_csvs}total-xs.csv")['mean']*cmm

    # Scattering Matrix
    scattering_P0_temp = pd.read_csv(f"{base_filepath_csvs}scatteringmatrix-xs.csv")['mean']*cmm
    scattering_P0=np.zeros([6,6])
    for i in range(0,6): 
        for j in range(0,6):
            scattering_P0[i][j] = scattering_P0_temp[i*6+j]

    # Scattering XS (array)
    scattering = pd.read_csv(f"{base_filepath_csvs}scattering-xs.csv")['mean']*cmm

    # Define off-diagonal sum (ODS) of scattering matrix terms for each group
    for i in range(0,5):
        scattering_ODS[i] = sum(scattering_P0[i][(i+1):5])
    scattering_ODS[5]=0

    """ scattering_ODS[0] = scattering_P0[0][1] + scattering_P0[0][2] + scattering_P0[0][3] + scattering_P0[0][4] + scattering_P0[0][5]
    scattering_ODS[1] = scattering_P0[1][2] + scattering_P0[1][3] + scattering_P0[1][4] + scattering_P0[1][5]
    scattering_ODS[2] = scattering_P0[2][3] + scattering_P0[2][4] + scattering_P0[2][5]
    scattering_ODS[3] = scattering_P0[3][4] + scattering_P0[3][5]
    scattering_ODS[4] = scattering_P0[4][5]
    scattering_ODS[5] = 0 """

    # Subtract ODS from scattering for Non-Transport-Corrected Scattering (NTC)?
    # TODO check this calculation/definition
    scattering_NTC = np.subtract(scattering, scattering_ODS)

    # Finally write sigmaDisapp
    for i in range(len(sigmaDisapp)):
        sigmaDisapp[i] = total[i] - scattering_NTC[i]

    # chi_prompt
    chi_prompt = pd.read_csv(f"{base_filepath_csvs}chi-prompt.csv")['mean']

    # chi_delayed
    chi_delayed = pd.read_csv(f"{base_filepath_csvs}chi-delayed.csv")['mean']

    # Beta (note this is not beta_eff yet, it is the default beta from OpenMC)
    beta_eff = pd.read_csv(f"{base_filepath_csvs}betaone.csv")['mean']

    # Lambda
    decay_const = pd.read_csv(f"{base_filepath_csvs}lambdaone.csv")['mean']

    # Start writing nuclearData with zones
    fileString="\nzones \n ( \n"

    for zone in (OF_NAME):
        fileString+="\n%s \n{ \n" % zone

        fileString+=" fuelFraction %s ; \n" % "{:.6e}".format(fuel_fraction) + "\n"
        fileString+=" IV nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in inv_velocity])) + "  );\n"

        # Groupwise diffusion coefficientsm nuSigma and fission powers
        fileString+="\n D nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in diffCoeff])) + "  );\n"
        fileString+="\n nuSigmaEff nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in nusigmaf])) + "  );\n"
        fileString+="\n sigmaPow nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in sigmaFissPow])) + "  );\n"

        # Scattering production matrix - P0 only
        fileString+="\n scatteringMatrixP0 %i %i (\n" % (energyGroups, energyGroups)

        # Energy groups
        for i in range(energyGroups):
            fileString+=" ("
            fileString+=str(' '.join(['{:.6e}'.format(x) for x in scattering_P0[i][:]])) + " )\n"
        fileString+=" );\n"

        # Other groupwise data
        fileString+="\n sigmaDisapp nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in sigmaDisapp])) + "  );\n"
        fileString+="\n chiPrompt nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in chi_prompt])) + "  );\n"
        fileString+="\n chiDelayed nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in chi_delayed])) + "  );\n"
        # Beta 
        fileString+="\n Beta nonuniform List<scalar> %i (" % precGroups + str(' '.join(['{:.6e}'.format(x) for x in beta_eff])) + "  );\n"
        # Decay constants for precursors
        fileString+="\n lambda nonuniform List<scalar> %i (" % precGroups + str(' '.join(['{:.6e}'.format(x) for x in decay_const])) + "  );\n"
        # Disc factors - use generic unity matrix
        fileString+="\n discFactor nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in generic_6])) + "  );\n"
        # Print integral fluxes 
        fileString+="\n integralFlux nonuniform List<scalar> %i (" % energyGroups + str(' '.join(['{:.6e}'.format(x) for x in generic_6])) + "  );\n"
        fileString+="\n }\n"

    fileString+="\n); \n"
    f.write(fileString)
    f.close()

    """ # TODO delete
    if args.copyToAllRequiredFiles:
        dataFileList=["Ksalt_pertdens", "Ksalt_perttemp", 
                      "nuclearData", "nuclearDataFuelTemp",
                      "nuclearDataRhoCool", "nuclearDataTCool"]
        os.mkdir(f"{baseFilePath}/GF_neutroRegion_data")
        originalFile = f"{baseFilePath}{fileName}"
        for dataFileName in dataFileList:
            shutil.copyfile(originalFile, f"{baseFilePath}/GF_neutroRegion_data/{dataFileName}") """

# Export the perturbed value files such as Ksalt_pertdens and nuclearDataRhoCool
def export_perturbed(perturbation=""):
    fnWriteString=""
    if perturbation == 'rho':
        fnWriteString += f"rhoCoolRef {original_value} ;\n"
        fnWriteString += f"rhoCoolPerturbed {perturbed_value};\n"
    elif perturbation == 'temp':
        fnWriteString += f"TfuelRef {original_value} ;\n" 
        fnWriteString += f"TfuelPerturbed {perturbed_value} ; \n"
    else:
        raise Exception("Unhandled perturbation")
    return fnWriteString
    



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
# parser.add_argument('-copyToAllRequiredFiles', help="Copy the produced GF data file to all files required for a transient run", default=True)
parser.add_argument('-i', '--mat_id', help='The material ID used for temperature and density perturbation, a non-zero argument here will enable perturbation', default=False)

args=parser.parse_args()


if args.mgxs_run:
    run_MGXS()
if args.CSV2GF_run:
    CSV_to_GF()
else:
    print("Choose -m (mgxs run) or -c (CSV2GF run) option (or both)")
    # print("To enable density and temperature perturbations, use the -mat or --material flag")

if args.mat_id:
    perts_list = ['rho', 'temp']
    for perts in perts_list:
        if args.mgxs_run:
            run_MGXS(perturbation=perts)
        if args.CSV2GF_run:
            CSV_to_GF(perturbation=perts)
            # export_perturbed(perturbation=perts)
    
print("DONE")
