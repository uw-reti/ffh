# %matplotlib inline
from fileinput import filename
from operator import inv
import numpy as np
import matplotlib.pyplot as plt
import os

import openmc
import sys
import openmc.mgxs as mgxs

# Create and mix materials
# fuel
fuel = openmc.Material(1,"fuel", temperature=873.15)
fuel.set_density('g/cm3', 1.62134)
fuel.add_nuclide('Np237', 0.01534, 'wo')
fuel.add_nuclide('Am241', 0.158667, 'wo')
fuel.add_nuclide('Am242', 0.0010575, 'wo')
fuel.add_nuclide('Am243', 0.110917, 'wo')
fuel.add_nuclide('Cm242', 0.0000847333, 'wo')
fuel.add_nuclide('Cm243', 0.00068, 'wo')
fuel.add_nuclide('Cm244', 0.06545, 'wo')
fuel.add_nuclide('Cm245', 0.0714333, 'wo')
fuel.add_nuclide('Cm246', 0.0012405, 'wo')

# By at%
# fuel.add_nuclide('Np237', 0.01534)
# fuel.add_nuclide('Am241', 0.158667)
# fuel.add_nuclide('Am242', 0.0010575)
# fuel.add_nuclide('Am243', 0.110917)
# fuel.add_nuclide('Cm242', 0.0000847333)
# fuel.add_nuclide('Cm243', 0.00068)
# fuel.add_nuclide('Cm244', 0.06545)
# fuel.add_nuclide('Cm245', 0.0714333)
# fuel.add_nuclide('Cm246', 0.0012405)

li_salt = openmc.Material(2, "li_salt", temperature=873.15) # these temps don't matter
li_salt.add_nuclide('Li7', 0.9241)
li_salt.add_nuclide('Li6', 0.0759)
# liRat=float(sys.argv[1])
# li_salt.add_nuclide('Li7', (1-liRat))
# li_salt.add_nuclide('Li6', liRat)

chlor = openmc.Material(3, "chlor", temperature=873.15)
# chlor.add_element('Cl', 1)
chlor.add_nuclide('Cl35', 0.05)
chlor.add_nuclide('Cl37', 0.95)

k_salt = openmc.Material(4, "k_salt", temperature=873.15)
k_salt.add_nuclide('K39', 0.932581)
k_salt.add_nuclide('K41', 0.067302)

be_salt = openmc.Material(14, "be_salt", temperature=873.15)
be_salt.add_element('Be', 1)

f_salt = openmc.Material(114, "f_salt", temperature=873.15)
f_salt.add_element('F', 1)

pb_salt = openmc.Material(1114, "pb_salt", temperature=900)
pb_salt.add_element('Pb', 1)

# For defining different values of Ratio
#Ratio = float(sys.argv[1])

# Constant Ratio - this refers to amount of MA relative to Li and K
Ratio = 0.2

# Eutectic is ~45 KCl 55 LiCl
MA = 1

Li = (MA/Ratio) * .55
K = (MA/Ratio) * .45
Cl = MA*3.5 + Li + K
tot = Li + K + MA + Cl
mixList = [MA/tot, Li/tot, K/tot, Cl/tot]
mixList[3] = mixList[3] + (1-sum(mixList)) # This step simply adds the amount necessary to Cl to get 100% in the mixList
fuel_salt = openmc.Material.mix_materials([fuel, li_salt, k_salt, chlor], mixList, 'ao') # default is 'ao'

# Define fuel salt temperature for density calculation
FStemp = 900 # Kelvin
FStempC = FStemp - 273.15

# Set overall fuel salt density

# Perturbed density value (arbitrary) 
# nominal density multiplied by 0.828848 (from normal GeN-Foam ULOF case)
fuel_salt.set_density('g/cm3', 2.185414447)

# New density Calculation
# From Kumar/Anderson paper
LiKCL_eq = 2.0286-(5.2676*10**-4)*FStempC

# Correction factor from Kumar/Anderson paper. Accounts for Trivalent elements in the mix
rho_cf = 0.82601

# From molten salt references paper https://moltensalt.org/references/static/downloads/pdf/element-salt-densities.pdf
LiCl_eq = 1.502 - 0.000432 * FStempC
KCl_eq = 1.527 - 0.000583 * FStempC
AmCl3_eq = rho_cf*6.6197 - 0.001522* (FStempC + 273.15)

# find weight percents based on "Ratio"
M_likcl = (6.941*0.6+39.10*0.4+37) * (1-Ratio)
M_amCl3 = (37*3+241) * Ratio

wp_likcl = M_likcl/(M_likcl+M_amCl3)
wp_amCl3 = M_amCl3/(M_likcl+M_amCl3)

# print(wp_amCl3)

# Calculate density from above information
rho_salt = 1/((wp_likcl/LiKCL_eq) + (wp_amCl3/AmCl3_eq)) # final density in g/cm^3

print('2.185414447')

#fuel_salt.set_density('g/cm3', rho_salt)

# for constant density at perturbed temp (900K density)
#fuel_salt.set_density('g/cm3', 2.6366890513280508)

# Set overall fuel salt temperature - nominal
fuel_salt.temperature = 900 # Original Value = 900 K



# Add materials 
mats = openmc.Materials()
mats += [fuel, li_salt, k_salt, chlor, fuel_salt]
mats.export_to_xml()

# Instantiate boundary Planes (default bounds were 0.63)
min_x = openmc.XPlane(boundary_type='reflective', x0=-100)
max_x = openmc.XPlane(boundary_type='reflective', x0=100)
min_y = openmc.YPlane(boundary_type='reflective', y0=-100)
max_y = openmc.YPlane(boundary_type='reflective', y0=100)
min_z = openmc.ZPlane(boundary_type='reflective', z0=-100)
max_z = openmc.ZPlane(boundary_type='reflective', z0=100)

# Instantiate a Cell
cell = openmc.Cell(cell_id=1, name='cell')

# Register bounding Surfaces with the Cell
cell.region = +min_x & -max_x & +min_y & -max_y & +min_z & -max_z

# Fill the Cell with the Material
cell.fill = fuel_salt

# Create root universe
root_universe = openmc.Universe(name='root universe', cells=[cell])

# Create Geometry and set root Universe
openmc_geometry = openmc.Geometry(root_universe)

# Export to "geometry.xml"
openmc_geometry.export_to_xml()

# OpenMC simulation parameters
#batches = 60
#inactive = 10
#particles = 5000

# Very accurate parameters
batches = 550
inactive = 50
particles = 1200000

# Instantiate a Settings object
settings_file = openmc.Settings()
settings_file.batches = batches
settings_file.inactive = inactive
settings_file.particles = particles
settings_file.output = {'tallies': True}

# Extrapolate Temperature Data
settings_file.temperature = {'method': 'interpolation'}

# Create an initial uniform spatial source distribution over fissionable zones
bounds = [-100, -100, -100, 100, 100, 100] # default bounds were 0.63
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.Source(space=uniform_dist)

# Export to "settings.xml"
settings_file.export_to_xml()

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

# Instantiate more than a few different sections
total = mgxs.TotalXS(domain=cell, energy_groups=groups)
absorption = mgxs.AbsorptionXS(domain=cell, energy_groups=groups)
scattering = mgxs.ScatterXS(domain=cell, energy_groups=groups)
fission = mgxs.FissionXS(domain=cell, energy_groups=groups, nu=False)
capture = mgxs.CaptureXS(domain=cell, energy_groups=groups)
nu_fission_matrix = mgxs.NuFissionMatrixXS(domain=cell, energy_groups=groups)
scattering_matrix = mgxs.ScatterMatrixXS(domain=cell, energy_groups=groups)
inverse_velocity = mgxs.InverseVelocity(domain=cell, energy_groups=groups)
chi_prompt=mgxs.Chi(domain=cell, energy_groups=groups, prompt=True)
chi=mgxs.Chi(domain=cell, energy_groups=groups)
chi_delayed=mgxs.ChiDelayed(domain=cell, energy_groups=groups)
nuSigmaEff=mgxs.FissionXS(domain=cell, energy_groups=groups, nu=True)
diffusion_coefficient=mgxs.DiffusionCoefficient(domain=cell, energy_groups=groups) # should nu be true or false? - doesn't seem to make a difference
sigmaPow=mgxs.KappaFissionXS(domain=cell, energy_groups=groups)
beta=mgxs.Beta(domain=cell, energy_groups=one_group, delayed_groups=delayed_groups)
lambda1=mgxs.DecayRate(domain=cell, energy_groups=one_group, delayed_groups=delayed_groups)


# Instantiate an empty Tallies object
tallies_file = openmc.Tallies()

# Add total tallies to the tallies file
tallies_file += total.tallies.values()

# Add absorption tallies to the tallies file
tallies_file += absorption.tallies.values()

# Add scattering tallies to the tallies file
tallies_file += scattering.tallies.values()

# Add fission tallies to the tallies file
tallies_file += fission.tallies.values()

# Add capture tallies to the tallies file
tallies_file += capture.tallies.values()

# Add nu_fission_matrix tallies to the tallies file
tallies_file += nu_fission_matrix.tallies.values()

# Add scattering_matrix tallies to the tallies file
tallies_file += scattering_matrix.tallies.values()

# Add inverse velocity tallies to the tallies file
tallies_file += inverse_velocity.tallies.values()

# Add nuSigmaEff tallies to the tallies file
tallies_file += nuSigmaEff.tallies.values()

# Add chi prompt tallies to the tallies file
tallies_file += chi_prompt.tallies.values()

# Add chi tallies to the tallies file
tallies_file += chi.tallies.values()

# Add chi delayed tallies to the tallies file
tallies_file += chi_delayed.tallies.values()

# Add diffusion coefficient tallies to the tallies file
tallies_file += diffusion_coefficient.tallies.values()

# Add sigmaPow tallies to the tallies file
tallies_file += sigmaPow.tallies.values()

# Add beta tallies to the tallies file
tallies_file += beta.tallies.values()

# Add lambda tallies to the tallies file
tallies_file += lambda1.tallies.values()

# Export to "tallies.xml"
tallies_file.export_to_xml()

# Run OpenMC
openmc.run(mpi_args=['mpirun','-np','16'])

# Rename the statepoint file
os.rename('statepoint.550.h5', 'MAstatepoint.550.h5')

# Load the last statepoint file
sp = openmc.StatePoint('MAstatepoint.550.h5') #change this with the simulation parameters 'statepoint.n.h5' n=number of batches

# Load the tallies from the statepoint into each MGXS object
total.load_from_statepoint(sp)
absorption.load_from_statepoint(sp)
scattering.load_from_statepoint(sp)
fission.load_from_statepoint(sp)
capture.load_from_statepoint(sp)
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
total.export_xs_data(filename='total-xs', format='csv')
absorption.export_xs_data(filename='absorption-xs', format='csv')
scattering.export_xs_data(filename='scattering-xs', format='csv')
fission.export_xs_data(filename='fission-xs', format='csv')
capture.export_xs_data(filename='capture-xs', format='csv')
nu_fission_matrix.export_xs_data(filename='nufissionmatrix-xs', format='csv')
scattering_matrix.export_xs_data(filename='scatteringmatrix-xs', format='csv')
total.export_xs_data(filename='total-xs', format='csv')
inverse_velocity.export_xs_data(filename='inverse-velocity', format='csv')
chi_prompt.export_xs_data(filename='chi-prompt', format='csv')
chi.export_xs_data(filename='chi', format='csv')
chi_delayed.export_xs_data(filename='chi-delayed', format='csv')
nuSigmaEff.export_xs_data(filename='nu-SigmaEff', format='csv')
diffusion_coefficient.export_xs_data(filename='diffusion-coefficient', format='csv')
sigmaPow.export_xs_data(filename='sigmaPow', format='csv')
beta.export_xs_data(filename='betaone', format='csv')
lambda1.export_xs_data(filename='lambdaone', format='csv')

total.build_hdf5_store(filename='mgxs', append=True)
absorption.build_hdf5_store(filename='mgxs', append=True)
scattering.build_hdf5_store(filename='mgxs', append=True)
fission.build_hdf5_store(filename='mgxs', append=True)
capture.build_hdf5_store(filename='mgxs', append=True)
nu_fission_matrix.build_hdf5_store(filename='mgxs', append=True)
scattering_matrix.build_hdf5_store(filename='mgxs', append=True)
inverse_velocity.build_hdf5_store(filename='mgxs', append=True)
chi_prompt.build_hdf5_store(filename='mgxs', append=True)
chi.build_hdf5_store(filename='mgxs', append=True)
chi_delayed.build_hdf5_store(filename='mgxs', append=True)
nuSigmaEff.build_hdf5_store(filename='mgxs', append=True)
diffusion_coefficient.build_hdf5_store(filename='mgxs', append=True)
sigmaPow.build_hdf5_store(filename='mgxs', append=True)
beta.build_hdf5_store(filename='mgxs', append=True)
lambda1.build_hdf5_store(filename='mgxs', append=True)
