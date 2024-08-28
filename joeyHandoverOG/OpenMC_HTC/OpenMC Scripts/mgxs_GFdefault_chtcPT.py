# %matplotlib inline
from fileinput import filename
from operator import inv
import numpy as np
import matplotlib.pyplot as plt

import openmc
import sys
import openmc.mgxs as mgxs

# Create and mix materials
li_salt = openmc.Material(2, "li_salt", temperature=873.15) # these temps don't matter
li_salt.add_nuclide('Li7', 0.9241)
li_salt.add_nuclide('Li6', 0.0759)

f_salt = openmc.Material(114, "f_salt", temperature=873.15)
f_salt.add_element('F', 1)

th_salt = openmc.Material(5, "th_salt", temperature=873.15)
th_salt.add_nuclide('Th232', 1)

u_salt = openmc.Material(6, "u_salt", temperature=873.15)
u_salt.add_nuclide('U233', 1)


# Nominal Composition is 77.5% LiF, 20% ThF4, 2.5% UF4
#fuel_salt = openmc.Material.mix_materials([li_salt, th_salt, u_salt, f_salt], [0.060137, 0.518889, 0.065138, 0.355836], 'wo') # default is 'ao'

# Try increasing amount of U233 to get closer to keff=1
fuel_salt = openmc.Material.mix_materials([li_salt, th_salt, u_salt, f_salt], [0.060137, 0.46, 0.124027, 0.355836], 'wo')

# Define fuel salt temperature for density calculation
FStemp = 900 # Kelvin
FStempC = FStemp - 273.15

# Density Calculation from https://www.sciencedirect.com/science/article/pii/S0306454913004118
#rho_salt = (4094 - (8.82*10**-1)*(FStemp-1008))*10**-3

#fuel_salt.set_density('g/cm3', rho_salt)

# for constant density at perturbed temp (900K density)
fuel_salt.set_density('g/cm3', 4.189256)

# Set overall fuel salt temperature
#fuel_salt.temperature = FStemp # Original Value = 900 K

# for constant temp at perturbed density
fuel_salt.temperature = 1100


# Add materials 
mats = openmc.Materials()
mats += [li_salt, th_salt, u_salt, f_salt, fuel_salt]
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
bounds = [-100, -100, -100, 100, 100, 100] 
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings_file.source = openmc.Source(space=uniform_dist)

# Export to "settings.xml"
settings_file.export_to_xml()

# Instantiate a 6-group EnergyGroups object
groups = mgxs.EnergyGroups()
groups.group_edges = np.array([0., 748.5, 5531.0, 24790.0, 497900.0, 2.231e6, 20.0e6]) #need to change these for thermal/epithermal reactor?

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

# Load the last statepoint file
sp = openmc.StatePoint('statepoint.550.h5') #change this with the simulation parameters 'statepoint.n.h5' n=number of batches

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
