import openmc


#### Geometry Section ####

## Surfaces ##

# Heat Exchanger
HX_top=openmc.ZPlane(z0=0.85*100, boundary_type='reflective', albedo=0.000001)
HX_bottom=openmc.ZPlane(z0=-0.85*100, boundary_type='reflective', albedo=0.000001)
HX_outside=openmc.ZCylinder(r=2.27*100, boundary_type='reflective', albedo=0.000001)
HX_inside=openmc.ZCylinder(r=1.835*100, boundary_type='reflective', albedo=0.000001)

# Box above and below HX
box_top=openmc.ZPlane(z0=0.95*100, boundary_type='reflective', albedo=0.000001)
box_bottom=openmc.ZPlane(z0=-0.95*100, boundary_type='reflective', albedo=0.000001)
box_outside=openmc.ZCylinder(r=2.17*100, boundary_type='reflective', albedo=0.000001)
box_inside=openmc.ZCylinder(r=1.935*100, boundary_type='reflective', albedo=0.000001)

# Curve above the HX
inner_curve_top=openmc.ZTorus(a=1.785*100,x0=0,y0=0,z0=0.95*100,b=0.15*100,c=0.15*100, boundary_type='reflective', albedo=0.000001)
outer_curve_top=openmc.ZTorus(a=1.785*100,x0=0,y0=0,z0=0.95*100,b=0.385*100,c=0.385*100, boundary_type='reflective', albedo=0.000001)
# use box_top as bottom of curve
curve_to_slant_bound=openmc.ZCylinder(r=1.785*100)

# Curve below the HX
inner_curve_bottom=openmc.ZTorus(a=1.785*100,x0=0,y0=0,z0=-0.95*100,b=0.15*100,c=0.15*100, boundary_type='reflective', albedo=0.000001)
outer_curve_bottom=openmc.ZTorus(a=1.785*100,x0=0,y0=0,z0=-0.95*100,b=0.385*100,c=0.385*100, boundary_type='reflective', albedo=0.000001)

# Core planes
top_wall=openmc.ZCone(x0=0,y0=0,z0=0.8*100,r2=11.56, boundary_type='reflective', albedo=0.666666)
top_cone_midplane=openmc.ZPlane(z0=0.8*100)
bottom_wall=openmc.ZCone(x0=0,y0=0,z0=-0.8*100,r2=11.56, boundary_type='reflective', albedo=0.666666)
bottom_cone_midplane=openmc.ZPlane(z0=-0.8*100)
core_outerCyl=openmc.ZTorus(a=2.32*100,x0=0,y0=0,z0=0,b=1.265*100,c=1.265*100, boundary_type='reflective', albedo=0.666666)
# Small planes parallel to walls
top_small=openmc.ZCone(x0=0,y0=0,z0=0.575*100,r2=11.56, boundary_type='reflective', albedo=0.666666)
top_small_cone_midplane=openmc.ZPlane(z0=0.575*100)
bottom_small=openmc.ZCone(x0=0,y0=0,z0=-0.575*100,r2=11.56, boundary_type='reflective', albedo=0.666666)
bottom_small_cone_midplane=openmc.ZPlane(z0=-0.575*100)
slant_to_core_plane=openmc.ZCylinder(r=1.585*100)
midplane=openmc.ZPlane(z0=0)

## Regions/Volumes ##
HX_region = -HX_top & +HX_inside & -HX_outside & +HX_bottom 
box_region = -box_top & +box_inside & -box_outside & +box_bottom
curve_top_region = +inner_curve_top & -outer_curve_top & +curve_to_slant_bound & +box_top
curve_bottom_region = +inner_curve_bottom & -outer_curve_bottom & +curve_to_slant_bound & -box_bottom
slant_top_region = +slant_to_core_plane & -curve_to_slant_bound & (+top_wall | -top_cone_midplane) & (-top_small & +top_small_cone_midplane) & +midplane
slant_bottom_region = +slant_to_core_plane & -curve_to_slant_bound & (+bottom_wall | +bottom_cone_midplane) & (-bottom_small & -bottom_small_cone_midplane) & -midplane
core_region=(+top_wall | -top_cone_midplane) & (+bottom_wall | +bottom_cone_midplane) & -slant_to_core_plane & +core_outerCyl
# full_region = HX_region & box_region & core_region & curve_top_region & slant_top_region & curve_bottom_region & slant_bottom_region
full_region = HX_region | box_region | core_region | curve_top_region | slant_top_region | curve_bottom_region | slant_bottom_region

## Material ## 

# Create and mix materials
# fuel
fuel = openmc.Material(1,"fuel", temperature=873.15)
fuel.set_density('g/cm3', 1.62134)
fuel.add_nuclide('Np237', 0.0642, 'wo')
fuel.add_nuclide('Pu238', 0.0318, 'wo')
fuel.add_nuclide('Pu239', 0.4393, 'wo')
fuel.add_nuclide('Pu240', 0.2127, 'wo')
fuel.add_nuclide('Pu241', 0.1352, 'wo')
fuel.add_nuclide('Pu242', 0.0788, 'wo')
fuel.add_nuclide('Am241', 0.0055, 'wo')
fuel.add_nuclide('Am243', 0.0233, 'wo')
fuel.add_nuclide('Cm244', 0.0092, 'wo')

li_salt = openmc.Material(2, "li_salt", temperature=873.15) # these temps don't matter
li_salt.add_nuclide('Li7', 0.9241)
li_salt.add_nuclide('Li6', 0.0759)

be_salt = openmc.Material(14, "be_salt", temperature=873.15)
be_salt.add_element('Be', 1)

f_salt = openmc.Material(114, "f_salt", temperature=873.15)
f_salt.add_element('F', 1)
# Ratio to get keff=1
Ratio=0.2
# Ratio=0.07
# MOSART paper lists ~73 LiF 27 BeF2
MA = 1
Li = (MA/Ratio) * .73
Be = (MA/Ratio) * .27
F = MA*3 + Li + Be # This comes from chemical formula
tot = Li + Be + MA + F
mixList = [MA/tot, Li/tot, Be/tot, F/tot]
mixList[3] = mixList[3] + (1-sum(mixList)) 
fuel_salt = openmc.Material.mix_materials([fuel, li_salt, be_salt, f_salt], mixList, 'ao') # default is 'ao'
# Define fuel salt temperature for density calculation
FStemp = 900 # Kelvin
# Define fuel salt density from (11) in https://www.oecd-nea.org/pt/docs/iem/jeju02/session3/SessionIII-12.pdf
rho_salt = 3.970-(7.20*10**-4)*FStemp
print(rho_salt)
fuel_salt.set_density('g/cm3', rho_salt)
# Set overall fuel salt temperature
fuel_salt.temperature = FStemp # Original Value = 900 K

# Add materials 
mats = openmc.Materials()
mats += [fuel, li_salt, be_salt, f_salt, fuel_salt]
mats.export_to_xml()

## Cells ## 

full_MSFR_3D = openmc.Cell(name='3D_MSFR')
full_MSFR_3D.region = full_region
# full_MSFR_3D.region = slant_top_region | slant_bottom_region
# full_MSFR_3D.region = core_region
full_MSFR_3D.fill = fuel_salt

## Geometry Model Setup ##

root_universe = openmc.Universe(cells=[full_MSFR_3D])
geometry = openmc.Geometry()
geometry.root_universe = root_universe
geometry.export_to_xml()

## Plots ##
plot = openmc.Plot()
plot.basis = 'yz'
plot.origin = (0,0,0)
plot.width = (5*100,5*100)
plot.pixels = (400,400)
plots = openmc.Plots([plot])
plots.export_to_xml()

## Run Settings ##

# Set the source/k-code settings
settings=openmc.Settings()
# Set the source type and location
point = openmc.stats.Point((0,0,0))
source = openmc.Source(space=point)
settings.source = source
# settings.batches = 250
settings.batches = 400
settings.inactive = 100
settings.particles = 10000
# settings.particles = 5000
## Use previous yaml (from Josh)
# settings.run_mode = 'fixed source'
settings.temperature = {'method' : 'interpolation'}
settings.export_to_xml()

tally2 = openmc.Tally(2)
# tally2.filters = [cell_filter]
tally2.scores = ['heating', 'fission']
# Export tallies
# tallies = openmc.Tallies([tally1, tally2])
tallies = openmc.Tallies([tally2])
tallies.export_to_xml()

openmc.run()