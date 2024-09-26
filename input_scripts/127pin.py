import openmc
import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Create materials for the problem

# Enrichment values from the assembly file (in weight percent)
enr1 = 12.77416936655832  # Enrichment for Fuel 1
enr2 = 6.971434136971359   # Enrichment for Fuel 2
enr3 = 2.8587658927920434  # Enrichment for Fuel 3
enr4 = 7.390724339693737   # Enrichment for Fuel 4

# Atom fraction of uranium in the UB2_diluted material from the assembly file
atom_fraction_U = 0.9006208156483683

# Fuel density
fuel_density = 10.4  # g/cm³

# Create fuel materials with updated enrichments and density
f1 = openmc.Material(name='Fuel 1')
f1.set_density('g/cm3', fuel_density)
f1.add_element('U', 1.0, enrichment=enr1)
f1.add_element('O', 2.0)

f2 = openmc.Material(name='Fuel 2')
f2.set_density('g/cm3', fuel_density)
f2.add_element('U', 1.0, enrichment=enr2)
f2.add_element('O', 2.0)

f3 = openmc.Material(name='Fuel 3')
f3.set_density('g/cm3', fuel_density)
f3.add_element('U', 1.0, enrichment=enr3)
f3.add_element('O', 2.0)

f4 = openmc.Material(name='Fuel 4')
f4.set_density('g/cm3', fuel_density)
f4.add_element('U', 1.0, enrichment=enr4)
f4.add_element('O', 2.0)

# Create the UB2_diluted material for the coating
atom_fraction_UB2 = 1.0 - atom_fraction_U

# Calculate the fractions of U-235 and U-238 in the UB2_diluted material
# Adjusted fractions for clarity and correctness
U235_fraction = ((enr4 / 100) * atom_fraction_U) + ((enr4 / 100) * atom_fraction_UB2)
U238_fraction = ((1.0 - (enr4 / 100)) * atom_fraction_U) + ((1.0 - (enr4 / 100)) * atom_fraction_UB2)

# Boron isotopes fractions assuming natural boron in UB₂
B10_fraction = 0.199 * atom_fraction_UB2 * 2.0  # Multiply by 2 for B₂
B11_fraction = 0.801 * atom_fraction_UB2 * 2.0

# Normalize fractions so they sum to 1
total_fraction = U235_fraction + U238_fraction + B10_fraction + B11_fraction

UB2_diluted = openmc.Material(name='UB2_diluted')
UB2_diluted.add_nuclide('U235', U235_fraction / total_fraction, 'ao')
UB2_diluted.add_nuclide('U238', U238_fraction / total_fraction, 'ao')
UB2_diluted.add_nuclide('B10', B10_fraction / total_fraction, 'ao')
UB2_diluted.add_nuclide('B11', B11_fraction / total_fraction, 'ao')

# Calculate the density of the UB2_diluted material
density_UB2 = 12.71  # g/cm³ for UB₂
density_U = 19.1     # g/cm³ for uranium
density_UB2_diluted = (atom_fraction_UB2 * density_UB2) + (atom_fraction_U * density_U)
UB2_diluted.set_density('g/cm3', density_UB2_diluted)

# Helium gas for the gap
helium = openmc.Material(name='Helium for gap')
helium.set_density('g/cm3', 0.000178)  # g/cm³
helium.add_element('He', 1.0)

# Zirconium cladding
zirconium = openmc.Material(name='zirconium')
zirconium.set_density('g/cm3', 6.6)
zirconium.add_element('Zr', 1.0)

# Water moderator (non-borated)
water = openmc.Material(name='water')
water.set_density('g/cm3', 1.0)
water.add_element('H', 2.0)
water.add_element('O', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')

# Collect materials together and export to XML
materials = openmc.Materials([f1, f2, f3, f4, UB2_diluted, helium, zirconium, water])
materials.export_to_xml()

###############################################################################
# Define problem geometry

# Fuel radii calculations from the assembly file parameters
t_urbo = 0.12530432387687737 * 1.0  # Thickness of the UB2_diluted layer (in cm)

# Fuel radii constants from the assembly file (in cm)
FUEL_AND_urbo_OUTER_RADIUS = 0.415 * 100  # Convert from meters to cm
FR1 = 0.20475 * 100
FR2 = 0.28956 * 100
FR3 = 0.35464 * 100
FR4 = 0.4095 * 100

# Function to calculate fuel radii
def get_fuel_radius_values(urbo_thickness):
    fr4_new = FUEL_AND_urbo_OUTER_RADIUS - urbo_thickness
    alpha = fr4_new / FR4
    fr1 = FR1 * alpha
    fr2 = fr1 + (FR2 - FR1) * alpha
    fr3 = fr2 + (FR3 - FR2) * alpha
    fr4 = fr4_new
    return fr1, fr2, fr3, fr4

# Calculate the radii (in cm)
fr1, fr2, fr3, fr4 = get_fuel_radius_values(t_urbo)
fuel_or_radius = fr4  # Outer radius of the fuel (in cm)

# Cladding dimensions (in cm)
clad_ir_radius = 0.4177 * 100  # Inner radius of cladding
clad_or_radius = 0.47490 * 100  # Outer radius of cladding

# Create cylindrical surfaces for fuel regions
r1 = openmc.ZCylinder(r=fr1)
r2 = openmc.ZCylinder(r=fr2)
r3 = openmc.ZCylinder(r=fr3)
r4 = openmc.ZCylinder(r=fr4)
fuel_or = openmc.ZCylinder(r=fuel_or_radius, name='Fuel OR')

# Define outer radius (includes UB2_diluted layer)
urbo_radius = FUEL_AND_urbo_OUTER_RADIUS
fuel_urbo_or = openmc.ZCylinder(r=urbo_radius, name='Fuel + UB2_diluted OR')

# Create surfaces for cladding
clad_ir = openmc.ZCylinder(r=clad_ir_radius, name='Clad IR')
clad_or = openmc.ZCylinder(r=clad_or_radius, name='Clad OR')

# Create a square prism boundary for the pin cell
pitch = 1.21 * 100  # Convert pitch to cm
box = openmc.model.rectangular_prism(pitch, pitch, boundary_type='reflective')

# Define regions and cells
fuel1_region = -r1
fuel2_region = +r1 & -r2
fuel3_region = +r2 & -r3
fuel4_region = +r3 & -r4
urbo_region = +r4 & -fuel_or
gap_region = +fuel_or & -clad_ir
clad_region = +clad_ir & -clad_or
water_region = +clad_or & box

# Create cells
fuel1_cell = openmc.Cell(name='Fuel1', fill=f1, region=fuel1_region)
fuel2_cell = openmc.Cell(name='Fuel2', fill=f2, region=fuel2_region)
fuel3_cell = openmc.Cell(name='Fuel3', fill=f3, region=fuel3_region)
fuel4_cell = openmc.Cell(name='Fuel4', fill=f4, region=fuel4_region)
urbo_cell = openmc.Cell(name='UB2_diluted coating', fill=UB2_diluted, region=urbo_region)
gap_cell = openmc.Cell(name='Helium gap', fill=helium, region=gap_region)
clad_cell = openmc.Cell(name='Cladding', fill=zirconium, region=clad_region)
water_cell = openmc.Cell(name='Moderator', fill=water, region=water_region)

# Create root universe
root_universe = openmc.Universe(cells=[fuel1_cell, fuel2_cell, fuel3_cell, fuel4_cell,
                                       urbo_cell, gap_cell, clad_cell, water_cell])

# Define geometry and export to XML
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

###############################################################################
# Define simulation settings

settings = openmc.Settings()
settings.batches = 150
settings.inactive = 50
settings.particles = 5000

# Define a point neutron source at the center of the fuel pin
point = openmc.stats.Point((0.0, 0.0, 0.0))
source = openmc.Source(space=point)
settings.source = source

# Enable Shannon entropy calculation for source convergence checks
entropy_mesh = openmc.RegularMesh()
entropy_mesh.dimension = [10, 10, 1]
entropy_mesh.lower_left = (-pitch / 2.0, -pitch / 2.0, -1.0)
entropy_mesh.upper_right = (pitch / 2.0, pitch / 2.0, 1.0)
settings.entropy_mesh = entropy_mesh

# Export settings to XML
settings.export_to_xml()

###############################################################################
# Define tallies

# Define a mesh for tallying fission reactions across the fuel diameter
mesh = openmc.RegularMesh()
mesh.dimension = [100, 1, 1]
mesh.lower_left = [-fuel_or_radius, -0.05, -0.05]
mesh.upper_right = [fuel_or_radius, 0.05, 0.05]

# Create a mesh filter for the tally
mesh_filter = openmc.MeshFilter(mesh)

# Create a tally to record fission reaction rates
mesh_tally = openmc.Tally(name='Fission Rate Mesh Tally')
mesh_tally.filters = [mesh_filter]
mesh_tally.scores = ['fission']

# Export tallies to XML
tallies = openmc.Tallies([mesh_tally])
tallies.export_to_xml()

###############################################################################
# Define plots (optional for visualization)

# Plot configuration to visualize the materials and geometry
plot = openmc.Plot()
plot.filename = 'pin_cell_plot'
plot.origin = (0, 0, 0)
plot.width = (pitch, pitch)
plot.pixels = (400, 400)
plot.color_by = 'material'
plot.colors = {
    f1: 'red',
    f2: 'orange',
    f3: 'yellow',
    f4: 'green',
    UB2_diluted: 'purple',
    helium: 'blue',
    zirconium: 'gray',
    water: 'lightblue'
}

# Export plots to XML
plots = openmc.Plots([plot])
plots.export_to_xml()

# Generate the plot (uncomment the following line if you want to generate the plot before running)
# openmc.plot_geometry()

###############################################################################
# Run the OpenMC simulation

# Run OpenMC using the generated XML files
openmc.run()

###############################################################################
# Process tally results and plot the radial distribution of fission rate

# Open the statepoint file (update the filename if necessary)
sp = openmc.StatePoint('statepoint.150.h5')

# Get the tally by name
tally = sp.get_tally(name='Fission Rate Mesh Tally')

# Extract the tally data
fission_rate = tally.mean.ravel()
fission_rate_sd = tally.std_dev.ravel()

# Get the mesh grid positions
mesh_filter = tally.find_filter(openmc.MeshFilter)
mesh = mesh_filter.mesh
# x_grid contains the boundaries of the mesh bins in x-direction
x_grid = np.linspace(mesh.lower_left[0], mesh.upper_right[0], num=mesh.dimension[0]+1)
# Calculate the center positions of each mesh bin
x_centers = 0.5 * (x_grid[:-1] + x_grid[1:])

# Plot the radial distribution of fission rate
plt.figure(figsize=(8,6))
plt.errorbar(x_centers, fission_rate, yerr=fission_rate_sd, fmt='o-', ecolor='red', capsize=2)
plt.xlabel('Radial Position (cm)')
plt.ylabel('Fission Rate (per unit volume)')
plt.title('Radial Distribution of Fission Rate')
plt.grid(True)
plt.tight_layout()

# Save the plot to a file (you can change the filename as needed)
plt.savefig('fission_rate_radial_distribution.png')
