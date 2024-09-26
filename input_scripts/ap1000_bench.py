import openmc
import openmc.deplete

#2.35% , IFBA version (no pyrex)

# Define materials
f_235 = openmc.Material(name='Fuel 2.35%')
f_235.add_element('U', 1.0, enrichment=2.35)
f_235.add_element('O', 2.0)
f_235.set_density('g/cc', 10.4)

f_158 = openmc.Material(name='Fuel 1.58%')
f_158.add_element('U', 1.0, enrichment=1.58)
f_158.add_element('O', 2.0)
f_158.set_density('g/cc', 10.4)

helium = openmc.Material(name='Helium')
helium.add_element('He', 1.0)
helium.set_density('g/cm3', 0.178e-3)

# Update the cladding material to Zircaloy with MCNP composition
zircaloy = openmc.Material(name='Zircaloy-4')
zircaloy.add_element('Zr', 0.9811, percent_type='wo')
zircaloy.add_element('Sn', 0.0067, percent_type='wo')
zircaloy.add_element('Nb', 0.01, percent_type='wo')
zircaloy.add_element('O', 0.0012, percent_type='wo')
zircaloy.set_density('g/cm3', 6.55)

water = openmc.Material(name='Water')
water.add_nuclide('H1', 2.0)
water.add_nuclide('O16', 1.0)
water.set_density('g/cm3', 1.0)
water.add_s_alpha_beta('c_H_in_H2O')

materials = openmc.Materials([f_235, f_158, helium, zircaloy, water])
materials.export_to_xml()

# Define geometry
fuel_outer_radius = openmc.ZCylinder(r=0.4095)
ifba_radius = openmc.ZCylinder(r=0.4010)
clad_inner_radius = openmc.ZCylinder(r=0.4177)
clad_outer_radius = openmc.ZCylinder(r=0.4749)

lower_cutback_height = openmc.ZPlane(z0=-426.0/2.0+15)
upper_cutback_height = openmc.ZPlane(z0=+426.0/2.0-15)
fuel_bottom = openmc.ZPlane(z0=-426.0/2.0, boundary_type='reflective')
fuel_top = openmc.ZPlane(z0=+426.0/2.0, boundary_type='reflective')

# Fuel regions
fuel_region = -fuel_outer_radius & +lower_cutback_height & -upper_cutback_height
lower_cutback_region = -fuel_outer_radius & -lower_cutback_height & +fuel_bottom
upper_cutback_region = -fuel_outer_radius & +upper_cutback_height & -fuel_top

# Define cells
fuel = openmc.Cell(name='Fuel')
fuel.fill = f_235
fuel.region = fuel_region

lower_cutback = openmc.Cell(name='Lower Cutback')
lower_cutback.fill = f_158
lower_cutback.region = lower_cutback_region

upper_cutback = openmc.Cell(name='Upper Cutback')
upper_cutback.fill = f_158
upper_cutback.region = upper_cutback_region

gap_region = +fuel_outer_radius & -clad_inner_radius & +fuel_bottom & -fuel_top
gap = openmc.Cell(name='Gap')
gap.region = gap_region
gap.fill = helium

clad_region = +clad_inner_radius & -clad_outer_radius & +fuel_bottom & -fuel_top
clad = openmc.Cell(name='Cladding')
clad.fill = zircaloy
clad.region = clad_region

moderator_region = +clad_outer_radius & +fuel_bottom & -fuel_top
moderator = openmc.Cell(name='Moderator')
moderator.fill = water
moderator.region = moderator_region

# Create pin cell universe
pin_cell_universe = openmc.Universe(cells=[fuel, lower_cutback, upper_cutback, gap, clad, moderator])

# Define guide tube geometry
gt_inner_radius = openmc.ZCylinder(r=0.5715)
gt_outer_radius = openmc.ZCylinder(r=0.612)

gt_water_region = -gt_inner_radius & +fuel_bottom & -fuel_top
gt_water = openmc.Cell(name='Guide Tube Water')
gt_water.fill = water
gt_water.region = gt_water_region

gt_clad_region = +gt_inner_radius & -gt_outer_radius & +fuel_bottom & -fuel_top
gt_clad = openmc.Cell(name='Guide Tube Cladding')
gt_clad.fill = zircaloy
gt_clad.region = gt_clad_region

gt_moderator_region = +gt_outer_radius & +fuel_bottom & -fuel_top
gt_moderator = openmc.Cell(name='Guide Tube Moderator')
gt_moderator.fill = water
gt_moderator.region = gt_moderator_region

gt_universe = openmc.Universe(cells=[gt_water, gt_clad, gt_moderator])

# Lattice parameters
pitch = 1.26
assembly = openmc.RectLattice(name='Fuel Assembly')
assembly.pitch = (pitch, pitch)
assembly.lower_left = [-pitch * 17 / 2] * 2

# Define lattice pattern
pin = pin_cell_universe
gt = gt_universe

assembly.universes = [
    [pin]*17,
    [pin]*17,
    [pin]*5 + [gt] + [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*5,
    [pin]*3 + [gt] + [pin]*9 + [gt] + [pin]*3,
    [pin]*17,
    [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*2,
    [pin]*17,
    [pin]*17,
    [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*2,
    [pin]*17,
    [pin]*17,
    [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*2,
    [pin]*17,
    [pin]*3 + [gt] + [pin]*9 + [gt] + [pin]*3,
    [pin]*5 + [gt] + [pin]*2 + [gt] + [pin]*2 + [gt] + [pin]*5,
    [pin]*17,
    [pin]*17,
]

# Root cell and universe
min_x = openmc.XPlane(x0=-10.71, boundary_type='reflective')
max_x = openmc.XPlane(x0=+10.71, boundary_type='reflective')
min_y = openmc.YPlane(y0=-10.71, boundary_type='reflective')
max_y = openmc.YPlane(y0=+10.71, boundary_type='reflective')

assembly_cell = openmc.Cell()
assembly_cell.fill = assembly
assembly_cell.region = +min_x & -max_x & +min_y & -max_y & +fuel_bottom & -fuel_top

root_universe = openmc.Universe(cells=[assembly_cell])

# Create geometry and export
geometry = openmc.Geometry(root_universe)
geometry.export_to_xml()

# Define settings
point = openmc.stats.Point((0, 0, 0))
source = openmc.Source(space=point)

settings = openmc.Settings()
settings.source = source
settings.batches = 150
settings.inactive = 50
settings.particles = 50000
settings.export_to_xml()

# Plotting (optional)
plot = openmc.Plot()
plot.filename = 'XY_basis'
plot.basis = 'xy'
plot.width = (pitch * 17, pitch * 17)
plot.pixels = (1200, 1200)
plot.color_by = 'material'

plot1 = openmc.Plot()
plot1.filename = 'XZ_basis'
plot1.basis = 'xz'
plot1.width = (pitch * 17, 426)
plot1.pixels = (1200, 1200)
plot1.color_by = 'material'
plot1.colors = {f_158: 'yellow', f_235: 'red'}

plots = openmc.Plots([plot, plot1])
plots.export_to_xml()


#openmc.run()
