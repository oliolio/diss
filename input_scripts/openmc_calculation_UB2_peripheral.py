import os
import openmc
import numpy as np
import openmc.deplete

ROOT_DIR = "/home/hpcatta/openmc-0.13.3/please/last/"
DEPLETION_FILE = "/home/hpcatta/submit1/mchain_endfb80_pwr.xml"
CROSS_SECTION_FILE = '/home/hpcatta/openmc/endfb-viii.0-hdf5/cross_sections.xml'
ENRICHMENT_VALUES = "enrichments.txt"
RADIUS_VALUES = "fuel_radius.txt"


def get_cross_sections():
    os.environ['OPENMC_CROSS_SECTIONS'] = os.path.abspath(CROSS_SECTION_FILE)
    openmc.config['cross_sections'] = CROSS_SECTION_FILE

get_cross_sections()
# [8.598314932618374, 5.088872602860442, 5.834651344648201, 2.585123097889322, 0.1128986004645355]
# [12.987205938613467, 7.232301444035258, 3.3239793905082253, 5.2691647471320024, 0.1289319839497537, 0.9003274386310784

r1 = {{r1}}  #enrichment valuess
r2 = {{r2}}
r3 = {{r3}}
r4 = {{r4}}
t_urbo= {{t_urbo}}} #peripheral UB2 layer
atom_fraction_U= {{atom_fraction_U}}   #amount of B in UO2:UB2 layer


FUEL_AND_urbo_OUTER_RADIUS = 0.415
FR1 = 0.20475
FR2 = 0.28956
FR3 = 0.35464
FR4 = 0.4095

def save_radius_values(fr1, fr2, fr3, fr4, thickness_urbo):
    gen_pop = ['standalone', 'run']
    with open(os.path.join(ROOT_DIR, RADIUS_VALUES), 'a') as file:
        line = f"{gen_pop[0]} {gen_pop[1]} {fr1} {fr2} {fr3} {fr4} {FUEL_AND_urbo_OUTER_RADIUS} {thickness_urbo} {atom_fraction_U}\n"
        file.write(line)

def get_fuel_radius_values(urbo_thickness: float):
    fr4_new = FUEL_AND_urbo_OUTER_RADIUS - urbo_thickness
    alpha = fr4_new / FR4
    fr1 = FR1 * alpha
    fr2 = fr1 + (FR2 - FR1) * alpha
    fr3 = fr2 + (FR3 - FR2) * alpha
    fr4 = fr4_new
    return fr1, fr2, fr3, fr4

def run_calculation(r1__, r2__, r3__, r4__, thickness_urbo=None):
    if thickness_urbo is None:
        thickness_urbo = FUEL_AND_urbo_OUTER_RADIUS - FR4
    plot_pixels = 1200
    mesh_nr = 9

    def define_universe_components(r1_, r2_, r3_, r4_, no_of_pins, thickness_urbo_=FUEL_AND_urbo_OUTER_RADIUS - FR4, rho_UB2_diluted=12.71):
        fr1, fr2, fr3, fr4 = get_fuel_radius_values(thickness_urbo_)
        save_radius_values(fr1, fr2, fr3, fr4, thickness_urbo_)
        fuel_radius1 = openmc.ZCylinder(r=fr1)
        fuel_radius2 = openmc.ZCylinder(r=fr2)
        fuel_radius3 = openmc.ZCylinder(r=fr3)
        fuel_outer_radius = openmc.ZCylinder(r=fr4)
        urbo_radius = openmc.ZCylinder(r=FUEL_AND_urbo_OUTER_RADIUS)
        
        clad_inner_radius = openmc.ZCylinder(r=0.4177)
        clad_outer_radius = openmc.ZCylinder(r=0.4749)
        gt_inner_radius = openmc.ZCylinder(r=0.5715)
        gt_outer_radius = openmc.ZCylinder(r=0.612)
        lower_cutback_height = openmc.ZPlane(z0=-426.0 / 2.0 + 15)
        upper_cutback_height = openmc.ZPlane(z0=+426.0 / 2.0 - 15)
        fuel_bottom = openmc.ZPlane(z0=-426.0 / 2.0, boundary_type='reflective')
        fuel_top = openmc.ZPlane(z0=+426.0 / 2.0, boundary_type='reflective')
        
        materials = openmc.Materials([])
        
        f_1 = openmc.Material(name='Fuel 1')
        f_1.add_element('U', 1.0, enrichment=r1_)
        f_1.add_nuclide('O16', 2.0)
        f_1.set_density('g/cc', 10.4)
        f_1.volume = (np.pi * fuel_radius1.r ** 2) * 396 * no_of_pins
        materials += [f_1]
        
        f_2 = openmc.Material(name='Fuel 2')
        f_2.add_element('U', 1.0, enrichment=r2_)
        f_2.add_nuclide('O16', 2.0)
        f_2.set_density('g/cc', 10.4)
        f_2.volume = np.pi * (fuel_radius2.r ** 2 - fuel_radius1.r ** 2) * 396 * no_of_pins
        materials += [f_2]
        
        f_3 = openmc.Material(name='Fuel 3')
        f_3.add_element('U', 1.0, enrichment=r3_)
        f_3.add_nuclide('O16', 2.0)
        f_3.set_density('g/cc', 10.4)
        f_3.volume = np.pi * (fuel_radius3.r ** 2 - fuel_radius2.r ** 2) * 396 * no_of_pins
        materials += [f_3]
        
        f_4 = openmc.Material(name='Fuel 4')
        f_4.add_element('U', 1.0, enrichment=r4_)
        f_4.add_nuclide('O16', 2.0)
        f_4.set_density('g/cc', 10.4)
        f_4.volume = np.pi * (fuel_outer_radius.r ** 2 - fuel_radius3.r ** 2) * 396 * no_of_pins
        materials += [f_4]
        
        r_cutback = r1_ - 0.5
        
        f_cutback = openmc.Material(name='Fuel Cutback')
        f_cutback.add_element('U', 1.0, enrichment=r_cutback)
        f_cutback.add_nuclide('O16', 2.0)
        f_cutback.set_density('g/cc', 10.4)
        f_cutback.volume = (np.pi * fuel_outer_radius.r ** 2) * 30 * no_of_pins
        materials += [f_cutback]
        
        atom_fraction_UB2 = 1 - atom_fraction_U
        
        UB2_diluted = openmc.Material(name='UB2_diluted')
        UB2_diluted.add_nuclide('U235', (r4_/100) * atom_fraction_U + (r4_/100) * atom_fraction_UB2, 'ao')
        UB2_diluted.add_nuclide('U238', (1 - (r4_/100)) * atom_fraction_UB2 + (1 - (r4_/100)) * atom_fraction_U, 'ao')
        UB2_diluted.add_nuclide('B10', 0.16 * atom_fraction_UB2 * 2, 'ao')
        UB2_diluted.add_nuclide('B11', 0.84 * atom_fraction_UB2 * 2, 'ao')
        density_UB2 = 12.71
        density_U = 19.1
        density_UB2_diluted = (atom_fraction_UB2 * density_UB2) + (atom_fraction_U * density_U)
        UB2_diluted.set_density('g/cm3', density_UB2_diluted)
        UB2_diluted.depletable = True
        UB2_diluted.volume = np.pi * (urbo_radius.r ** 2 - fuel_outer_radius.r ** 2) * 426 * no_of_pins
        
        materials += [UB2_diluted]
        
        zirconium = openmc.Material(name="zirconium")
        zirconium.add_element('Zr', 1.0)
        zirconium.set_density('g/cm3', 6.6)
        
        helium = openmc.Material(name='Helium')
        helium.add_element('He', 1.0)
        helium.set_density('g/cm3', 0.178e-3)
        
        water = openmc.Material(name='water')
        water.add_nuclide('H1', 2.0)
        water.add_nuclide('O16', 1.0)
        water.set_density('g/cm3', 1.0)
        water.add_s_alpha_beta('c_H_in_H2O')
        
        materials += [helium, zirconium, water]
        
        fuel1_region = -fuel_radius1 & +lower_cutback_height & -upper_cutback_height
        fuel2_region = +fuel_radius1 & -fuel_radius2 & +lower_cutback_height & -upper_cutback_height
        fuel3_region = +fuel_radius2 & -fuel_radius3 & +lower_cutback_height & -upper_cutback_height
        fuel4_region = +fuel_radius3 & -fuel_outer_radius & +lower_cutback_height & -upper_cutback_height
        
        cutback_region = (-fuel_outer_radius & -lower_cutback_height & +fuel_bottom) | (-fuel_outer_radius & +upper_cutback_height & -fuel_top)
        
        urbo_region = +fuel_outer_radius & -urbo_radius & +fuel_bottom & -fuel_top
        gap_region = +urbo_radius & -clad_inner_radius & +fuel_bottom & -fuel_top
        clad_region = +clad_inner_radius & -clad_outer_radius & +fuel_bottom & -fuel_top
        
        gt_water_region = -gt_inner_radius & +fuel_bottom & -fuel_top
        gt_clad_region = +gt_inner_radius & -gt_outer_radius & +fuel_bottom & -fuel_top
        
        fuel1 = openmc.Cell(name='fuel1')
        fuel1.fill = f_1
        fuel1.region = fuel1_region
        
        fuel2 = openmc.Cell(name='fuel2')
        fuel2.fill = f_2
        fuel2.region = fuel2_region
        
        fuel3 = openmc.Cell(name='fuel3')
        fuel3.fill = f_3
        fuel3.region = fuel3_region
        
        fuel4 = openmc.Cell(name='fuel4')
        fuel4.fill = f_4
        fuel4.region = fuel4_region
        
        cutback = openmc.Cell(name='Cutback')
        cutback.fill = f_cutback
        cutback.region = cutback_region
        
        urbo = openmc.Cell(name='urbo coating')
        urbo.fill = UB2_diluted
        urbo.region = urbo_region
        
        gap = openmc.Cell(name='air gap')
        gap.region = gap_region
        gap.fill = helium
        
        clad = openmc.Cell(name='clad')
        clad.fill = zirconium
        clad.region = clad_region
        
        moderator = openmc.Cell(name='Moderator', fill=water,
                                region=+clad_outer_radius & -fuel_top & +fuel_bottom)
        
        return (gt_outer_radius, fuel_bottom, fuel_top, f_1, f_2, f_3, f_4, f_cutback, UB2_diluted,
                zirconium, water, gt_water_region, gt_clad_region, fuel1, fuel2, fuel3, fuel4, cutback,
                gap, clad, moderator, urbo)
    
    (gt_outer_radius, fuel_bottom, fuel_top, f_1, f_2, f_3, f_4, f_cutback, UB2_diluted,
     zirconium, water, gt_water_region, gt_clad_region, fuel1, fuel2, fuel3, fuel4, cutback,
     gap, clad, moderator, urbo) = define_universe_components(r1__, r2__, r3__, r4__, thickness_urbo_=1E-6, no_of_pins=53.5)
    
    (gt_outer_radius_, fuel_bottom_, fuel_top_, f_1_, f_2_, f_3_, f_4_, f_cutback_, UB2_diluted,
     zirconium_, water_, gt_water_region_, gt_clad_region_, fuel1_, fuel2_, fuel3_, fuel4_, cutback_,
     gap_, clad_, moderator_, urbo_) = define_universe_components(r1__, r2__, r3__, r4__, thickness_urbo_=thickness_urbo, no_of_pins=12.5)
    
    pin_cell_universe_outer = openmc.Universe(name='Fuel Pin (Outer)', cells=[
        fuel1, fuel2, fuel3, fuel4, cutback, urbo, gap, clad, moderator])
    
    pin_cell_universe_inner = openmc.Universe(name='Fuel Pin (Inner)', cells=[
        fuel1_, fuel2_, fuel3_, fuel4_, cutback_, urbo_, gap_, clad_, moderator_])
    
    gt1 = openmc.Cell(name='Guide tube filled with water', fill=water, region=gt_water_region)
    gt2 = openmc.Cell(name='Guide tube clad', fill=zirconium, region=gt_clad_region)
    mode = openmc.Cell(name='Moderator', fill=water,
                       region=+gt_outer_radius & -fuel_top & +fuel_bottom)
    gt = openmc.Universe(name='Guide Tube', cells=[gt1, gt2, mode])
    
    pitch = 1.26
    
    assembly = openmc.RectLattice(name='Fuel - 0BA')
    assembly.pitch = (pitch, pitch)
    assembly.lower_left = [-pitch * mesh_nr / 2] * 2
    
    assembly.universes = [
        [pin_cell_universe_outer for _ in range(mesh_nr)],
        [pin_cell_universe_outer for _ in range(mesh_nr)],
        [pin_cell_universe_outer for _ in range(5)] + [gt] + [pin_cell_universe_inner for _ in range(2)] + [gt],
        [pin_cell_universe_outer for _ in range(3)] + [gt] + [pin_cell_universe_inner for _ in range(5)],
        [pin_cell_universe_outer for _ in range(mesh_nr)],
        [pin_cell_universe_outer for _ in range(2)] + [gt] + [pin_cell_universe_inner for _ in range(2)] + [gt] + [pin_cell_universe_inner for _ in range(3)],
        [pin_cell_universe_outer for _ in range(mesh_nr)],
        [pin_cell_universe_outer for _ in range(mesh_nr)],
        [pin_cell_universe_outer for _ in range(2)] + [gt] + [pin_cell_universe_inner for _ in range(2)] + [gt] + [pin_cell_universe_outer for _ in range(3)],
    ]
    
    root_cell = openmc.Cell(name='root cell', fill=assembly)
    min_x = openmc.XPlane(x0=-pitch * mesh_nr / 2, boundary_type='reflective')
    max_x = openmc.XPlane(x0=pitch * mesh_nr / 2, boundary_type='reflective')
    min_y = openmc.YPlane(y0=-pitch * mesh_nr / 2, boundary_type='reflective')
    max_y = openmc.YPlane(y0=pitch * mesh_nr / 2, boundary_type='reflective')
    
    root_cell.region = +min_x & -max_x & +min_y & -max_y & +fuel_bottom & -fuel_top
    
    root_universe = openmc.Universe(name='root universe')
    root_universe.add_cell(root_cell)
    
    geometry = openmc.Geometry(root_universe)
    geometry.export_to_xml()
    
    point = openmc.stats.Point((0, 0, 0))
    source = openmc.source.Source(space=point)
    
    settings = openmc.Settings()
    settings.source = source
    settings.batches = 300
    settings.inactive = 100
    settings.particles = 10000
    settings.export_to_xml()
    
    tallies = openmc.Tallies()
    mesh = openmc.RegularMesh(mesh_id=1)
    mesh.dimension = [mesh_nr, mesh_nr]
    mesh.lower_left = [-pitch * mesh_nr / 2] * 2
    mesh.width = [pitch, pitch]
    
    mesh_filter = openmc.MeshFilter(mesh)
    
    tally = openmc.Tally(name='mesh tally')
    tally.filters = [mesh_filter]
    tally.scores = ['fission']
    
    tallies.append(tally)
    tallies.export_to_xml()
    
    plot1 = openmc.Plot()
    plot1.filename = 'pinplot'
    plot1.origin = (0, 0, 0)
    plot1.width = (pitch * mesh_nr, pitch * mesh_nr)
    plot1.pixels = (plot_pixels, plot_pixels)
    plot1.color_by = 'material'
    plot1.colors = {f_1: 'red', f_2: 'orange', f_3: 'yellow', f_4: 'green', UB2_diluted: 'maroon', water: 'blue',
                    f_cutback: 'gray'}
    
    plot2 = openmc.Plot()
    plot2.filename = 'XY basis'
    plot2.basis = 'xy'
    plot2.origin = (0, 0, 0)
    plot2.width = (pitch * mesh_nr, pitch * mesh_nr)
    plot2.pixels = (plot_pixels, plot_pixels)
    plot2.color_by = 'material'
    plot2.colors = {f_1: 'red', f_2: 'orange', f_3: 'yellow', f_4: 'green', UB2_diluted: 'maroon', water: 'blue',
                    f_cutback: 'gray'}
    
    plot3 = openmc.Plot()
    plot3.filename = 'YZ basis'
    plot3.basis = 'yz'
    plot3.width = (pitch * mesh_nr, 426)
    plot3.pixels = (plot_pixels, plot_pixels)
    plot3.color_by = 'cell'
    plot3.colors = {f_1: 'red', f_2: 'orange', f_3: 'yellow', f_4: 'green', UB2_diluted: 'maroon', water: 'blue',
                    f_cutback: 'gray'}
    
    plots = openmc.Plots([plot1, plot2, plot3])
    plots.export_to_xml()
    
    model = openmc.Model(geometry=geometry, settings=settings)
    model.export_to_xml()
    operator = openmc.deplete.CoupledOperator(model, DEPLETION_FILE)
    power = 21.75E06 / 4
    time_steps = [0.5, 1.5, 1, 27, 60, 90, 90, 90, 90, 30, 67.5]    
    integrator = openmc.deplete.CECMIntegrator(operator, time_steps, power, timestep_units='d')
    integrator.integrate()

def save_enrichments():
    """enrichments saved into file."""
    try:
        gen_pop = os.path.split(os.path.abspath(os.curdir))[-1].split("_")
        gen = gen_pop[0]
        pop = gen_pop[1]
    except IndexError:
        gen = "unknown"
        pop = "unknown"
    
    with open(os.path.join(ROOT_DIR, ENRICHMENT_VALUES), 'a') as file:
        line = f"{gen} {pop} {r1} {r2} {r3} {r4} {t_urbo}\n"
        file.write(line)

save_enrichments()
run_calculation(r1, r2, r3, r4, thickness_urbo=t_urbo)