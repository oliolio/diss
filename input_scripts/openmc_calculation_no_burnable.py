import os
import openmc
import numpy as np
import openmc.deplete

# User-defined paths and filenames
ROOT_DIR = "/home/hpcatta/openmc-0.13.3/please/pretty"
DEPLETION_FILE = "/home/hpcatta/openmc-0.13.3/please/chain_endfb80_pwr.xml"
CROSS_SECTION_FILE = '/home/hpcatta/openmc/endfb-viii.0-hdf5/cross_sections.xml'

ENRICHMENT_VALUES = "enrichments.txt"
RADIUS_VALUES = "fuel_radius.txt"

def get_cross_sections():
    os.environ['OPENMC_CROSS_SECTIONS'] = os.path.abspath(CROSS_SECTION_FILE)
    openmc.config['cross_sections'] = CROSS_SECTION_FILE

get_cross_sections()

#hof (too hot): 4.038040123104094, 8.393284293720084, 3.7770687867723876, 1.6373403037006096, 7.161185801977469, 10.536799150808212, 19.63676155271593, 3.30441685661253, 18.747002613624964, 12.288070345018621, 16.510965397018822, 3.689010461088275

r1 = {{r1}}
r2 = {{r2}}
r3 = {{r3}}
r4 = {{r4}}
r5 = {{r5}}
r6 = {{r6}}
r7 = {{r7}}
r8 = {{r8}}
r9 = {{r9}}
r10 = {{r10}}
r11 = {{r11}}
r12 = {{r12}}





t_urbo = 1E-8  #treating as 0
atom_fraction_U = 1  #keep it here

FUEL_AND_urbo_OUTER_RADIUS = 0.415  # Radius of the 4 fuels plus the urbo
# Default values from the original script
FR1 = 0.20475
FR2 = 0.28956
FR3 = 0.35464
FR4 = 0.4095


def save_radius_values(fr1, fr2, fr3, fr4, thickness_urbo):
    """Saves radius values into the file."""
    gen_pop = ['standalone', 'run']  # Replace this with meaningful values if necessary

    with open(os.path.join(ROOT_DIR, RADIUS_VALUES), 'a') as file:
        line = f"{gen_pop[0]} {gen_pop[1]} {fr1} {fr2} {fr3} {fr4} {FUEL_AND_urbo_OUTER_RADIUS} {thickness_urbo} {atom_fraction_U}\n"
        file.write(line)

def get_fuel_radius_values(urbo_thickness: float):
    """Scales the fuel ring radius values according to the change of the
    urbo thickness keeping the total radius fixed"""
    fr4_new = FUEL_AND_urbo_OUTER_RADIUS - urbo_thickness
    alpha = fr4_new / FR4
    fr1 = FR1 * alpha
    fr2 = fr1 + (FR2 - FR1) * alpha
    fr3 = fr2 + (FR3 - FR2) * alpha
    fr4 = fr4_new
    return fr1, fr2, fr3, fr4

def run_calculation(r1__, r2__, r3__, r4__,
                    r5__, r6__, r7__, r8__,
                    r9__, r10__, r11__, r12__,
                    thickness_urbo=None):
    if thickness_urbo is None:
        thickness_urbo = FUEL_AND_urbo_OUTER_RADIUS - FR4
    plot_pixels = 1200
    mesh_nr = 9  # Set to 9 for quarter assembly

    # Small delta to prevent geometry overlaps
    delta = 1e-5  # cm

    def define_universe_components(r1_, r2_, r3_, r4_, r5_, r6_, r7_, r8_,
                                   r9_, r10_, r11_, r12_,
                                   no_of_pins, thickness_urbo_=FUEL_AND_urbo_OUTER_RADIUS - FR4,
                                   rho_UB2_diluted=12.71):
        # Fuel radii calculations
        fr1, fr2, fr3, fr4 = get_fuel_radius_values(thickness_urbo_)
        save_radius_values(fr1, fr2, fr3, fr4, thickness_urbo_)
        fuel_radius1 = openmc.ZCylinder(r=fr1)
        fuel_radius2 = openmc.ZCylinder(r=fr2)
        fuel_radius3 = openmc.ZCylinder(r=fr3)
        fuel_outer_radius = openmc.ZCylinder(r=fr4 - delta)  # Adjusted radius
        urbo_radius = openmc.ZCylinder(r=FUEL_AND_urbo_OUTER_RADIUS)
        # Other cylinders and planes
        clad_inner_radius = openmc.ZCylinder(r=0.4177)
        clad_outer_radius = openmc.ZCylinder(r=0.4749)
        gt_inner_radius = openmc.ZCylinder(r=0.5715)
        gt_outer_radius = openmc.ZCylinder(r=0.612)
        lower_cutback_height = openmc.ZPlane(z0=-426.0 / 2.0 + 15)
        upper_cutback_height = openmc.ZPlane(z0=+426.0 / 2.0 - 15)
        fuel_bottom = openmc.ZPlane(z0=-426.0 / 2.0, boundary_type='reflective')
        fuel_top = openmc.ZPlane(z0=+426.0 / 2.0, boundary_type='reflective')

        # Initialize materials list
        materials = []

        # Define materials with isotopic compositions
        def calculate_mass_fractions(enrichment_percent):
            # Enrichment as fraction
            enrichment = enrichment_percent / 100.0
            # Approximate U-234 content: U-234 mass fraction is about 0.01 times U-235 mass fraction
            mass_frac_U235 = enrichment
            mass_frac_U234 = mass_frac_U235 * 0.01
            mass_frac_U238 = 1.0 - mass_frac_U235 - mass_frac_U234
            # Return mass fractions
            return mass_frac_U234, mass_frac_U235, mass_frac_U238

        # Fuel materials for main axial fuel region
        # Fuel 1
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r1_)
        f_1 = openmc.Material(name='Fuel 1')
        f_1.add_nuclide('U234', mass_frac_U234)
        f_1.add_nuclide('U235', mass_frac_U235)
        f_1.add_nuclide('U238', mass_frac_U238)
        f_1.add_nuclide('O16', 2.0)
        f_1.set_density('g/cm3', 10.4)
        f_1.depletable = True
        f_1.volume = (np.pi * fuel_radius1.r ** 2) * 396 * no_of_pins
        materials.append(f_1)

        # Fuel 2
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r2_)
        f_2 = openmc.Material(name='Fuel 2')
        f_2.add_nuclide('U234', mass_frac_U234)
        f_2.add_nuclide('U235', mass_frac_U235)
        f_2.add_nuclide('U238', mass_frac_U238)
        f_2.add_nuclide('O16', 2.0)
        f_2.set_density('g/cm3', 10.4)
        f_2.depletable = True
        f_2.volume = np.pi * (fuel_radius2.r ** 2 - fuel_radius1.r ** 2) * 396 * no_of_pins
        materials.append(f_2)

        # Fuel 3
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r3_)
        f_3 = openmc.Material(name='Fuel 3')
        f_3.add_nuclide('U234', mass_frac_U234)
        f_3.add_nuclide('U235', mass_frac_U235)
        f_3.add_nuclide('U238', mass_frac_U238)
        f_3.add_nuclide('O16', 2.0)
        f_3.set_density('g/cm3', 10.4)
        f_3.depletable = True
        f_3.volume = np.pi * (fuel_radius3.r ** 2 - fuel_radius2.r ** 2) * 396 * no_of_pins
        materials.append(f_3)

        # Fuel 4
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r4_)
        f_4 = openmc.Material(name='Fuel 4')
        f_4.add_nuclide('U234', mass_frac_U234)
        f_4.add_nuclide('U235', mass_frac_U235)
        f_4.add_nuclide('U238', mass_frac_U238)
        f_4.add_nuclide('O16', 2.0)
        f_4.set_density('g/cm3', 10.4)
        f_4.depletable = True
        f_4.volume = np.pi * (fuel_outer_radius.r ** 2 - fuel_radius3.r ** 2) * 396 * no_of_pins
        materials.append(f_4)

        # Fuel materials for lower cutback region
        # Fuel 1b
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r5_)
        f_1b = openmc.Material(name='Fuel 1b')
        f_1b.add_nuclide('U234', mass_frac_U234)
        f_1b.add_nuclide('U235', mass_frac_U235)
        f_1b.add_nuclide('U238', mass_frac_U238)
        f_1b.add_nuclide('O16', 2.0)
        f_1b.set_density('g/cm3', 10.4)
        f_1b.depletable = True
        f_1b.volume = (np.pi * fuel_radius1.r ** 2) * 15 * no_of_pins
        materials.append(f_1b)

        # Fuel 2b
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r6_)
        f_2b = openmc.Material(name='Fuel 2b')
        f_2b.add_nuclide('U234', mass_frac_U234)
        f_2b.add_nuclide('U235', mass_frac_U235)
        f_2b.add_nuclide('U238', mass_frac_U238)
        f_2b.add_nuclide('O16', 2.0)
        f_2b.set_density('g/cm3', 10.4)
        f_2b.depletable = True
        f_2b.volume = np.pi * (fuel_radius2.r ** 2 - fuel_radius1.r ** 2) * 15 * no_of_pins
        materials.append(f_2b)

        # Fuel 3b
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r7_)
        f_3b = openmc.Material(name='Fuel 3b')
        f_3b.add_nuclide('U234', mass_frac_U234)
        f_3b.add_nuclide('U235', mass_frac_U235)
        f_3b.add_nuclide('U238', mass_frac_U238)
        f_3b.add_nuclide('O16', 2.0)
        f_3b.set_density('g/cm3', 10.4)
        f_3b.depletable = True
        f_3b.volume = np.pi * (fuel_radius3.r ** 2 - fuel_radius2.r ** 2) * 15 * no_of_pins
        materials.append(f_3b)

        # Fuel 4b
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r8_)
        f_4b = openmc.Material(name='Fuel 4b')
        f_4b.add_nuclide('U234', mass_frac_U234)
        f_4b.add_nuclide('U235', mass_frac_U235)
        f_4b.add_nuclide('U238', mass_frac_U238)
        f_4b.add_nuclide('O16', 2.0)
        f_4b.set_density('g/cm3', 10.4)
        f_4b.depletable = True
        f_4b.volume = np.pi * (fuel_outer_radius.r ** 2 - fuel_radius3.r ** 2) * 15 * no_of_pins
        materials.append(f_4b)

        # Fuel materials for upper cutback region
        # Fuel 1u
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r9_)
        f_1u = openmc.Material(name='Fuel 1u')
        f_1u.add_nuclide('U234', mass_frac_U234)
        f_1u.add_nuclide('U235', mass_frac_U235)
        f_1u.add_nuclide('U238', mass_frac_U238)
        f_1u.add_nuclide('O16', 2.0)
        f_1u.set_density('g/cm3', 10.4)
        f_1u.depletable = True
        f_1u.volume = (np.pi * fuel_radius1.r ** 2) * 15 * no_of_pins
        materials.append(f_1u)

        # Fuel 2u
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r10_)
        f_2u = openmc.Material(name='Fuel 2u')
        f_2u.add_nuclide('U234', mass_frac_U234)
        f_2u.add_nuclide('U235', mass_frac_U235)
        f_2u.add_nuclide('U238', mass_frac_U238)
        f_2u.add_nuclide('O16', 2.0)
        f_2u.set_density('g/cm3', 10.4)
        f_2u.depletable = True
        f_2u.volume = np.pi * (fuel_radius2.r ** 2 - fuel_radius1.r ** 2) * 15 * no_of_pins
        materials.append(f_2u)

        # Fuel 3u
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r11_)
        f_3u = openmc.Material(name='Fuel 3u')
        f_3u.add_nuclide('U234', mass_frac_U234)
        f_3u.add_nuclide('U235', mass_frac_U235)
        f_3u.add_nuclide('U238', mass_frac_U238)
        f_3u.add_nuclide('O16', 2.0)
        f_3u.set_density('g/cm3', 10.4)
        f_3u.depletable = True
        f_3u.volume = np.pi * (fuel_radius3.r ** 2 - fuel_radius2.r ** 2) * 15 * no_of_pins
        materials.append(f_3u)

        # Fuel 4u
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r12_)
        f_4u = openmc.Material(name='Fuel 4u')
        f_4u.add_nuclide('U234', mass_frac_U234)
        f_4u.add_nuclide('U235', mass_frac_U235)
        f_4u.add_nuclide('U238', mass_frac_U238)
        f_4u.add_nuclide('O16', 2.0)
        f_4u.set_density('g/cm3', 10.4)
        f_4u.depletable = True
        f_4u.volume = np.pi * (fuel_outer_radius.r ** 2 - fuel_radius3.r ** 2) * 15 * no_of_pins
        materials.append(f_4u)

        ### Urbo layer material
        UB2_diluted = openmc.Material(name='UB2_diluted')
        # Calculate mass fractions for uranium in the UB2_diluted material
        mass_frac_U234, mass_frac_U235, mass_frac_U238 = calculate_mass_fractions(r12_)
        total_uranium_mass_fraction = atom_fraction_U
        UB2_diluted.add_nuclide('U234', mass_frac_U234 * total_uranium_mass_fraction)
        UB2_diluted.add_nuclide('U235', mass_frac_U235 * total_uranium_mass_fraction)
        UB2_diluted.add_nuclide('U238', mass_frac_U238 * total_uranium_mass_fraction)
        # Add boron isotopes
        UB2_diluted.add_nuclide('B10', 0.16 * (1 - atom_fraction_U) * 2)
        UB2_diluted.add_nuclide('B11', 0.84 * (1 - atom_fraction_U) * 2)
        # Densities of UB2 and U
        density_UB2 = 12.71  # g/cm3
        density_U = 19.1  # g/cm3
        # Calculate the weighted average density
        density_UB2_diluted = (1 - atom_fraction_U) * density_UB2 + atom_fraction_U * density_U
        UB2_diluted.set_density('g/cm3', density_UB2_diluted)
        UB2_diluted.depletable = True
        UB2_diluted.volume = np.pi * (urbo_radius.r ** 2 - (fuel_outer_radius.r + delta) ** 2) * 426 * no_of_pins
        materials.append(UB2_diluted)

        zirconium = openmc.Material(name="zirconium")
        zirconium.add_element('Zr', 1.0)
        zirconium.set_density('g/cm3', 6.6)
        materials.append(zirconium)

        helium = openmc.Material(name='Helium')
        helium.add_element('He', 1.0)
        helium.set_density('g/cm3', 0.178e-3)
        materials.append(helium)

        water = openmc.Material(name='water')
        water.add_nuclide('H1', 2.0)
        water.add_nuclide('O16', 1.0)
        water.set_density('g/cm3', 1.0)
        water.add_s_alpha_beta('c_H_in_H2O')
        materials.append(water)

        # Define regions
        fuel1_region = -fuel_radius1 & +lower_cutback_height & -upper_cutback_height
        fuel2_region = +fuel_radius1 & -fuel_radius2 & +lower_cutback_height & -upper_cutback_height
        fuel3_region = +fuel_radius2 & -fuel_radius3 & +lower_cutback_height & -upper_cutback_height
        fuel4_region = +fuel_radius3 & -fuel_outer_radius & +lower_cutback_height & -upper_cutback_height

        # Cutback regions
        lower_cutback_region1 = -fuel_radius1 & -lower_cutback_height & +fuel_bottom
        lower_cutback_region2 = +fuel_radius1 & -fuel_radius2 & -lower_cutback_height & +fuel_bottom
        lower_cutback_region3 = +fuel_radius2 & -fuel_radius3 & -lower_cutback_height & +fuel_bottom
        lower_cutback_region4 = +fuel_radius3 & -fuel_outer_radius & -lower_cutback_height & +fuel_bottom

        upper_cutback_region1 = -fuel_radius1 & +upper_cutback_height & -fuel_top
        upper_cutback_region2 = +fuel_radius1 & -fuel_radius2 & +upper_cutback_height & -fuel_top
        upper_cutback_region3 = +fuel_radius2 & -fuel_radius3 & +upper_cutback_height & -fuel_top
        upper_cutback_region4 = +fuel_radius3 & -fuel_outer_radius & +upper_cutback_height & -fuel_top

        # Urbo region with adjusted boundaries to prevent overlaps
        urbo_region = +fuel_outer_radius & -urbo_radius & +fuel_bottom & -fuel_top

        # Gap and clad regions
        gap_region = +urbo_radius & -clad_inner_radius & +fuel_bottom & -fuel_top
        clad_region = +clad_inner_radius & -clad_outer_radius & +fuel_bottom & -fuel_top

        # Guide tube regions (if applicable)
        gt_water_region = -gt_inner_radius & +fuel_bottom & -fuel_top
        gt_clad_region = +gt_inner_radius & -gt_outer_radius & +fuel_bottom & -fuel_top

        # Define cells
        fuel1 = openmc.Cell(name='fuel1', fill=f_1, region=fuel1_region)
        fuel2 = openmc.Cell(name='fuel2', fill=f_2, region=fuel2_region)
        fuel3 = openmc.Cell(name='fuel3', fill=f_3, region=fuel3_region)
        fuel4 = openmc.Cell(name='fuel4', fill=f_4, region=fuel4_region)

        lower_cutback1 = openmc.Cell(name='Lower cutback 1', fill=f_1b, region=lower_cutback_region1)
        lower_cutback2 = openmc.Cell(name='Lower cutback 2', fill=f_2b, region=lower_cutback_region2)
        lower_cutback3 = openmc.Cell(name='Lower cutback 3', fill=f_3b, region=lower_cutback_region3)
        lower_cutback4 = openmc.Cell(name='Lower cutback 4', fill=f_4b, region=lower_cutback_region4)

        upper_cutback1 = openmc.Cell(name='Upper cutback 1', fill=f_1u, region=upper_cutback_region1)
        upper_cutback2 = openmc.Cell(name='Upper cutback 2', fill=f_2u, region=upper_cutback_region2)
        upper_cutback3 = openmc.Cell(name='Upper cutback 3', fill=f_3u, region=upper_cutback_region3)
        upper_cutback4 = openmc.Cell(name='Upper cutback 4', fill=f_4u, region=upper_cutback_region4)

        urbo = openmc.Cell(name='urbo coating', fill=UB2_diluted, region=urbo_region)
        gap = openmc.Cell(name='air gap', fill=helium, region=gap_region)
        clad = openmc.Cell(name='clad', fill=zirconium, region=clad_region)
        moderator = openmc.Cell(name='Moderator', fill=water,
                                region=+clad_outer_radius & -fuel_top & +fuel_bottom)

        # Return necessary components, including the materials list
        return (gt_outer_radius, fuel_bottom, fuel_top, f_1, f_2, f_3, f_4, f_1b, f_2b,
                f_3b, f_4b, f_1u, f_2u, f_3u, f_4u,
                UB2_diluted, zirconium, water, gt_water_region, gt_clad_region,
                fuel1, fuel2, fuel3, fuel4,
                lower_cutback1, lower_cutback2, lower_cutback3, lower_cutback4,
                upper_cutback1, upper_cutback2, upper_cutback3, upper_cutback4,
                gap, clad, moderator, urbo, materials)

    # Outer pins with minimal urbo thickness
    (gt_outer_radius_1, fuel_bottom_1, fuel_top_1, f_1, f_2, f_3, f_4,
     f_1b, f_2b, f_3b, f_4b, f_1u, f_2u, f_3u, f_4u,
     UB2_diluted, zirconium, water, gt_water_region, gt_clad_region,
     fuel1, fuel2, fuel3, fuel4,
     lower_cutback1, lower_cutback2, lower_cutback3, lower_cutback4,
     upper_cutback1, upper_cutback2, upper_cutback3, upper_cutback4,
     gap, clad, moderator, urbo, materials_outer) = define_universe_components(
        r1__, r2__, r3__, r4__, r5__, r6__, r7__, r8__,
        r9__, r10__, r11__, r12__, no_of_pins=53.5, thickness_urbo_=1e-6)

    # Define outer pin cell universe
    pin_cell_universe_outer = openmc.Universe(name='Fuel Pin (Outer)', cells=[
        fuel1, fuel2, fuel3, fuel4,
        lower_cutback1, lower_cutback2, lower_cutback3, lower_cutback4,
        upper_cutback1, upper_cutback2, upper_cutback3, upper_cutback4,
        urbo, gap, clad, moderator])

    # Inner pins with specified urbo thickness
    (gt_outer_radius_2, fuel_bottom_2, fuel_top_2, f_1_, f_2_, f_3_, f_4_,
     f_1b_, f_2b_, f_3b_, f_4b_, f_1u_, f_2u_, f_3u_, f_4u_,
     UB2_diluted_, zirconium_, water_, gt_water_region_, gt_clad_region_,
     fuel1_, fuel2_, fuel3_, fuel4_,
     lower_cutback1_, lower_cutback2_, lower_cutback3_, lower_cutback4_,
     upper_cutback1_, upper_cutback2_, upper_cutback3_, upper_cutback4_,
     gap_, clad_, moderator_, urbo_, materials_inner) = define_universe_components(
        r1__, r2__, r3__, r4__, r5__, r6__, r7__, r8__,
        r9__, r10__, r11__, r12__, no_of_pins=12.5, thickness_urbo_=thickness_urbo)

    # Define inner pin cell universe
    pin_cell_universe_inner = openmc.Universe(name='Fuel Pin (Inner)', cells=[
        fuel1_, fuel2_, fuel3_, fuel4_,
        lower_cutback1_, lower_cutback2_, lower_cutback3_, lower_cutback4_,
        upper_cutback1_, upper_cutback2_, upper_cutback3_, upper_cutback4_,
        urbo_, gap_, clad_, moderator_])

    # Combine materials from both outer and inner pins, avoiding duplicates
    all_materials = materials_outer + materials_inner

    # Remove duplicates based on material IDs
    unique_materials_dict = {material.id: material for material in all_materials}
    unique_materials = list(unique_materials_dict.values())

    # Create Materials object and export
    materials_file = openmc.Materials(unique_materials)
    materials_file.export_to_xml()

    # Define guide tube universe
    gt1 = openmc.Cell(name='Guide tube filled with water', fill=water, region=gt_water_region)
    gt2 = openmc.Cell(name='Guide tube clad', fill=zirconium, region=gt_clad_region)
    mode = openmc.Cell(name='Moderator', fill=water,
                       region=+gt_outer_radius_1 & -fuel_top_1 & +fuel_bottom_1)
    gt_universe = openmc.Universe(name='Guide Tube', cells=[gt1, gt2, mode])

    pitch = 1.26  # cm
    assembly = openmc.RectLattice(name='Fuel - 0BA')
    assembly.pitch = (pitch, pitch)
    assembly.lower_left = (0.0, 0.0)  # For quarter-core symmetry

    # Define the assembly universes for a 9x9 quarter assembly
    assembly.universes = [
        [pin_cell_universe_outer for _ in range(mesh_nr)],  # Row 1
        [pin_cell_universe_outer for _ in range(mesh_nr)],  # Row 2
        [pin_cell_universe_outer for _ in range(5)] + [gt_universe] + [pin_cell_universe_inner for _ in range(2)] + [gt_universe],  # Row 3
        [pin_cell_universe_outer for _ in range(3)] + [gt_universe] + [pin_cell_universe_inner for _ in range(5)],  # Row 4
        [pin_cell_universe_outer for _ in range(mesh_nr)],  # Row 5
        [pin_cell_universe_outer for _ in range(2)] + [gt_universe] + [pin_cell_universe_inner for _ in range(2)] + [gt_universe] + [pin_cell_universe_inner for _ in range(3)],  # Row 6
        [pin_cell_universe_outer for _ in range(mesh_nr)],  # Row 7
        [pin_cell_universe_outer for _ in range(mesh_nr)],  # Row 8
        [pin_cell_universe_outer for _ in range(2)] + [gt_universe] + [pin_cell_universe_inner for _ in range(2)] + [gt_universe] + [pin_cell_universe_outer for _ in range(3)],  # Row 9
    ]

    # Verify row lengths and element types
    for i, row in enumerate(assembly.universes, start=1):
        assert len(row) == mesh_nr, f"Row {i} length {len(row)} does not match mesh_nr {mesh_nr}"
        for j, universe in enumerate(row, start=1):
            assert isinstance(universe, openmc.Universe), f"Element at row {i}, column {j} is not a Universe."

    # Root cell and universe
    root_cell = openmc.Cell(name='root cell', fill=assembly)

    # Define boundary planes with symmetry considerations
    # Define boundary planes with symmetry considerations
    min_x = openmc.XPlane(x0=0.0, boundary_type='reflective')  # Symmetry plane (left)
    max_x = openmc.XPlane(x0=pitch * mesh_nr, boundary_type='reflective')  # Change to 'reflective'
    min_y = openmc.YPlane(y0=0.0, boundary_type='reflective')  # Symmetry plane (bottom)
    max_y = openmc.YPlane(y0=pitch * mesh_nr, boundary_type='reflective')  # Change to 'reflective'

    # Add boundary planes to root cell
    root_cell.region = +min_x & -max_x & +min_y & -max_y & +fuel_bottom_1 & -fuel_top_1

    # Create root Universe
    root_universe = openmc.Universe(name='root universe')
    root_universe.add_cell(root_cell)

    geometry = openmc.Geometry(root_universe)
    geometry.export_to_xml()

    # Settings
    settings = openmc.Settings()
    settings.source = openmc.Source(space=openmc.stats.Box(
        (0.0, 0.0, fuel_bottom_1.z0),
        (pitch * mesh_nr, pitch * mesh_nr, fuel_top_1.z0),
        only_fissionable=True
    ))
    settings.batches = 70
    settings.inactive = 40
    settings.particles = 5000
    settings.check_overlaps = True  # Enable overlap checking
    settings.export_to_xml()

    # Tallies
    tallies = openmc.Tallies()
    mesh = openmc.RegularMesh(mesh_id=1)
    mesh.dimension = [mesh_nr, mesh_nr]
    mesh.lower_left = [0.0, 0.0]
    mesh.width = [pitch, pitch]

    mesh_filter = openmc.MeshFilter(mesh)

    tally = openmc.Tally(name='mesh tally')
    tally.filters = [mesh_filter]
    tally.scores = ['fission']

    tallies.append(tally)
    tallies.export_to_xml()

    # Plots
    plot1 = openmc.Plot()
    plot1.filename = 'pinplot'
    plot1.origin = (pitch * mesh_nr / 2, pitch * mesh_nr / 2, 0)
    plot1.width = (pitch * mesh_nr, pitch * mesh_nr)
    plot1.pixels = (plot_pixels, plot_pixels)
    plot1.color_by = 'material'
    plot1.colors = {
        f_1: 'red', f_2: 'orange', f_3: 'yellow', f_4: 'green',
        UB2_diluted: 'maroon', water: 'blue',
        f_1b: 'lightcoral', f_2b: 'gold', f_3b: 'khaki', f_4b: 'lightgreen',
        f_1u: 'darkred', f_2u: 'chocolate', f_3u: 'olive', f_4u: 'darkgreen',
    }

    plots = openmc.Plots([plot1])
    plots.export_to_xml()
    time_step_nr = 5

    model = openmc.Model(geometry=geometry, settings=settings)
    model.export_to_xml()
    operator = openmc.deplete.CoupledOperator(model, DEPLETION_FILE)
    power = 21.75E06 / 4
    time_steps = [109.5] * time_step_nr
    integrator = openmc.deplete.CECMIntegrator(operator, time_steps, power, timestep_units='d')
    integrator.integrate()

def save_enrichments():
    """enrichments shall be saved into file."""
    try:
        gen_pop = os.path.split(os.path.abspath(os.curdir))[-1].split("_")
        gen = gen_pop[0]
        pop = gen_pop[1]
    except IndexError:
        gen = "unknown"
        pop = "unknown"
    
    with open(os.path.join(ROOT_DIR, ENRICHMENT_VALUES), 'a') as file:
        line = f"{gen} {pop} {r1} {r2} {r3} {r4} {r5} {r6} {r7} {r8} {r9} {r10} {r11} {r12} {t_urbo}\n"
        file.write(line)


save_enrichments()


# Run the calculation
run_calculation(r1, r2, r3, r4,
                r5, r6, r7, r8,
                r9, r10, r11, r12,
                thickness_urbo=t_urbo)