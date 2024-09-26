import os
import re
import openmc
import numpy as np
import openmc.deplete
import matplotlib.pyplot as plt

ROOT_DIR = "/home/hpcatta/submit1/beyond_100"
RESULT_FILE_NAME = "keff_ppf.txt"

def get_simulation_files():
    def extract_number_from_filename(filepath):
        # Extract the numeric part of the filenames
        match = re.search(r'_n(\d+)\.h5$', filepath)
        if match:
            return int(match.group(1))
        return float('inf')  # In case of missing number push it to the end

    def sort_file_paths(file_paths):
        return sorted(file_paths, key=extract_number_from_filename)

    all_files = os.listdir(os.curdir)
    simulation_files = [os.path.join(os.path.abspath(os.curdir), f) for f in all_files if
                        f.startswith("openmc_simulation_") and f.endswith(".h5")]
    simulation_files = sort_file_paths(simulation_files)
    return simulation_files

def get_keff_ppf(simulation_file: str, mesh_nr = 17):
    sp = openmc.StatePoint(simulation_file)
    tally = sp.get_tally(scores=['fission'])
    tally.mean.shape = (mesh_nr, mesh_nr)
    a = np.array(tally.mean)
    average = np.sum(a) / np.count_nonzero(a)
    normalized_array = a / average
    ppf = np.max(normalized_array)
    keff = sp.keff.nominal_value
    return keff, ppf

def extract_eol_overall_ppf_eol_keff(simulation_files):
    ppfs = []
    keff_values = []
    bol_keff, eol_keff, eol_ppf = 0., 0., 0.
    if len(simulation_files) < 1:
        return bol_keff, eol_keff, 0, 0, eol_ppf, keff_values
    bol_keff, _ = get_keff_ppf(simulation_files[0])
    for file in simulation_files:
        keff, ppf = get_keff_ppf(file)
        keff_values.append(keff)  # Store keff values
        ppfs.append(ppf)
        eol_keff, eol_ppf = keff, ppf

    ppfs = np.array(ppfs)
    overall_ppf_index = np.argmax(ppfs)
    return bol_keff, eol_keff, overall_ppf_index, ppfs[overall_ppf_index], eol_ppf, keff_values

def save_quantities(bol_keff, eol_keff, overall_ppf_index, overall_ppf, eol_ppf):
    gen_pop = os.path.split(os.path.abspath(os.curdir))[-1].split("_")
    with open(os.path.join(ROOT_DIR, RESULT_FILE_NAME), 'a') as file:
        line = f"{gen_pop[0]} {gen_pop[1]} {bol_keff} {eol_keff} {overall_ppf_index} {overall_ppf} {eol_ppf}\n"
        file.write(line)

# def plot_keff_evolution(keff_values, save_path='keff_evolution.png'):
#     plt.figure(figsize=(8, 6))
#     plt.plot(range(len(keff_values)), keff_values, marker='o', linestyle='-', label='keff Evolution')
#     plt.xlabel('Depletion Stage')
#     plt.ylabel('keff')
#     plt.title('keff Evolution Over Time')
#     plt.legend()
#     plt.grid(True)
#     plt.savefig(save_path)  # Save the plot to a file
#     plt.show()  # Optionally, show the plot if running in an environment that supports it

# Main execution
simulation_files = get_simulation_files()

bol_keff, eol_keff, overall_ppf_index, overall_ppf, eol_ppf, keff_values = extract_eol_overall_ppf_eol_keff(simulation_files)
save_quantities(bol_keff, eol_keff, overall_ppf_index, overall_ppf, eol_ppf)

# Plot the keff evolution
# plot_keff_evolution(keff_values)

print({"bol_keff": bol_keff, "eol_keff": eol_keff, "overall_ppf": overall_ppf})
