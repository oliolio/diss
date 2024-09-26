import os
import pandas as pd
import matplotlib.pyplot as plt

# Define file paths
ROOT_DIR = "."
RESULT_FILE_NAME = "keff_ppf.txt"

# Function to read keff and ppf data
def read_keff_ppf_file(file_path):
    df = pd.read_csv(file_path, delim_whitespace=True, header=None)
    keffppf_column_names = ['gen', 'pop', 'bol_keff', 'eol_keff', 'overall_ppf_index', 'overall_ppf', 'eol_ppf']
    df.columns = keffppf_column_names
    return df

# Reading the file
keffppf = read_keff_ppf_file(os.path.join(ROOT_DIR, RESULT_FILE_NAME))
print(keffppf.head())

# Plotting PPFs
def plot_ppf(keffppf):
    plt.figure()
    plt.plot(keffppf['gen'], keffppf['eol_ppf'], 'o', alpha=0.5, label='EOL PPF')
    plt.plot(keffppf['gen'], keffppf['overall_ppf'], 'x', alpha=0.5, label='Overall PPF')
    plt.xlabel('Generation')
    plt.ylabel('Power Peaking Factor (PPF)')
    plt.title('EOL and Overall PPF over Generations')
    plt.legend()
    plt.savefig("eol_overall_ppf.png", dpi=300, bbox_inches='tight')
    plt.show()

# Plot the graph
plot_ppf(keffppf)
