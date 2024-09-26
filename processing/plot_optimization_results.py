import os
import pandas as pd
from matplotlib import pyplot as plt


def read_fixed_structure_file(file_path, column_names):
    """
    Reads a fixed structure text file into a pandas DataFrame.

    Parameters:
    - file_path: str, path to the input file
    - column_names: list of str, list of column names to assign to the DataFrame

    Returns:
    - df: pandas DataFrame with the parsed data and renamed columns
    """
    df = pd.read_csv(file_path, sep=" ", header=None)
    df.columns = column_names
    return df


def plot_scatter(df, x_col, y_col, xlabel=None, ylabel=None, title=None):
    """
    Plots a 2D scatter plot for the specified columns in the DataFrame.

    Parameters:
    - df: pandas DataFrame, the data to plot
    - x_col: str, column name to use for the x-axis
    - y_col: str, column name to use for the y-axis
    - xlabel: str, label for the x-axis (optional)
    - ylabel: str, label for the y-axis (optional)
    - title: str, title for the plot (optional)
    """
    plt.plot(df[x_col], df[y_col], "o", alpha=0.5, label=y_col)
    # Set labels and title if provided
    plt.xlabel(xlabel if xlabel else x_col)
    plt.ylabel(ylabel if ylabel else y_col)
    plt.title(title if title else f'Scatter Plot of {x_col} vs {y_col}')


# ROOT_DIR = "<path to optimization root>"
ROOT_DIR = "."
RESULT_FILE_NAME = "keff_ppf.txt"
ENRICHMENT_VALUES = "enrichments.txt"
RADIUS_VALUES = "fuel_radius.txt"

# Parse the optimization outputs
keffppf_file_path = os.path.join(ROOT_DIR, RESULT_FILE_NAME)
keffppf_column_names = ['gen', 'pop', 'bol_keff', 'eol_keff', 'overall_ppf_index', 'overall_ppf', 'eol_ppf']
keffppf = read_fixed_structure_file(keffppf_file_path, keffppf_column_names)
print(keffppf.head())

enrichment_file_path = os.path.join(ROOT_DIR, ENRICHMENT_VALUES)
enrichment_column_names = ['gen', 'pop'] + [f"r{i}" for i in range(13)]
enrichment = read_fixed_structure_file(enrichment_file_path, enrichment_column_names)
print(enrichment.head())


#0.20475 0.28956 0.35464 0.4095 0.415 0.005500000000000005
radius_file_path = os.path.join(ROOT_DIR, RADIUS_VALUES)
radius_column_names = ['gen', 'pop'] + [f"fuel_r{i}" for i in range(1, 6)] + ["ifba_thickness"]
radius = read_fixed_structure_file(radius_file_path, radius_column_names)
print(radius.head())

# Calculate the average value over all the populations within a generation
keffppf = keffppf.groupby('gen').mean().reset_index()


# ---------------------------------------------------------------------
def plot_eol_keff(keffppf):
    plt.figure()
    x_col = 'gen'
    y_col = 'eol_keff'
    plot_scatter(keffppf, x_col, y_col, title="End of Life Keff. over the generations")
    plt.legend()
    # plt.show()
    plt.savefig("eol_keff.png", dpi=300, bbox_inches='tight')


def plot_ppf(keffppf):
    plt.figure()
    x_col = 'gen'
    y_col = 'eol_ppf'
    y_col2 = 'overall_ppf'
    plot_scatter(keffppf, x_col, y_col)
    plot_scatter(keffppf, x_col, y_col2, ylabel="PPF", title="End of Life and overall highest PPF.")

    plt.legend()
    # plt.show()
    plt.savefig("eol_overall_ppf.png", dpi=300, bbox_inches='tight')


def plot_enrichments(enrichment):
    plt.figure()
    enrichment_averaged = enrichment.groupby('gen').mean().reset_index()
    x_col = 'gen'
    for r_index in [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11][:12]:
        plot_scatter(enrichment_averaged, x_col, f'r{r_index}', xlabel=x_col,
                     ylabel=r'Enrichments [$\%$]', title='Fuel enrichments over the optimization')

    plt.legend()
    # plt.show()
    plt.savefig("enrichments.png", dpi=300, bbox_inches='tight')


def plot_ifba_thickness(enrichment):
    plt.figure()
    enrichment_averaged = enrichment.groupby('gen').mean().reset_index()
    x_col = 'gen'
    plot_scatter(enrichment_averaged, x_col, f'r{12}', xlabel=x_col, ylabel=r'Thickness of ifba layer [cm]',
                 title='-')

    plt.legend()
    # plt.show()
    plt.savefig("ifba_thickness.png", dpi=300, bbox_inches='tight')


def plot_radius(radius):
    plt.figure()
    radius = radius.iloc[1:][::2]
    radius_averaged = radius.groupby('gen').mean().reset_index()
    x_col = 'gen'
    for r_index in range(1, 6):
        plot_scatter(radius_averaged, x_col, f'fuel_r{r_index}', xlabel=x_col,
                     ylabel=r'Enrichments [$\%$]', title='Fuel enrichments over the optimization')

    plt.legend()
    # plt.show()
    plt.savefig("fuel_radius.png", dpi=300, bbox_inches='tight')


plot_eol_keff(keffppf)
plot_ppf(keffppf)
plot_enrichments(enrichment)
plot_ifba_thickness(enrichment)
# Plot the fuel radius over the generations
plot_radius(radius)
