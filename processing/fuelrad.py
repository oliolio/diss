import os
import pandas as pd
import matplotlib.pyplot as plt

# Function to read fixed structure files
def read_fixed_structure_file(file_path, column_names):
    df = pd.read_csv(file_path, sep=" ", header=None)
    df.columns = column_names
    return df

# Function for plotting line charts
def plot_line(df, x_col, y_col, xlabel=None, ylabel=None, title=None, label=None):
    plt.plot(df[x_col], df[y_col], "-", alpha=0.7, label=label if label else y_col)
    plt.xlabel(xlabel if xlabel else x_col)
    plt.ylabel(ylabel if ylabel else y_col)
    plt.title(title if title else f'Line Plot of {x_col} vs {y_col}')

# Function for dual-axis plotting
def plot_dual_axis(df, x_col, y_col1, y_col2, xlabel=None, ylabel1=None, ylabel2=None, title=None, label1=None, label2=None):
    fig, ax1 = plt.subplots()
    ax1.plot(df[x_col], df[y_col1], "-", alpha=0.7, color="blue", label=label1 if label1 else y_col1)
    ax1.set_xlabel(xlabel if xlabel else x_col)
    ax1.set_ylabel(ylabel1 if ylabel1 else y_col1, color="blue")
    ax1.tick_params(axis='y', labelcolor="blue")
    ax2 = ax1.twinx()
    ax2.plot(df[x_col], df[y_col2], "-", alpha=0.7, color="green", label=label2 if label2 else y_col2)
    ax2.set_ylabel(ylabel2 if ylabel2 else y_col2, color="green")
    ax2.tick_params(axis='y', labelcolor="green")
    plt.title(title if title else f'{y_col1} and {y_col2} vs {x_col}')
    fig.tight_layout()
    plt.savefig("t_urbo_vs_atom_fraction_U_dual_axis.png", dpi=300, bbox_inches='tight')

# Function to read the radius file
def read_radius_file(file_path):
    before_change, after_change = None, None
    try:
        before_change = pd.read_csv(
            file_path, delim_whitespace=True, header=None, nrows=846, 
            names=['gen', 'pop', 'fuel_r1', 'fuel_r2', 'fuel_r3', 'fuel_r4', 'fuel_r5', 't_urbo']
        )
    except pd.errors.ParserError as e:
        print(f"Error reading lines before the change in structure: {e}")

    try:
        after_change = pd.read_csv(
            file_path, delim_whitespace=True, header=None, skiprows=847, 
            names=['gen', 'pop', 'fuel_r1', 'fuel_r2', 'fuel_r3', 'fuel_r4', 'fuel_r5', 't_urbo', 'atom_fraction_U']
        )
    except pd.errors.ParserError as e:
        print(f"Error reading lines after the change in structure: {e}")

    if before_change is not None and after_change is not None:
        df_combined = pd.concat([before_change, after_change], ignore_index=True)
    elif before_change is not None:
        df_combined = before_change
    elif after_change is not None:
        df_combined = after_change
    else:
        print("Failed to read any data.")
        return None

    df_combined = df_combined[df_combined['t_urbo'] != 1e-6].copy()
    df_combined.loc[:, 't_urbo'] *= 1000  # Convert to microns

    return df_combined

# Example usage for UB2 Radius vs U fraction
file_path = "fuel_radius.txt"
radius = read_radius_file(file_path)

if radius is not None:
    radius_filtered = radius.iloc[1:][::2]
    radius_averaged = radius_filtered.groupby('gen').mean().reset_index()
    plot_dual_axis(radius_averaged, 'gen', 't_urbo', 'atom_fraction_U', xlabel='Generation',
                   ylabel1='UB2 Radius (µm)', ylabel2='Atom Fraction U:B', title='UB2 Layer Radius vs U atom fraction in UB2',
                   label1='UB2 Radius (µm)', label2='Fraction of U in UB2')

# Plot functions for EOL keff and EOL PPF, and fuel enrichment plots here
# Function to plot EOL keff
def plot_eol_keff(df):
    plt.figure()
    plt.plot(df['gen'], df['eol_keff'], "-", alpha=0.7, label='EOL keff')
    plt.xlabel('Generation')
    plt.ylabel('End of Life Keff')
    plt.title('End of Life Keff over Generations')
    plt.legend()
    plt.savefig("eol_keff.png", dpi=300, bbox_inches='tight')

# Function to plot EOL PPF
def plot_ppf(df):
    plt.figure()
    plt.plot(df['gen'], df['eol_ppf'], "-", alpha=0.7, label='EOL PPF')
    plt.plot(df['gen'], df['overall_ppf'], "-", alpha=0.7, label='Overall PPF')
    plt.xlabel('Generation')
    plt.ylabel('PPF')
    plt.title('End of Life and Overall PPF over Generations')
    plt.legend()
    plt.savefig("eol_overall_ppf.png", dpi=300, bbox_inches='tight')

# Function to plot Fuel Enrichments
def plot_enrichments(enrichment, start_index, end_index, filename, labels):
    plt.figure()
    enrichment_averaged = enrichment.groupby('gen').mean().reset_index()
    x_col = 'gen'
    for idx, r_index in enumerate(range(start_index, end_index + 1)):
        plot_line(enrichment_averaged, x_col, f'r{r_index}', xlabel='Generation',
                  ylabel='Enrichment Values [%]', title='Fuel Enrichment Values Over Generations',
                  label=labels[idx])
    plt.legend()
    plt.xticks(range(0, int(enrichment['gen'].max()) + 1, 15))
    plt.savefig(filename, dpi=300, bbox_inches='tight')

# Read the keff and PPF data
keffppf_file_path = "keff_ppf.txt"
keffppf_column_names = ['gen', 'pop', 'bol_keff', 'eol_keff', 'overall_ppf_index', 'overall_ppf', 'eol_ppf']
keffppf = read_fixed_structure_file(keffppf_file_path, keffppf_column_names)
plot_eol_keff(keffppf)
plot_ppf(keffppf)

# Read the enrichment data
enrichment_file_path = "enrichments.txt"
enrichment_column_names = ['gen', 'pop'] + [f"r{i}" for i in range(1, 14)]
enrichment = read_fixed_structure_file(enrichment_file_path, enrichment_column_names)

# Plot fuel zones 1 to 4
plot_enrichments(enrichment, 1, 4, "enrichments_zones_1_4.png", 
                 labels=["Innermost Pellet Zone", "Pellet Zone 2", "Pellet Zone 3", "Outermost Pellet Zone"])

# Plot fuel zones 5 to 8
plot_enrichments(enrichment, 5, 8, "enrichments_zones_5_8.png", 
                 labels=["Lower Innermost Pellet Zone", "Lower Pellet Zone 2", "Lower Pellet Zone 3", "Lower Outermost Pellet Zone"])



def read_radius_file(file_path):
    """
    Reads the radius file and handles different column structures before and after the change.
    """
    before_change, after_change = None, None

    # Attempt to read lines before the change
    try:
        before_change = pd.read_csv(
            file_path, delim_whitespace=True, header=None, nrows=847, 
            names=['gen', 'pop', 'fuel_r1', 'fuel_r2', 'fuel_r3', 'fuel_r4', 'fuel_r5', 't_urbo'],
            engine='python'
        )
    except pd.errors.ParserError as e:
        print(f"Error reading lines before the change in structure: {e}")

    # Attempt to read lines after the change
    try:
        after_change = pd.read_csv(
            file_path, delim_whitespace=True, header=None, skiprows=847, 
            names=['gen', 'pop', 'fuel_r1', 'fuel_r2', 'fuel_r3', 'fuel_r4', 'fuel_r5', 't_urbo', 'atom_fraction_U'],
            engine='python'
        )
    except pd.errors.ParserError as e:
        print(f"Error reading lines after the change in structure: {e}")

    # Check if either dataframe is None before concatenation
    if before_change is not None and after_change is not None:
        df_combined = pd.concat([before_change, after_change], ignore_index=True)
    elif before_change is not None:
        df_combined = before_change
    elif after_change is not None:
        df_combined = after_change
    else:
        print("Failed to read any data.")
        return None

    # Filter out rows where t_urbo is 1e-6
    df_combined = df_combined[df_combined['t_urbo'] != 1e-6]
    
    return df_combined

def plot_fuel_radius(radius):
    plt.figure()
    radius_averaged = radius.groupby('gen').mean().reset_index()
    x_col = 'gen'
    for r_index in range(1, 5):
        plt.plot(radius_averaged[x_col], radius_averaged[f'fuel_r{r_index}'] * 1000, "-", alpha=0.7, label=f'Zone {r_index} Radius (µm)')
    plt.plot(radius_averaged[x_col], radius_averaged['t_urbo'] * 1000, "-", alpha=0.7, label='UB2 Layer (µm)')
    plt.xlabel('Generation')
    plt.ylabel('Radius (µm)')
    plt.title('Fuel Radii vs Generation')
    plt.legend()
    plt.savefig("fuel_radius_plot.png", dpi=300, bbox_inches='tight')

# Example usage
file_path = "fuel_radius.txt"
radius = read_radius_file(file_path)

if radius is not None:
    plot_fuel_radius(radius)
