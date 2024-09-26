import os
import pandas as pd
from matplotlib import pyplot as plt

def read_fixed_structure_file(file_path, column_names):
    """
    Reads a fixed structure text file into a pandas DataFrame.
    """
    df = pd.read_csv(file_path, sep=" ", header=None)
    df.columns = column_names
    return df

def plot_line(df, x_col, y_col, xlabel=None, ylabel=None, title=None, label=None):
    """
    Plots a 2D line plot for the specified columns in the DataFrame.
    """
    plt.plot(df[x_col], df[y_col], "-", alpha=0.7, label=label if label else y_col)
    plt.xlabel(xlabel if xlabel else x_col)
    plt.ylabel(ylabel if ylabel else y_col)
    plt.title(title if title else f'Line Plot of {x_col} vs {y_col}')

ROOT_DIR = "."
ENRICHMENT_VALUES = "enrichments.txt"

# Parse the enrichment file
enrichment_file_path = os.path.join(ROOT_DIR, ENRICHMENT_VALUES)
enrichment_column_names = ['gen', 'pop'] + [f"r{i}" for i in range(1, 14)]
enrichment = read_fixed_structure_file(enrichment_file_path, enrichment_column_names)
print(enrichment.head())

def plot_enrichments(enrichment, start_index, end_index, filename, labels, legend_position='upper left'):
    plt.figure()
    enrichment_averaged = enrichment.groupby('gen').mean().reset_index()
    x_col = 'gen'
    for idx, r_index in enumerate(range(start_index, end_index + 1)):
        plot_line(enrichment_averaged, x_col, f'r{r_index}', xlabel='Generation',
                  ylabel='Enrichment Values', title='Fuel Enrichment Values Over Generations',
                  label=labels[idx])
    legend = plt.legend(loc=legend_position, bbox_to_anchor=(1, 0.3), frameon=False, fontsize='small', fancybox=True, framealpha=0.5)
    for text in legend.get_texts():
        text.set_alpha(0.4)  # Adjust text transparency
    plt.xticks(range(0, int(enrichment['gen'].max()) + 1, 15))  # Show every 15 generations on x-axis
    plt.savefig(filename, dpi=300, bbox_inches='tight')



# Plot fuel zones 1 to 4 with specific labels
plot_enrichments(enrichment, 1, 4, "enrichments_zones_1_4.png", 
                 labels=["Innermost Pellet Zone", "Pellet Zone 2", "Pellet Zone 3", "Outermost Pellet Zone"])

# Plot fuel zones 5 to 8 with specific labels, legend in the upper right corner
plot_enrichments(enrichment, 5, 8, "enrichments_zones_5_8.png", 
                 labels=["Lower Innermost Zone", "Lower Zone 2", "Lower Zone 3", "Lower Outermost Zone"], 
                 legend_position='lower right')
