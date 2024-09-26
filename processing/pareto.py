import pandas as pd
from scipy.spatial import ConvexHull
import numpy as np

# Load the OpenMC results data
file_path = 'keff_ppf.txt'
column_names = ['Generation', 'Individual', 'bol_keff', 'eol_keff', 'ppf_index', 'overall_ppf', 'eol_ppf']
openmc_data = pd.read_csv(file_path, delim_whitespace=True, header=None, names=column_names)

# Convert columns to numeric to avoid any type issues
for col in ['bol_keff', 'eol_keff', 'overall_ppf']:
    openmc_data[col] = pd.to_numeric(openmc_data[col], errors='coerce')

# Normalize the data columns
def normalize(data, column):
    return (data[column] - data[column].min()) / (data[column].max() - data[column].min())

openmc_data['Normalized BOL keff'] = normalize(openmc_data, 'bol_keff')
openmc_data['Normalized EOL keff'] = normalize(openmc_data, 'eol_keff')
openmc_data['Normalized Overall PPF'] = normalize(openmc_data, 'overall_ppf')

# Assign weights
weights = {'Normalized BOL keff': 1, 'Normalized EOL keff': 1.2, 'Normalized Overall PPF': -0.8}

# Calculate the weighted fitness for each individual
openmc_data['Weighted Fitness'] = (
    weights['Normalized BOL keff'] * openmc_data['Normalized BOL keff'] +
    weights['Normalized EOL keff'] * openmc_data['Normalized EOL keff'] +
    weights['Normalized Overall PPF'] * openmc_data['Normalized Overall PPF']
)

# Display the data with weighted fitness
print("Data with Weighted Fitness:")
print(openmc_data[['Generation', 'Individual', 'Weighted Fitness']].head())

# Determine the Pareto front
fitness_data = openmc_data[['Normalized BOL keff', 'Normalized EOL keff', 'Normalized Overall PPF']].values

# Compute the Pareto front using Convex Hull
hull = ConvexHull(fitness_data)

# Identify the points that are part of the Pareto front
pareto_indices = hull.vertices
pareto_front = openmc_data.iloc[pareto_indices]

# Display the Pareto front individuals
print("\nPareto Front Individuals:")
print(pareto_front[['Generation', 'Individual', 'Weighted Fitness']])

# Save the Pareto front individuals to a CSV file
pareto_front.to_csv('pareto_front_results.csv', index=False)
print("\nPareto front individuals saved to 'pareto_front_results.csv'")
