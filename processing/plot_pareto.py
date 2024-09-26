import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Load Pareto front results
pareto_results = pd.read_csv('pareto_front_results.csv')
# Extract normalized values
bol_keff = pareto_results['Normalized BOL keff']
eol_keff = pareto_results['Normalized EOL keff']
overall_ppf = pareto_results['Normalized Overall PPF']

# Plotting the 3D scatter plot
fig = plt.figure(figsize=(10, 7))
ax = fig.add_subplot(111, projection='3d')

# Scatter plot of Pareto front solutions
sc = ax.scatter(bol_keff, eol_keff, overall_ppf, c='blue', marker='o', s=50)

# Label the axes
ax.set_xlabel('Normalized BOL keff')
ax.set_ylabel('Normalized EOL keff')
ax.set_zlabel('Normalized Overall PPF')
ax.set_title('Pareto Front: Multi-Objective Optimization Results')

plt.savefig('pareto_front_plot.png', dpi=300)

print("Pareto front plot saved as 'pareto_front_plot.png'.")