import pickle
import csv

# Load the checkpoint file
checkpoint_file = "checkpoint.pkl"
with open(checkpoint_file, "rb") as cp_file:
    cp = pickle.load(cp_file)

# Prepare the output CSV file
output_file = "fitness_values.csv"

# Open the CSV file for writing
with open(output_file, mode='w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    
    # Write the header row
    csv_writer.writerow(['Generation', 'Individual', 'BOL keff', 'EOL keff', 'PPF'])

    # Iterate over the generations and their populations
    results = cp["all"]
    for i, pop in enumerate(results["populations"]):
        print(f"Generation {i} Fitness Values:")
        for j, ind in enumerate(pop):
            # Extract fitness values
            bol_keff, eol_keff, ppf = ind.fitness.values  # Adjust if these are named differently in your context
            # Write the fitness values to the CSV file
            csv_writer.writerow([i, j + 1, bol_keff, eol_keff, ppf])
            # Optionally, you can also print to confirm
            print(f"Generation {i}, Individual {j + 1}: BOL keff={bol_keff}, EOL keff={eol_keff}, PPF={ppf}")

print(f"Fitness values saved to {output_file}")
