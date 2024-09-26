import pickle
import numpy as np
import matplotlib.pyplot as plt

# Load the checkpoint file
checkpoint_file = "checkpoint.pkl"
with open(checkpoint_file, "rb") as cp_file:
    cp = pickle.load(cp_file)

# Extract results
results = cp["all"]
generations = list(range(len(results["populations"])))
avg_fitness_obj1, avg_fitness_obj2, avg_fitness_obj3 = [], [], []
avg_overall_fitness = []

for pop in results["populations"]:
    fitness_obj1 = [ind.fitness.values[0] for ind in pop]
    fitness_obj2 = [ind.fitness.values[1] for ind in pop]
    fitness_obj3 = [ind.fitness.values[2] for ind in pop]

    avg_fitness_obj1.append(np.mean(fitness_obj1))
    avg_fitness_obj2.append(np.mean(fitness_obj2))
    avg_fitness_obj3.append(np.mean(fitness_obj3))

    # Calculate overall fitness
    overall_fitness = [f1 + f2 + f3 for f1, f2, f3 in zip(fitness_obj1, fitness_obj2, fitness_obj3)]
    avg_overall_fitness.append(np.mean(overall_fitness))

# Plot separate objective fitnesses
plt.figure(figsize=(10, 6))
plt.plot(generations, avg_fitness_obj1, label='BoL Keff')
plt.plot(generations, avg_fitness_obj2, label='EoL Keff')
plt.plot(generations, avg_fitness_obj3, label='PPF')
plt.title('Convergence of Objectives')
plt.xlabel('Generation')
plt.ylabel('Average Fitness')
plt.legend()
plt.grid(True)
plt.savefig('convergence_separate_objectives.png')

# Plot overall fitness
plt.figure(figsize=(10, 6))
plt.plot(generations, avg_overall_fitness, label='Overall Fitness', color='purple')
plt.title('Overall Fitness Convergence')
plt.xlabel('Generation')
plt.ylabel('Average Overall Fitness')
plt.legend()
plt.grid(True)
plt.savefig('convergence_overall_fitness.png')
