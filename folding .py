import random
import numpy as np
import matplotlib.pyplot as plt

# Define the amino acid sequence of the protein
amino_acids = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
protein_length = 20
protein_sequence = [random.choice(amino_acids) for i in range(protein_length)]

# Define the chromosome representation
chromosome_length = protein_length * 3
chromosome_bounds = np.array([[0, 360]]*chromosome_length)

# Define the fitness function
def fitness_function(chromosome):
    protein_structure = decode_chromosome(chromosome)
    # Calculate the energy of the protein structure using a scoring function
    energy = calculate_energy(protein_structure)
    return -energy

# Define the genetic operators
def selection(population):
    # Select two individuals randomly from the population
    ind1 = random.choice(population)
    ind2 = random.choice(population)
    # Return the fitter individual
    return ind1 if fitness_function(ind1) > fitness_function(ind2) else ind2

def crossover(parent1, parent2):
    # Perform single-point crossover
    crossover_point = random.randint(1, chromosome_length - 1)
    child1 = np.concatenate((parent1[:crossover_point], parent2[crossover_point:]))
    child2 = np.concatenate((parent2[:crossover_point], parent1[crossover_point:]))
    return child1, child2

def mutation(chromosome):
    # Mutate each gene with a probability of 0.1
    mutated_chromosome = chromosome.copy()
    for i in range(chromosome_length):
        if random.random() < 0.1:
            mutated_chromosome[i] = random.uniform(chromosome_bounds[i][0], chromosome_bounds[i][1])
    return mutated_chromosome

# Define the decode function
def decode_chromosome(chromosome):
    # Convert the chromosome to a protein structure
    protein_structure = []
    for i in range(protein_length):
        phi, psi, omega = chromosome[i*3:(i+1)*3]
        amino_acid = protein_sequence[i]
        residue = {'phi': phi, 'psi': psi, 'omega': omega, 'amino_acid': amino_acid}
        protein_structure.append(residue)
    return protein_structure

# Define the calculate_energy function
def calculate_energy(protein_structure):
    # Calculate the energy of the protein structure using a simple scoring function
    energy = 0
    for i in range(len(protein_structure)):
        phi = protein_structure[i]['phi']
        psi = protein_structure[i]['psi']
        omega = protein_structure[i]['omega']
        energy += np.sin(phi) + np.cos(psi) + 0.5 * np.sin(2 * omega)
    return energy

# Define the main function
def main():
    pop_size = 50
    num_generations = 100

    # Initialize the population
    population = []
    for i in range(pop_size):
        chromosome = np.random.uniform(chromosome_bounds[:, 0], chromosome_bounds[:, 1], chromosome_length)
        population.append(chromosome)

    # Record the best fitness value at each generation for visualization
    best_fitness_values = []

    for i in range(num_generations):
        # Apply selection to generate parents
        parents = [selection(population) for j in range(pop_size)]
        # Apply crossover to generate offspring
        offspring = [crossover(parents[j], parents[j+1]) for j in range(0, pop_size, 2)]
        offspring = [child for pair in offspring for child in pair]
        # Apply mutation to the offspring
        offspring = [mutation(child) for child in offspring]
        # Evaluate the fitness of the offspring
        offspring_fitness = [fitness_function(child) for child in offspring]
        # Combine the population and the offspring
        combined_population = population + offspring
        combined_fitness = [fitness_function(child) for child in combined_population]
        # Select the fittest individuals for the next generation
        population = [combined_population[j] for j in sorted(range(len(combined_fitness)), key=lambda k: combined_fitness[k], reverse=True)[:pop_size]]
        # Record the best fitness value at each generation for visualization
        best_fitness_values.append(-fitness_function(population[0]))

    # Decode the fittest chromosome to obtain the protein structure
    fittest_chromosome = max(population, key=fitness_function)
    fittest_protein_structure = decode_chromosome(fittest_chromosome)

    # Plot the best fitness value at each generation
    plt.plot(best_fitness_values)
    plt.xlabel('Generation')
    plt.ylabel('Fitness')
    plt.title('Protein Folding Using Genetic Algorithms')
    plt.show()

    return fittest_protein_structure

if __name__ == '__main__':
    main()