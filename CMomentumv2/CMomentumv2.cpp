// CMomentum.cpp : definisce il punto di ingresso dell'applicazione console.
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <functional>
#include <map>
#include <memory>
#include <unordered_map>

int random_int(int min, int max) {
	return min + rand() / (RAND_MAX / (max - min) + 1);
}

float random_float() {
	return (float)rand() / (float)RAND_MAX;
}

template<typename T>
class Chromosome;

template<typename T>
struct Parent;

template<typename T>
class GeneticAlgorithm {
public:

	const std::function<Chromosome<T>(int, std::map<std::string, float>*)> initialization_ = NULL;
	const std::function<std::pair<Chromosome<T>, Chromosome<T>>(std::vector<Chromosome<T>>*, std::map<std::string, float>*)> selection_ = NULL;
	const std::function<void(Chromosome<T>*, Chromosome<T>*, std::map<std::string, float>*)> crossover_ = NULL;
	const std::function<void(Chromosome<T>*, std::map<std::string, float>*)> mutation_ = NULL;
	const std::function<void(Chromosome<T>*, const Parent<T>, float, std::map<std::string, float>*)> recombination_ = NULL;
	const std::function<float(Chromosome<T>, std::map<std::string, float>*)> fitness_function_ = NULL;

	std::map<std::string, float> additional_parameters_ = std::map<std::string, float>();

	int chromosome_length_ = 64;
	int population_size_ = 300;
	float mutation_probability_ = 0.1f;
	float crossover_probability_ = 0.75f;
	float elitism_rate_ = 0.1f;

	float recombination_rate_ = 0.5f;

	int max_fitness_evaluations_ = INT32_MAX;
	int max_generation_ = INT32_MAX;
	int target_fitness_ = INT32_MAX;

	GeneticAlgorithm(
		std::function<Chromosome<T>(int, std::map<std::string, float>*)> initialization,
		std::function<std::pair<Chromosome<T>, Chromosome<T>>(std::vector<Chromosome<T>>*, std::map<std::string, float>*)> selection,
		std::function<void(Chromosome<T>*, Chromosome<T>*, std::map<std::string, float>*)> crossover,
		std::function<void(Chromosome<T>*, std::map<std::string, float>*)> mutation,
		std::function<void(Chromosome<T>*, const Parent<T>, float, std::map<std::string, float>*)> recombination,
		std::function<float(Chromosome<T>, std::map<std::string, float>*)> fitness_function)
		: initialization_(initialization),
		selection_(selection),
		crossover_(crossover),
		mutation_(mutation),
		recombination_(recombination),
		fitness_function_(fitness_function)
	{
		if (initialization == NULL) {
			throw std::invalid_argument("The parameter \"initialization\" must be defined.");
		}
		if (selection == NULL) {
			throw std::invalid_argument("The parameter \"selection\" must be defined.");
		}
		if (crossover == NULL) {
			throw std::invalid_argument("The parameter \"crossover\" must be defined.");
		}
		if (mutation == NULL) {
			throw std::invalid_argument("The parameter \"mutation\" must be defined.");
		}
		if (recombination == NULL) {
			throw std::invalid_argument("The parameter \"recombination\" must be defined.");
		}
		if (fitness_function == NULL) {
			throw std::invalid_argument("The parameter \"fitness_function\" must be defined.");
		}
	}

	int RunAlgorithm() {
		CheckParameters();
		int generation = 1;
		int fitness_evaluations = 0;
		float best_fitness = INT32_MIN;

		//Create the initial population
		std::vector<Chromosome<T>> population = std::vector<Chromosome<T>>();

		for (int i = 0; i < population_size_; i++) {
			Chromosome<T> chromosome = initialization_(chromosome_length_, &additional_parameters_);
			chromosome.fitness_ = fitness_function_(chromosome, &additional_parameters_);
			chromosome.fitness_is_valid_ = true;
			population.push_back(chromosome);
		}

		//Sort the population by fitness in descending order
		sort(population.begin(), population.end(),
			[](const auto& lhs, const auto& rhs) {
			return lhs.fitness_ > rhs.fitness_;
		});

		while (best_fitness < target_fitness_ && fitness_evaluations < max_fitness_evaluations_ && generation < max_generation_) {

			std::vector<Chromosome<T>> offspring = std::vector<Chromosome<T>>();

			//Compute the elitism size (flooring to multiples of two) 
			int elitism_size = population_size_ * elitism_rate_;

			if (elitism_size % 2 == 1) {
				elitism_size--;
			}

			for (int i = 0; i < elitism_size; i++) {
				offspring.push_back(population[i]);
			}

			for (int i = 0; i < population_size_ - elitism_size; i += 2) {
				std::pair<Chromosome<T>, Chromosome<T>> parents = selection_(&population, &additional_parameters_);

				//Copy the parents
				Chromosome<T> child1 = parents.first;// std::shared_ptr<Chromosome<T>>(new Chromosome<T>(*parents.first.get()));
				Chromosome<T> child2 = parents.second;//std::shared_ptr<Chromosome<T>>(new Chromosome<T>(*parents.second.get()));

				if (random_float() < crossover_probability_) {
					crossover_(&child1, &child2, &additional_parameters_);

					child1.parent1_ = Parent<T>(parents.first.genes_, parents.first.fitness_);
					child1.parent2_ = Parent<T>(parents.second.genes_, parents.second.fitness_);
					child2.parent1_ = Parent<T>(parents.first.genes_, parents.first.fitness_);
					child2.parent2_ = Parent<T>(parents.second.genes_, parents.second.fitness_);

					child1.has_parents_ = true;
					child2.has_parents_ = true;

					child1.fitness_is_valid_ = false;
					child2.fitness_is_valid_ = false;
				}

				offspring.push_back(child1);
				offspring.push_back(child2);
			}

			for (auto& chromosome : offspring)
			{
				if (random_float() < mutation_probability_) {
					mutation_(&chromosome, &additional_parameters_);
					chromosome.fitness_is_valid_ = false;
				}

				if (!chromosome.fitness_is_valid_) {
					chromosome.fitness_ = fitness_function_(chromosome, &additional_parameters_);
					chromosome.fitness_is_valid_ = true;

					fitness_evaluations++;

					if (chromosome.has_parents_) {
						//If both parents are eligible, pick one randomly
						if (chromosome.fitness_ < chromosome.parent1_.fitness_ && chromosome.fitness_ < chromosome.parent2_.fitness_)
						{
							if (random_float() < 0.5f) {
								recombination_(&chromosome, chromosome.parent1_, recombination_rate_, &additional_parameters_);
							}
							else {
								recombination_(&chromosome, chromosome.parent2_, recombination_rate_, &additional_parameters_);
							}
						}
						else if (chromosome.fitness_ < chromosome.parent1_.fitness_) {
							recombination_(&chromosome, chromosome.parent1_, recombination_rate_, &additional_parameters_);
						}
						else if (chromosome.fitness_ < chromosome.parent2_.fitness_) {
							recombination_(&chromosome, chromosome.parent2_, recombination_rate_, &additional_parameters_);
						}
					}
				}
			}

			population = offspring;//è sufficiente?

			//Sort the population by fitness in descending order
			sort(population.begin(), population.end(),
				[](const auto& lhs, const auto& rhs) {
				return lhs.fitness_ > rhs.fitness_;
			});

			best_fitness = population[0].fitness_;

			generation++;
		}

		return fitness_evaluations;
	}

	std::string DumpParameters() {
		std::string parameters = "";

		parameters += std::string("Chromosome Length: ") + std::to_string(chromosome_length_) + std::string("\n");
		parameters += std::string("Population Size: ") + std::to_string(population_size_) + std::string("\n");
		parameters += std::string("Mutation Probability: ") + std::to_string(mutation_probability_) + std::string("\n");
		parameters += std::string("Crossover Probability: ") + std::to_string(crossover_probability_) + std::string("\n");
		parameters += std::string("Elitism Rate: ") + std::to_string(elitism_rate_) + std::string("\n");
		parameters += std::string("Recombination Rate: ") + std::to_string(recombination_rate_) + std::string("\n");

		for (auto it = additional_parameters_->begin(); it != additional_parameters_->end(); it++) {
			parameters += it->first + std::to_string(it->second) + std::string("\n");
		}

		return parameters;
	}
private:
	void CheckParameters() {
		if (population_size_ % 2 != 0 || population_size_ < 2) {
			throw std::invalid_argument("The parameter \"population_size_\" must be even number bigger or equal to 2.");
		}
		if (target_fitness_ == INT32_MAX && max_fitness_evaluations_ == INT32_MAX && max_generation_ == INT32_MAX) {
			throw std::invalid_argument("At least one of the following parameters must not be INT32_MAX: \"target_fitness_\", \"max_fitness_evaluations_\", \"max_generation_\".");
		}
	}
};

template<typename T>
struct Parent {
	std::vector<T> genes_;
	float fitness_ = 0;

	Parent() {

	}

	Parent(std::vector<T> genes, float fitness)
		: genes_(genes),
		fitness_(fitness) {

	}
};

template<typename T>
class Chromosome {
public:

	Parent<T> parent1_;
	Parent<T> parent2_;
	float fitness_ = 0;
	bool fitness_is_valid_ = false;
	bool has_parents_ = false;
	std::vector<T> genes_;

	Chromosome() {

	}

	Chromosome(const Chromosome& c) = default;
};

void BitFlipMutation(Chromosome<bool>* chromosome, std::map<std::string, float>* additional_parameters) {
	float gene_mutation_rate = (*additional_parameters)["gene_mutation_rate"];

	for (int i = 0; i < chromosome->genes_.size(); i++) {
		if (random_float() < gene_mutation_rate) {
			chromosome->genes_[i] = !chromosome->genes_[i];
		}
	}
}

template<typename T>
void TwoPointCrossover(Chromosome<T>* parent1, Chromosome<T>* parent2, std::map<std::string, float>* additional_parameters) {
	int pos1 = random_int(0, parent1->genes_.size());
	int pos2 = random_int(0, parent2->genes_.size());

	//Make sure that pos2 > pos1
	if (pos1 > pos2)
	{
		std::swap(pos1, pos2);
	}

	for (int i = pos1; i <= pos2; i++) {
		std::swap(parent1->genes_[i], parent2->genes_[i]);
	}
}

Chromosome<bool> BinaryInitialization(int length, std::map<std::string, float>* additional_parameters) {
	Chromosome<bool> c = Chromosome<bool>();
	for (int i = 0; i < length; i++) {
		c.genes_.push_back(random_int(0, 2));
	}
	return c;
}

void BinaryRecombinate(Chromosome<bool>* current, const Parent<bool> parent, float recombination_rate, std::map<std::string, float>* additional_parameters) {

	std::vector<int> diff = std::vector<int>();

	for (int i = 0; i < current->genes_.size(); i++) {
		if (current->genes_[i] != parent.genes_[i]) {
			diff.push_back(i);
		}
	}

	int recombinations = static_cast<int>(recombination_rate * diff.size());

	for (int i = 0; i < recombinations; i++) {
		int random_index = random_int(0, diff.size());
		int diff_index = diff[random_index];

		current->genes_[diff_index] = !current->genes_[diff_index];

		diff.erase(diff.begin() + random_index);
	}
}

template<typename T>
std::pair<Chromosome<T>, Chromosome<T>> TournamentSelection(std::vector<Chromosome<T>>* population, std::map<std::string, float>* additional_parameters) {

	int tournament_size = static_cast<int>((*additional_parameters)["tournament_size"]);

	return std::pair<Chromosome<T>, Chromosome<T>>(RunTournament(population, tournament_size), RunTournament(population, tournament_size));
}

template<typename T>
Chromosome<T> RunTournament(std::vector<Chromosome<T>>* population, int tournament_size) {
	Chromosome<T> winner;

	for (int i = 0; i < tournament_size; i++) {
		int random_index = random_int(0, population->size());
		if (i == 0 || (*population)[random_index].fitness_ > winner.fitness_) {
			winner = (*population)[random_index];
		}
	}

	return winner;
}

float OneMaxFitness(Chromosome<bool> chromosome, std::map<std::string, float>* additional_parameters) {
	float fitness = 0;
	for each(bool gene in chromosome.genes_) {
		if (!gene) {
			fitness--;
		}
	}
	return fitness;
}

float Median(std::vector<int> vector) {
	float median = 0;

	std::sort(vector.begin(), vector.end());

	if (vector.size() % 2 == 0) {
		median = (vector[(vector.size() - 1) / 2] + vector[(vector.size() - 1) / 2 + 1]) / 2;
	}
	else {
		median = vector[(vector.size() - 1) / 2];
	}

	return median;
}

float RunTest(GeneticAlgorithm<bool> ga, int test_size) {

	std::vector<int> evaluations = std::vector<int>();

	for (int i = 0; i < test_size; i++) {
		evaluations.push_back(ga.RunAlgorithm());
	}

	return Median(evaluations);
}

void erase_writeln(std::string text) {
	std::cout << text;
	std::cout << std::string(text.length(), '\b');
}

template<typename T>
std::vector<std::pair<GeneticAlgorithm<T>, float>> RunParallelTests(std::vector<GeneticAlgorithm<T>> gas, float recombination_rate, int test_size) {

	std::vector<std::pair<GeneticAlgorithm<T>, float>> results = std::vector<std::pair<GeneticAlgorithm<T>, float>>();

	int gas_size = gas.size();

	int completed_tests = 0;

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < gas.size(); i++) {
			//std::cout << i << std::endl;
			GeneticAlgorithm<T> ga = gas[i];
			ga.recombination_rate_ = recombination_rate;

			results.push_back(std::make_pair(ga, RunTest(ga, test_size)));
			completed_tests++;
			if (completed_tests % std::max(1, gas_size / 100) == 0) {
				erase_writeln("Progress: " + std::to_string(completed_tests * 100 / gas_size));
			}
		}
	}
	return results;
}

//TODO: Memory leak per additional_parameters?

std::vector<GeneticAlgorithm<bool>> MakeGas() {

	std::vector<GeneticAlgorithm<bool>> gas = std::vector<GeneticAlgorithm<bool>>();

	std::vector<int> population_sizes = { 100/*, 200, 300*/ };
	std::vector<float> mutation_probabilities = { 0.05f/*, 0.1f, 0.15f*/ };
	std::vector<float> crossover_probabilities = { 0.5f/*, 0.6f, 0.7f, 0.8f*/ };
	std::vector<float> elitism_rates = { 0, /*0.05f, 0.1f, 0.15f*/ };

	std::vector<int> tournament_sizes = { 2/*, 3, 4, 5*/ };
	std::vector<float> gene_mutation_rates = { 0.05f/*, 0.1f, 0.15f*/ };

	for each(int population_size in population_sizes) {
		for each(float mutation_probability in mutation_probabilities) {
			for each(float crossover_probability in crossover_probabilities) {
				for each(float elitism_rate in elitism_rates) {
					for each(int tournament_size in tournament_sizes) {
						for each (float gene_mutation_rate in gene_mutation_rates) {

							GeneticAlgorithm<bool> ga = GeneticAlgorithm<bool>(
								BinaryInitialization,
								TournamentSelection<bool>,
								TwoPointCrossover<bool>,
								BitFlipMutation,
								BinaryRecombinate,
								OneMaxFitness);

							ga.target_fitness_ = 0;

							ga.population_size_ = population_size;
							ga.mutation_probability_ = mutation_probability;
							ga.crossover_probability_ = crossover_probability;
							ga.elitism_rate_ = elitism_rate;
							ga.additional_parameters_["tournament_size"] = tournament_size;
							ga.additional_parameters_["gene_mutation_rate"] = gene_mutation_rate;

							gas.push_back(ga);

							//std::cout << ga.DumpParameters();
						}
					}
				}
			}
		}
	}

	return gas;
}

int main()
{
	srand(time(NULL));
	{
	std::vector<GeneticAlgorithm<bool>> gas = MakeGas();

	int test_size = 1;

	std::vector<std::pair<GeneticAlgorithm<bool>, float>> normal_evaluations = RunParallelTests(gas, 0, test_size);
	std::vector<std::pair<GeneticAlgorithm<bool>, float>> optimised_evaluations = std::vector<std::pair<GeneticAlgorithm<bool>, float>>();
	for (float recombination_rate = 0.1f; recombination_rate < 0.2f/*1*/; recombination_rate += 0.1f) {
		std::vector<std::pair<GeneticAlgorithm<bool>, float>> optimised_results = RunParallelTests(gas, recombination_rate, test_size);
		for each(std::pair<GeneticAlgorithm<bool>, float> result in optimised_results) {
			optimised_evaluations.push_back(result);
		}
	}
	

	//std::sort(standard_evaluations.begin(), standard_evaluations.end());

	//float standard_median = Median(standard_evaluations);
	}
	//_CrtDumpMemoryLeaks();
	system("PAUSE");
	return 0;
}

