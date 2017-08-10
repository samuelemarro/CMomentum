// CMomentum.cpp : definisce il punto di ingresso dell'applicazione console.
//

#include "stdafx.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <functional>
#include <map>
#include <memory>
#include <unordered_map>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <sstream>

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

	std::function<Chromosome<T>(int, std::map<std::string, float>*)> initialization_ = NULL;
	std::function<std::pair<Chromosome<T>, Chromosome<T>>(std::vector<Chromosome<T>>*, std::map<std::string, float>*)> selection_ = NULL;
	std::function<void(Chromosome<T>*, Chromosome<T>*, std::map<std::string, float>*)> crossover_ = NULL;
	std::function<void(Chromosome<T>*, std::map<std::string, float>*)> mutation_ = NULL;
	std::function<void(Chromosome<T>*, const Parent<T>, float, std::map<std::string, float>*)> recombination_ = NULL;
	std::function<float(Chromosome<T>, std::map<std::string, float>*)> fitness_function_ = NULL;

	std::map<std::string, float> additional_parameters_ = std::map<std::string, float>();

	int chromosome_length_ = 64;
	int population_size_ = 300;
	float mutation_probability_ = 0.1f;
	float crossover_probability_ = 0.75f;
	float elitism_rate_ = 0.1f;

	float recombination_rate_ = 0.5f;

	int max_fitness_evaluations_ = -1;
	int max_generation_ = -1;
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
		CheckFunctions();
	}

	int RunAlgorithm() {
		CheckFunctions();
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

		while ((best_fitness < target_fitness_ || target_fitness_ == INT32_MAX) && 
			(fitness_evaluations < max_fitness_evaluations_ || max_fitness_evaluations_ == -1) && 
			(generation < max_generation_ || max_generation_ == -1)) {

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

		if (additional_parameters_.size() != 0) {
			parameters += "Additional Parameters:\n";
		}

		for (auto it = additional_parameters_.begin(); it != additional_parameters_.end(); it++) {
			parameters += "-" + it->first + ": " + std::to_string(it->second) + std::string("\n");
		}

		return parameters;
	}
private:
	void CheckFunctions() {
		if (initialization_ == NULL) {
			throw std::invalid_argument("The parameter \"initialization\" must be defined.");
		}
		if (selection_ == NULL) {
			throw std::invalid_argument("The parameter \"selection\" must be defined.");
		}
		if (crossover_ == NULL) {
			throw std::invalid_argument("The parameter \"crossover\" must be defined.");
		}
		if (mutation_ == NULL) {
			throw std::invalid_argument("The parameter \"mutation\" must be defined.");
		}
		if (recombination_ == NULL) {
			throw std::invalid_argument("The parameter \"recombination\" must be defined.");
		}
		if (fitness_function_ == NULL) {
			throw std::invalid_argument("The parameter \"fitness_function\" must be defined.");
		}
	}
	void CheckParameters() {
		if (population_size_ % 2 != 0 || population_size_ < 2) {
			throw std::invalid_argument("The parameter \"population_size_\" must be even number bigger or equal to 2.");
		}
		if (target_fitness_ == INT32_MAX && max_fitness_evaluations_ == -1 && max_generation_ == -1) {
			throw std::invalid_argument("At least one of the following stop conditions must be set: \"target_fitness_\", \"max_fitness_evaluations_\", \"max_generation_\".");
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

template<typename T>
float RunTest(GeneticAlgorithm<T> ga, int test_size, int cutoff_test, int cutoff_minimum) {

	std::vector<int> evaluations = std::vector<int>();

	for (int i = 1; i <= test_size; i++) {
		evaluations.push_back(ga.RunAlgorithm());

		if (i == cutoff_test && Median(evaluations) > cutoff_minimum) {
			break;
		}
	}

	return Median(evaluations);
}

void erase_writeln(std::string text) {
	std::cout << text;
	std::cout << std::string(text.length(), '\b');
}

template<typename T>
std::vector<std::pair<GeneticAlgorithm<T>, float>> RunParallelTests(std::vector<GeneticAlgorithm<T>> gas, int test_size, int cutoff_test, int cutoff_minimum) {

	std::vector<std::pair<GeneticAlgorithm<T>, float>> results = std::vector<std::pair<GeneticAlgorithm<T>, float>>();
	int total_tests = gas.size();
	int executed_tests = 0;
#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < gas.size(); i++) {
			GeneticAlgorithm<T> ga = gas[i];

			results.push_back(std::make_pair(ga, RunTest(ga, test_size, cutoff_test, cutoff_minimum)));
			executed_tests++;
			if (executed_tests % std::max(1, total_tests / 100) == 0) {
				erase_writeln("Progress: " + std::to_string(executed_tests * 100 / total_tests) + "%");
			}
		}
	}

	std::cout << "\n";

	return results;
}

std::vector<GeneticAlgorithm<bool>> MakeGas(bool optimised) {

	std::vector<GeneticAlgorithm<bool>> gas = std::vector<GeneticAlgorithm<bool>>();

	std::vector<int> population_sizes = { 100, 200, 300 };
	std::vector<float> mutation_probabilities = { 0.05f, 0.1f, 0.15f };
	std::vector<float> crossover_probabilities = { 0.5f, 0.6f, 0.7f, 0.8f };
	std::vector<float> elitism_rates = { 0, 0.05f, 0.1f, 0.15f };

	std::vector<int> tournament_sizes = { 2, 3, 4, 5 };
	std::vector<float> gene_mutation_rates = { 0.05f, 0.1f, 0.15f };

	for each(int population_size in population_sizes) {
		for each(float mutation_probability in mutation_probabilities) {
			for each(float crossover_probability in crossover_probabilities) {
				for each(float elitism_rate in elitism_rates) {
					for each(int tournament_size in tournament_sizes) {
						for each (float gene_mutation_rate in gene_mutation_rates) {

							if (optimised) {
								for (float recombination_rate = 0.1f; recombination_rate < 1; recombination_rate += 0.1f) {
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

									ga.recombination_rate_ = recombination_rate;

									ga.max_fitness_evaluations_ = 100000;

									gas.push_back(ga);
								}
							}
							else {
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

								ga.recombination_rate_ = 0;

								ga.max_fitness_evaluations_ = 100000;

								gas.push_back(ga);
							}
						}
					}
				}
			}
		}
	}

	return gas;
}

std::string FormatTime(std::time_t time_to_format, char* format) {
	std::ostringstream oss;

	std::tm tm = *std::localtime(&time_to_format);
	oss << std::put_time(&tm, format);

	return oss.str();
}

template<typename T>
std::pair<GeneticAlgorithm<T>, float> SelectBestConfiguration(std::vector<GeneticAlgorithm<T>> gas, int base_test_size, float test_size_increase_rate, float tournament_elimination, int cutoff_tests, int cutoff_minimum) {

	std::vector<std::pair<GeneticAlgorithm<bool>, float>> tournament = std::vector<std::pair<GeneticAlgorithm<bool>, float>>();
	int tournament_size = gas.size();
	int test_size = base_test_size;
	int tournament_number = 1;

	for each(GeneticAlgorithm<bool> ga in gas) {
		tournament.push_back(std::make_pair(ga, 0));
	}

	std::pair<GeneticAlgorithm<T>, float> empty_pair = std::make_pair(gas[0], 0);

	while (tournament_size > 0) {
		std::cout << "Running tournament n." << std::to_string(tournament_number) << ": " << std::to_string(test_size) << " tests for " << std::to_string(tournament_size) << " configurations" << std::endl;

		tournament = RunParallelTests(gas, test_size, (tournament_number == 1 ? cutoff_tests : -1), (tournament_number == 1 ? cutoff_minimum : 0));

		std::sort(tournament.begin(), tournament.end(),
			[](const auto& lhs, const auto& rhs) {
			return lhs.second < rhs.second;
		});

		tournament.resize(tournament_size, empty_pair);

		tournament_number++;
		tournament_size = std::fmin(static_cast<float>(tournament_size) * (1 - tournament_elimination), tournament_size - 1);
		test_size *= (1 + test_size_increase_rate);
	}
	return tournament[0];
}

int main(int argc, char **argv)
{
	srand(time(nullptr));

	std::string directory = argv[0];
	directory.erase(directory.find_last_of('\\') + 1);

	std::vector<GeneticAlgorithm<bool>> standard_gas = MakeGas(false);
	std::vector<GeneticAlgorithm<bool>> optimised_gas = MakeGas(true);

	int base_test_size = 20;
	int cutoff_test = 5;
	int cutoff_minimum = 15000;
	float test_size_increase_rate = 2;
	float elimination_rate = 0.9f;

	int final_test_size = 100000;

	int executed_tests = 0;
	std::cout << "Running Standard Test..." << std::endl;
	std::pair<GeneticAlgorithm<bool>, float> best_standard = SelectBestConfiguration(standard_gas, base_test_size, test_size_increase_rate, elimination_rate, cutoff_test, cutoff_minimum);
	
	float best_standard_evaluations = RunTest(best_standard.first, final_test_size, -1, 0);
	
	std::cout << "Running Optimised Test..." << std::endl;
	std::pair<GeneticAlgorithm<bool>, float> best_optimised = SelectBestConfiguration(optimised_gas, base_test_size, 2, 0.9f, cutoff_test, cutoff_minimum);

	float best_optimised_evaluations = RunTest(best_optimised.first, final_test_size, -1, 0);

	std::string finish_time = FormatTime(std::time(nullptr), "%d-%m-%y %H-%M-%S");
	std::string file_path = directory + "GA Results " + finish_time + ".txt";

	std::ofstream result_file(file_path);

	result_file << "BEST STANDARD (" << best_standard_evaluations << "):\n";
	result_file << best_standard.first.DumpParameters() << "\n";
	result_file << "BEST OPTIMISED(" << best_optimised_evaluations << "):\n";
	result_file << best_optimised.first.DumpParameters() << "\n";

	result_file.flush();
	result_file.close();

	system(("notepad.exe " + file_path).c_str());

	/*
	std::vector<std::pair<GeneticAlgorithm<bool>, float>> standard_evaluations = RunParallelTests(gas, test_size, &executed_tests, total_tests, cutoff_test, cutoff_minimum);

	std::vector<std::pair<GeneticAlgorithm<bool>, float>> optimised_evaluations = std::vector<std::pair<GeneticAlgorithm<bool>, float>>();

	for (float recombination_rate = 0.1f; recombination_rate < 1; recombination_rate += 0.1f) {
		std::vector<std::pair<GeneticAlgorithm<bool>, float>> optimised_results = RunParallelTests(gas, test_size, &executed_tests, total_tests, cutoff_test, cutoff_minimum);
		for each(std::pair<GeneticAlgorithm<bool>, float> result in optimised_results) {
			optimised_evaluations.push_back(result);
		}
	}

	std::sort(standard_evaluations.begin(), standard_evaluations.end(),
		[](const auto& lhs, const auto& rhs) {
		return lhs.second < rhs.second;
	});
	std::sort(optimised_evaluations.begin(), optimised_evaluations.end(),
		[](const auto& lhs, const auto& rhs) {
		return lhs.second < rhs.second;
	});
	*/

	/*std::string finish_time = FormatTime(std::time(nullptr), "%d-%m-%y %H-%M-%S");

	std::string file_path = directory + "GA Results " + finish_time + ".txt";

	std::ofstream result_file(file_path);

	result_file << "BEST STANDARD (" << standard_evaluations[0].second << "):\n";
	result_file << standard_evaluations[0].first.DumpParameters() << "\n";
	result_file << "BEST OPTIMISED(" << optimised_evaluations[0].second << "):\n";
	result_file << optimised_evaluations[0].first.DumpParameters() << "\n";

	result_file << "====STANDARD CONFIGURATIONS====\n\n";

	for (auto it = standard_evaluations.begin(); it != standard_evaluations.end(); ++it) {
		int index = it - standard_evaluations.begin() + 1;
		result_file << "Standard Configuration " << index << " (" << it->second << "):\n";
		result_file << it->first.DumpParameters() << "\n";
	}

	result_file << "====OPTIMISED CONFIGURATIONS====\n\n";

	for (auto it = optimised_evaluations.begin(); it != optimised_evaluations.end(); ++it) {
		int index = it - optimised_evaluations.begin() + 1;
		result_file << "Optimised Configuration " << index << " (" << it->second << "):\n";
		result_file << it->first.DumpParameters() << "\n";
	}
	
	result_file.flush();
	result_file.close();

	system(("notepad.exe " + file_path).c_str());
	*/
	system("PAUSE");
	return 0;
}

