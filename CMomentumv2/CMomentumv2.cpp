// CMomentum.cpp : definisce il punto di ingresso dell'applicazione console.
//

#include "stdafx.h"
#include "stdio.h"
#include "time.h"
#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <map>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <math.h>

#include "Core.h"
#include "fast_rand.h"
#include "TestSuite.h"

void BitFlipMutation(Chromosome<bool>& chromosome, std::map<std::string, float>& additional_parameters) {
	float gene_mutation_rate = additional_parameters["gene_mutation_rate"];

	for (int i = 0; i < chromosome.genes_.size(); i++) {
		if (FastRand::RandomFloat() < gene_mutation_rate) {
			chromosome.genes_[i] = !chromosome.genes_[i];
		}
	}
}

void RealValuedMutation(Chromosome<float>& chromosome, std::map<std::string, float>& additional_parameters) {
	float gene_mutation_rate = additional_parameters["gene_mutation_rate"];
	float range_min = additional_parameters["range_min"];
	float range_max = additional_parameters["range_max"];
	float mutation_size = additional_parameters["relative_mutation_size"] * (range_max - range_min);

	for (int i = 0; i < chromosome.genes_.size(); i++) {
		if (FastRand::RandomFloat() < gene_mutation_rate) {
			float added_value = FastRand::RandomFloat(-mutation_size, +mutation_size);
			chromosome.genes_[i] = std::min(std::max(chromosome.genes_[i] + added_value, range_min), range_max);
		}
	}
}

template<typename T>
void TwoPointCrossover(Chromosome<T>& parent1, Chromosome<T>& parent2, std::map<std::string, float>& additional_parameters) {
	int pos1 = FastRand::RandomInt(parent1.genes_.size());
	int pos2 = FastRand::RandomInt(parent2.genes_.size());

	//Make sure that pos2 > pos1
	if (pos1 > pos2)
	{
		std::swap(pos1, pos2);
	}

	for (int i = pos1; i <= pos2; i++) {
		std::swap(parent1.genes_[i], parent2.genes_[i]);
	}
}

void IntermediateCrossover(Chromosome<float>& parent1, Chromosome<float>& parent2, std::map<std::string, float>& additional_parameters) {
	float ratio = additional_parameters["crossover_ratio"];

	float random = FastRand::RandomFloat(-ratio, 1 + ratio);

	for (int i = 0; i < parent1.genes_.size(); i++) {
		float parent1_gene = parent1.genes_[i];
		float parent2_gene = parent2.genes_[i];

		parent1.genes_[i] = parent1_gene * random + parent2_gene * (1 - random);
		parent2.genes_[i] = parent2_gene * random + parent1_gene * (1 - random);
	}
}

Chromosome<bool> BinaryInitialization(int length, std::map<std::string, float>& additional_parameters) {
	Chromosome<bool> c = Chromosome<bool>();
	c.genes_.reserve(length);
	for (int i = 0; i < length; i++) {
		c.genes_.push_back(FastRand::RandomInt(2));
	}
	return c;
}

Chromosome<float> UniformInitialization(int length, std::map<std::string, float>& additional_parameters) {
	float range_min = additional_parameters["range_min"];
	float range_max = additional_parameters["range_max"];

	Chromosome<float> c = Chromosome<float>();
	c.genes_.reserve(length);

	for (int i = 0; i < length; i++) {
		c.genes_.push_back(FastRand::RandomFloat(range_min, range_max));
	}
	return c;
}

template<typename T>
std::vector<int> ComputeDiff(Chromosome<T>& current, const Parent<T>& parent) {
	std::vector<int> diff = std::vector<int>();
	int genes_size = current.genes_.size();
	diff.reserve(genes_size);
	for (int i = 0; i < genes_size; i++) {
		if (current.genes_[i] != parent.genes_[i]) {
			diff.push_back(i);
		}
	}
	return diff;
}

void BinaryRecombination(Chromosome<bool>& current, const Parent<bool>& parent, float recombination_rate, std::map<std::string, float>& additional_parameters) {
	if (recombination_rate != 0) {
		std::vector<int> diff = ComputeDiff(current, parent);

		if (diff.size() == 0) {
			return;
		}

		int recombinations = FastRand::PolinomialInt(recombination_rate, diff.size());

		for (int i = 0; i < recombinations; i++) {
			int random_index = FastRand::RandomInt(diff.size());
			int diff_index = diff[random_index];

			current.genes_[diff_index] = !current.genes_[diff_index];

			diff.erase(diff.begin() + random_index);
		}
	}
}

void RealValuedRecombination(Chromosome<float>& current, const Parent<float>& parent, float recombination_rate, std::map<std::string, float>& additional_parameters) {
	if (recombination_rate != 0) {
		float range_min = additional_parameters["range_min"];
		float range_max = additional_parameters["range_max"];
		float mutation_size = additional_parameters["relative_mutation_size"] * (range_max - range_min);

		float diff_sum = 0;
		int diff_count = 0;
		std::vector<float> diff = std::vector<float>();
		for (int i = 0; i < current.genes_.size(); i++) {
			float diff_element = std::abs(current.genes_[i] - parent.genes_[i]);
			diff.push_back(diff_element);
			diff_sum += diff_element;
			if (diff_element != 0) {
				diff_count++;
			}
		}

		for (int i = 0; i < current.genes_.size(); i++) {
			if (diff[i] != 0 && FastRand::RandomFloat() < recombination_rate * diff_count * diff[i] / diff_sum) {
				float added_value = FastRand::RandomFloat(-mutation_size, +mutation_size);
				current.genes_[i] = std::min(std::max(current.genes_[i] + added_value, range_min), range_max);
			}
		}
	}
}

template<typename T>
std::pair<Chromosome<T>*, Chromosome<T>*> TournamentSelection(std::vector<std::unique_ptr<Chromosome<T>>>& population, std::map<std::string, float>& additional_parameters) {

	int tournament_size = static_cast<int>(additional_parameters["tournament_size"]);

	return std::pair<Chromosome<T>*, Chromosome<T>*>(RunTournament(population, tournament_size), RunTournament(population, tournament_size));
}

template<typename T>
Chromosome<T>* RunTournament(std::vector<std::unique_ptr<Chromosome<T>>>& population, int tournament_size) {
	int population_size = population.size();

	std::unique_ptr<Chromosome<T>>* winner = &population[FastRand::RandomInt(population_size)];

	for (int i = 1; i < tournament_size; i++) {
		int random_index = FastRand::RandomInt(population_size);
		if (population[random_index]->fitness_ > (*winner)->fitness_) {
			winner = &population[random_index];
		}
	}

	return winner->get();
}

float OneMaxFitness(Chromosome<bool>& chromosome, std::map<std::string, float>& additional_parameters) {
	float fitness = 0;
	for each(bool gene in chromosome.genes_) {
		if (!gene) {
			fitness--;
		}
	}
	return fitness;
}

float SphereFitness(Chromosome<float>& chromosome, std::map<std::string, float>& additional_parameters) {
	float fitness = 0;
	for each(float gene in chromosome.genes_) {
		fitness -= gene * gene;
	}
	return fitness;
}

float GriewankFitness(Chromosome<float>& chromosome, std::map<std::string, float>& additional_parameters) {
	float sum = 0;
	float product = 1;

	for (int i = 0; i < chromosome.genes_.size(); i++) {
		sum += chromosome.genes_[i] * chromosome.genes_[i];
		product *= std::cos(chromosome.genes_[i] / std::sqrt(i + 1));
	}

	return -(sum / 4000 - product + 1);
}

float RastriginFitness(Chromosome<float>& chromosome, std::map<std::string, float>& additional_parameters) {
	float sum = 0;
	for (int i = 0; i < chromosome.genes_.size(); i++) {
		sum += chromosome.genes_[i] * chromosome.genes_[i] - 10 * std::cosf(2 * std::_Pi * chromosome.genes_[i]);
	}
	return -(10 * chromosome.genes_.size() + sum);
}

float RosenbrockFitness(Chromosome<float>& chromosome, std::map<std::string, float>& additional_parameters) {
	float sum = 0;
	for (int i = 0; i < chromosome.genes_.size() - 1; i++) {
		sum += 100 * (chromosome.genes_[i + 1] - chromosome.genes_[i] * chromosome.genes_[i]) * (chromosome.genes_[i + 1] - chromosome.genes_[i] * chromosome.genes_[i]) +
			(chromosome.genes_[i] - 1) * (chromosome.genes_[i] - 1);
	}
	return -sum;
}

const int knapsack_volumes[5] = { 1, 2, 3, 4, 5 };
const int knapsack_values[5] = { 1, 2, 3, 4, 5 };
const int max_volume = 10;
const int best_value = 10;

float KnapsackFitness(Chromosome<bool>& chromosome, std::map<std::string, float>& additional_parameters) {

	int volume = 0;
	int value = 0;
	do
	{
		for (int i = 0; i < chromosome.genes_.size(); i++) {
			if (chromosome.genes_[i]) {
				volume += knapsack_volumes[i];
				value += knapsack_values[i];
			}
		}

		if (volume > max_volume) {
			int random_index = FastRand::RandomInt(chromosome.genes_.size());
			if (chromosome.genes_[random_index]) {
				chromosome.genes_[random_index] = false;
			}
		}
	} while (volume > max_volume);

	return value - best_value;
}

std::vector<GeneticAlgorithm<bool>> MakeBinaryGas(std::function<float(Chromosome<bool>&, std::map<std::string, float>& additional_parameters)> fitness, bool optimised, int chromosome_length, float target, int max_fitness_evaluations, int max_stagnation) {

	std::vector<GeneticAlgorithm<bool>> gas = std::vector<GeneticAlgorithm<bool>>();

	std::vector<int> population_sizes = { 100, 200 };
	std::vector<float> mutation_probabilities = { 0.05f, 0.1f, 0.15f };
	std::vector<float> crossover_probabilities = { 0.6f, 0.65f, 0.7f, 0.75f, 0.8f };
	std::vector<float> elitism_rates = { 0.05f, 0.1f, 0.15f };

	std::vector<int> tournament_sizes = { 2, 3, 4, 5 };
	std::vector<float> gene_mutation_rates = { 0.05f, 0.1f, 0.15f };
	std::vector<float> recombination_rates;
	if (optimised) {
		recombination_rates = { 0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f };
	}
	else {
		recombination_rates = { 0 };
	}

	for each(int population_size in population_sizes) {
		for each(float mutation_probability in mutation_probabilities) {
			for each(float crossover_probability in crossover_probabilities) {
				for each(float elitism_rate in elitism_rates) {
					for each(int tournament_size in tournament_sizes) {
						for each (float gene_mutation_rate in gene_mutation_rates) {

							for each (float recombination_rate in recombination_rates) {
								GeneticAlgorithm<bool> ga = GeneticAlgorithm<bool>(
									BinaryInitialization,
									TournamentSelection<bool>,
									TwoPointCrossover<bool>,
									BitFlipMutation,
									BinaryRecombination,
									fitness);

								ga.target_fitness_ = 0;

								ga.chromosome_length_ = chromosome_length;

								ga.population_size_ = population_size;
								ga.mutation_probability_ = mutation_probability;
								ga.crossover_probability_ = crossover_probability;
								ga.elitism_rate_ = elitism_rate;

								ga.additional_parameters_["tournament_size"] = tournament_size;
								ga.additional_parameters_["gene_mutation_rate"] = gene_mutation_rate;

								ga.recombination_rate_ = recombination_rate;

								ga.max_fitness_evaluations_ = max_fitness_evaluations;
								ga.max_stagnation_ = max_stagnation;

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

std::vector<GeneticAlgorithm<float>> MakeRealValuedGas(std::function<float(Chromosome<float>&, std::map<std::string, float>& additional_parameters)> fitness, bool optimised, int chromosome_length, float bound, float target, int max_fitness_evaluations, int max_stagnation) {
	std::vector<GeneticAlgorithm<float>> gas = std::vector<GeneticAlgorithm<float>>();

	std::vector<int> population_sizes = { 100, 200 };
	std::vector<float> mutation_probabilities = { 0.05f, 0.1f, 0.15f };
	std::vector<float> crossover_probabilities = { 0.6f, 0.65f, 0.7f, 0.75f, 0.8f };
	std::vector<float> elitism_rates = { 0.05f, 0.1f, 0.15f };

	std::vector<int> tournament_sizes = { 2, 3, 4, 5 };
	std::vector<float> gene_mutation_rates = { 0.05f, 0.1f, 0.15f };
	std::vector<float> relative_mutation_sizes = { 0.05f, 0.1f, 0.2f };

	std::vector<float> crossover_ratios = { 0, 0.25f };

	std::vector<float> recombination_rates;
	if (optimised) {
		recombination_rates = { 0.01f, 0.02f, 0.05f, 0.1f, 0.2f, 0.3f/*, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f, 0.9f, 1*/ };
	}
	else {
		recombination_rates = { 0 };
	}

	for each(int population_size in population_sizes) {
		for each(float mutation_probability in mutation_probabilities) {
			for each(float crossover_probability in crossover_probabilities) {
				for each(float elitism_rate in elitism_rates) {
					for each(int tournament_size in tournament_sizes) {
						for each (float gene_mutation_rate in gene_mutation_rates) {
							for each(float relative_mutation_size in relative_mutation_sizes) {
								for each(float crossover_ratio in crossover_ratios) {
									for each (float recombination_rate in recombination_rates) {
										GeneticAlgorithm<float> ga = GeneticAlgorithm<float>(
											UniformInitialization,
											TournamentSelection<float>,
											IntermediateCrossover,
											RealValuedMutation,
											RealValuedRecombination,
											fitness);

										ga.target_fitness_ = target;

										ga.chromosome_length_ = chromosome_length;

										ga.population_size_ = population_size;
										ga.mutation_probability_ = mutation_probability;
										ga.crossover_probability_ = crossover_probability;
										ga.elitism_rate_ = elitism_rate;

										ga.additional_parameters_["range_min"] = -bound;
										ga.additional_parameters_["range_max"] = bound;

										ga.additional_parameters_["tournament_size"] = tournament_size;
										ga.additional_parameters_["gene_mutation_rate"] = gene_mutation_rate;
										ga.additional_parameters_["relative_mutation_size"] = relative_mutation_size;
										ga.additional_parameters_["crossover_ratio"] = crossover_ratio;

										ga.recombination_rate_ = recombination_rate;

										ga.max_fitness_evaluations_ = max_fitness_evaluations;
										ga.max_stagnation_ = max_stagnation;

										gas.push_back(ga);
									}
								}
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

std::string FloatFormat(float x, const int decimal_digits) {
	std::stringstream ss;
	ss << std::fixed;
	ss.precision(decimal_digits);
	ss << x;
	return ss.str();
}

std::string ScientificNotation(float x, const int decimal_digits) {
	std::stringstream ss;
	ss.precision(decimal_digits);
	ss << std::scientific << x;
	return ss.str();
}

template<typename T>
void RunCompleteTest(std::string name, std::vector<GeneticAlgorithm<T>> gas, int base_test_size, float test_size_increase_rate, float elimination_rate, int success_rates_test_size, int section_size, std::string directory) {
	std::cout << "Finding Best Configuration...\n";
	GeneticAlgorithm<T> best = TestSuite<T>::SelectBestConfiguration(gas, base_test_size, test_size_increase_rate, elimination_rate);
	std::cout << "Found Best Configuration! Winner:\n";
	std::cout << best.DumpParameters() << "\n";
	std::cout << "Running Success Rate Test...\n";
	std::pair<DataPoint, std::vector<float>> result = TestSuite<T>::SuccessRatesTest(best, success_rates_test_size, section_size);

	directory.erase(directory.find_last_of('\\') + 1);

	std::string file_name = name + " " + std::to_string(best.chromosome_length_);
	if (best.additional_parameters_.find("range_max") != best.additional_parameters_.end()) {
		file_name += " " + std::to_string(best.additional_parameters_["range_max"]);
	}

	if (best.target_fitness_ != 0) {
		file_name += ScientificNotation(best.target_fitness_, 2);
	}

	file_name += std::string(" ") + (best.recombination_rate_ == 0 ? "Standard" : "Optimised");

	std::string finish_time = FormatTime(std::time(nullptr), "%d-%m-%y %H-%M-%S");

	std::string file_path = directory + file_name + " - GA Results " + finish_time + ".txt";

	std::ofstream result_file(file_path);

	result_file << "EXECUTION PARAMETERS:\n";
	result_file << "Name: " << name << "\n";
	result_file << "Base Test Size: " << std::to_string(base_test_size) << "\n";
	result_file << "Test Size Increase Rate: " << std::to_string(test_size_increase_rate) << "\n";
	result_file << "Elimination Rate: " << std::to_string(elimination_rate) << "\n";
	result_file << "Success Rates Test Size: " << std::to_string(success_rates_test_size) << "\n";
	result_file << "Section Size: " << std::to_string(section_size) << "\n\n";

	result_file << "BEST:\n";
	result_file << best.DumpParameters() << "\n\n";

	result_file << "Success Rate: " << std::to_string(std::accumulate(result.second.begin(), result.second.end(), 0.0f)) << "\n\n";

	result_file << "SUCCESSFUL EXECUTIONS STATISTICS:\n";
	result_file << result.first.ToString() << "\n\n";

	result_file << "SUCCESSFUL EXECUTIONS BREAKDOWN:\n";

	int section_counter = section_size;
	for (int i = 0; i < result.second.size(); i++) {
		result_file << std::to_string(section_counter) << ";" << std::to_string(result.second[i]) << "\n";
		section_counter += section_size;
	}

	result_file.flush();
	result_file.close();
}

struct ParameterSet {
	std::vector<float> population_sizes;
	std::vector<float> crossover_probabilities;
	std::vector<float> mutation_probabilities;
	std::vector<float> elitism_rates;
	std::vector<float> recombination_rates;
	std::map<std::string, std::vector<float>> additional_parameters;
};

template<typename T>
bool future_is_ready(std::future<T>& t) {
	return t.wait_for(std::chrono::seconds(0)) == std::future_status::ready;
}


int main(int argc, char **argv)
{
	std::string directory = argv[0];

	//TODO:
	//Automatizzare
	//Pausa

	int base_test_size = 10;
	float test_size_increase_rate = 1;
	float elimination_rate = 0.9f;

	int final_test_size = 100000;

	std::vector<GeneticAlgorithm<float>> gas = MakeRealValuedGas(RastriginFitness, false, 5, 5.12f, -1e-5f, 100000, 250);
	RunCompleteTest("Rastrigin", gas, base_test_size, test_size_increase_rate, elimination_rate, final_test_size, 100, directory);
	
	return 0;
}