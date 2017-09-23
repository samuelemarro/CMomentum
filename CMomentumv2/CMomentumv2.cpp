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
		std::vector<int> diff = ComputeDiff(current, parent);
		float range_min = additional_parameters["range_min"];
		float range_max = additional_parameters["range_max"];
		float mutation_size = additional_parameters["relative_mutation_size"] * (range_max - range_min);
		int recombinations = FastRand::PolinomialInt(recombination_rate, diff.size());

		for (int i = 0; i < recombinations; i++) {
			int random_index = FastRand::RandomInt(diff.size());
			int diff_index = diff[random_index];
			float added_value = FastRand::RandomFloat(-mutation_size, +mutation_size);
			current.genes_[diff_index] = std::min(std::max(current.genes_[i] + added_value, range_min), range_max);

			diff.erase(diff.begin() + random_index);
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

std::vector<GeneticAlgorithm<bool>> MakeBinaryGas(std::function<float(Chromosome<bool>&, std::map<std::string, float>& additional_parameters)> fitness, bool optimised, int chromosome_length) {

	std::vector<GeneticAlgorithm<bool>> gas = std::vector<GeneticAlgorithm<bool>>();

	std::vector<int> population_sizes = { 100, 200, 300 };
	std::vector<float> mutation_probabilities = { 0.05f, 0.1f, 0.15f };
	std::vector<float> crossover_probabilities = { 0.6f, 0.65f, 0.7f, 0.75f, 0.8f };
	std::vector<float> elitism_rates = { 0.05f, 0.1f, 0.15f };

	std::vector<int> tournament_sizes = { 2, 3, 4, 5 };
	std::vector<float> gene_mutation_rates = { 0.05f, 0.1f, 0.15f };
	std::vector<float> recombination_rates;
	if (optimised) {
		recombination_rates = { 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f };
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

std::vector<GeneticAlgorithm<float>> MakeRealValuedGas(std::function<float(Chromosome<float>&, std::map<std::string, float>& additional_parameters)> fitness, bool optimised, int chromosome_length, float bound) {
	std::vector<GeneticAlgorithm<float>> gas = std::vector<GeneticAlgorithm<float>>();

	std::vector<int> population_sizes = { 100, 200, 300 };
	std::vector<float> mutation_probabilities = { 0.05f, 0.1f, 0.15f };
	std::vector<float> crossover_probabilities = { 0.6f, 0.65f, 0.7f, 0.75f, 0.8f };
	std::vector<float> elitism_rates = { 0.05f, 0.1f, 0.15f };

	std::vector<int> tournament_sizes = { 2, 3, 4, 5 };
	std::vector<float> gene_mutation_rates = { 0.05f, 0.1f, 0.15f };
	std::vector<float> relative_mutation_sizes = { 0.05f, 0.1f, 0.2f };

	std::vector<float> crossover_ratios = { 0, 0.25f };

	std::vector<float> recombination_rates;
	if (optimised) {
		recombination_rates = { 0.3f, 0.4f, 0.5f, 0.6f, 0.7f, 0.8f };
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

										ga.target_fitness_ = -0.001f;

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

										ga.max_fitness_evaluations_ = 100000;

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

template<typename T>
void RunCompleteTest(std::vector<GeneticAlgorithm<T>> standard_gas, std::vector<GeneticAlgorithm<T>> optimised_gas, int base_test_size, float test_size_increase_rate, float elimination_rate, int final_test_size, std::string directory) {
	std::cout << "Running Standard Exploration Test..." << std::endl;
	GeneticAlgorithm<T> best_standard = TestSuite<T>::SelectBestConfiguration(standard_gas, base_test_size, test_size_increase_rate, elimination_rate);

	std::cout << "Standard Exploration Test completed! Winner: " << std::endl;
	std::cout << best_standard.DumpParameters() << std::endl;
	std::cout << "Running Standard Main Test..." << std::endl;

	float best_standard_evaluations = TestSuite<T>::RunBaseTest(best_standard, final_test_size);

	std::cout << "Standard Main Test finished! Evaluations: " << std::to_string(best_standard_evaluations) << std::endl;

	std::cout << "Running Optimised Exploration Test..." << std::endl;
	GeneticAlgorithm<T> best_optimised = TestSuite<T>::SelectBestConfiguration(optimised_gas, base_test_size, test_size_increase_rate, elimination_rate);

	std::cout << "Optimised Exploration Test completed! Winner: " << std::endl;
	std::cout << best_optimised.DumpParameters() << std::endl;
	std::cout << "Running Optimised Main Test..." << std::endl;

	float best_optimised_evaluations = TestSuite<T>::RunBaseTest(best_optimised, final_test_size);

	std::cout << "Optimised Main Test finished! Evaluations: " << std::to_string(best_optimised_evaluations) << std::endl;

	directory.erase(directory.find_last_of('\\') + 1);
	std::string finish_time = FormatTime(std::time(nullptr), "%d-%m-%y %H-%M-%S");
	std::string file_path = directory + "GA Results " + finish_time + ".txt";

	std::ofstream result_file(file_path);

	result_file << "BEST STANDARD (" << best_standard_evaluations << "):\n";
	result_file << best_standard.DumpParameters() << "\n";
	result_file << "BEST OPTIMISED(" << best_optimised_evaluations << "):\n";
	result_file << best_optimised.DumpParameters() << "\n";

	result_file.flush();
	result_file.close();

	system(("notepad.exe " + file_path).c_str());
}

std::string DataPointToCSVLine(DataPoint datapoint, std::string separator) {
	return std::to_string(datapoint.average) + separator + std::to_string(datapoint.standard_deviation) + separator +
		std::to_string(datapoint.min) + separator +
		std::to_string(datapoint.q1) + separator +
		std::to_string(datapoint.median) + separator +
		std::to_string(datapoint.q3) + separator +
		std::to_string(datapoint.max);
}

void SaveDetailedTest(std::pair<std::vector<DataPoint>, std::vector<DataPoint>> result, int snapshot_period, std::string extra_info, std::string directory) {
	directory.erase(directory.find_last_of('\\') + 1);
	std::string finish_time = FormatTime(std::time(nullptr), "%d-%m-%y %H-%M-%S");
	std::string file_path = directory + "GA Detailed Results " + finish_time + ".txt";

	std::ofstream result_file(file_path);

	result_file << extra_info << "\n\n";
	result_file << "FITNESS STATS:\n";
	result_file << "Evaluations;Average;Standard Deviation;Min;Q1;Median(Q2);Q3;Max\n";

	int snapshot_evaluations_fitness = 0;
	for (DataPoint datapoint : result.first) {
		result_file << snapshot_evaluations_fitness << ";" << DataPointToCSVLine(datapoint, ";") << "\n";
		snapshot_evaluations_fitness += snapshot_period;
	}

	result_file << extra_info << "\n\n";
	result_file << "DIVERSITY STATS:\n";
	result_file << "Evaluations;Average;Standard Deviation;Min;Q1;Median(Q2);Q3;Max\n";

	int snapshot_evaluations_diversity = 0;
	for (DataPoint datapoint : result.second) {
		result_file << snapshot_evaluations_diversity << ";" << DataPointToCSVLine(datapoint, ";") << "\n";
		snapshot_evaluations_diversity += snapshot_period;
	}

	result_file.flush();
	result_file.close();

	system(("notepad.exe " + file_path).c_str());
}

template<typename T>
void RunCompleteDetailedTest(GeneticAlgorithm<T> ga, int test_size, int snapshot_period, std::string directory) {
	std::pair<std::vector<DataPoint>, std::vector<DataPoint>> result = TestSuite<T>::RunDetailedTest(ga, test_size, snapshot_period);
	SaveDetailedTest(result, snapshot_period,
		ga.DumpParameters() + "\n\nTest size: " + std::to_string(test_size) + "\nSnapshot Period: " + std::to_string(snapshot_period) + "\nEvaluations: " + std::to_string(ga.max_fitness_evaluations_) + "\n",
		directory);
}

GeneticAlgorithm<bool> MakeOneMax(int max_evaluations, bool optimised) {
	GeneticAlgorithm<bool> ga = GeneticAlgorithm<bool>(BinaryInitialization, TournamentSelection<bool>, TwoPointCrossover<bool>, BitFlipMutation, BinaryRecombination, OneMaxFitness);
	ga.max_fitness_evaluations_ = max_evaluations;
	ga.chromosome_length_ = 80;
	ga.population_size_ = optimised ? 100 : 200;
	ga.mutation_probability_ = 0.15f;
	ga.crossover_probability_ = 0.8f;
	ga.elitism_rate_ = 0.15f;
	ga.recombination_rate_ = optimised ? 0.5f : 0;
	ga.additional_parameters_["gene_mutation_rate"] = 0.05f;
	ga.additional_parameters_["tournament_size"] = optimised ? 2 : 4;

	return ga;
}

GeneticAlgorithm<float> MakeSphere(int max_evaluations, bool optimised) {
	GeneticAlgorithm<float> ga = GeneticAlgorithm<float>(
		UniformInitialization,
		TournamentSelection<float>,
		IntermediateCrossover,
		RealValuedMutation,
		RealValuedRecombination,
		SphereFitness);

	ga.max_fitness_evaluations_ = max_evaluations;
	ga.chromosome_length_ = 20;

	ga.population_size_ = 100;
	ga.mutation_probability_ = optimised ? 0.05f : 0.15f;
	ga.crossover_probability_ = 0.6f;
	ga.elitism_rate_ = optimised ? 0.05f : 0.15f;
	ga.recombination_rate_ = optimised ? 0.3f : 0;

	ga.additional_parameters_["range_min"] = -5.12f;
	ga.additional_parameters_["range_max"] = 5.12f;

	ga.additional_parameters_["tournament_size"] = optimised ? 3 : 5;
	ga.additional_parameters_["gene_mutation_rate"] = 0.05f;
	ga.additional_parameters_["relative_mutation_size"] = 0.05f;
	ga.additional_parameters_["crossover_ratio"] = 0.25f;

	return ga;
}

int main(int argc, char **argv)
{
	std::string directory = argv[0];

	/*std::vector<GeneticAlgorithm<float>> standard_gas = MakeRealValuedGas(RosenbrockFitness, false, 5, 2.048f);
	std::vector<GeneticAlgorithm<float>> optimised_gas = MakeRealValuedGas(RosenbrockFitness, true, 5, 2.048f);

	int base_test_size = 20;
	float test_size_increase_rate = 2;
	float elimination_rate = 0.85f;

	int final_test_size = 100000;*/

	//RunCompleteTest(standard_gas, optimised_gas, base_test_size, test_size_increase_rate, elimination_rate, final_test_size, directory);
	RunCompleteDetailedTest(MakeOneMax(10000, false), 100000, 200, directory);
	system("PAUSE");
	return 0;
}