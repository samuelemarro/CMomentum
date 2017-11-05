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
#include "json.hpp"

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

								GeneticAlgorithmParameters gap;

								gap.target_fitness_ = 0;

								gap.chromosome_length_ = chromosome_length;

								gap.population_size_ = population_size;
								gap.mutation_probability_ = mutation_probability;
								gap.crossover_probability_ = crossover_probability;
								gap.elitism_rate_ = elitism_rate;

								gap.additional_parameters_["tournament_size"] = tournament_size;
								gap.additional_parameters_["gene_mutation_rate"] = gene_mutation_rate;

								gap.recombination_rate_ = recombination_rate;

								gap.max_fitness_evaluations_ = max_fitness_evaluations;
								gap.max_stagnation_ = max_stagnation;

								GeneticAlgorithm<bool> ga = GeneticAlgorithm<bool>(
									gap,
									BinaryInitialization,
									TournamentSelection<bool>,
									TwoPointCrossover<bool>,
									BitFlipMutation,
									BinaryRecombination,
									fitness);

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

										GeneticAlgorithmParameters gap = GeneticAlgorithmParameters();

										gap.target_fitness_ = target;

										gap.chromosome_length_ = chromosome_length;

										gap.population_size_ = population_size;
										gap.mutation_probability_ = mutation_probability;
										gap.crossover_probability_ = crossover_probability;
										gap.elitism_rate_ = elitism_rate;

										gap.additional_parameters_["range_min"] = -bound;
										gap.additional_parameters_["range_max"] = bound;

										gap.additional_parameters_["tournament_size"] = tournament_size;
										gap.additional_parameters_["gene_mutation_rate"] = gene_mutation_rate;
										gap.additional_parameters_["relative_mutation_size"] = relative_mutation_size;
										gap.additional_parameters_["crossover_ratio"] = crossover_ratio;

										gap.recombination_rate_ = recombination_rate;

										gap.max_fitness_evaluations_ = max_fitness_evaluations;
										gap.max_stagnation_ = max_stagnation;

										GeneticAlgorithm<float> ga = GeneticAlgorithm<float>(
											gap,
											UniformInitialization,
											TournamentSelection<float>,
											IntermediateCrossover,
											RealValuedMutation,
											RealValuedRecombination,
											fitness);

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
	std::string str = ss.str();
	str.erase(str.find_last_not_of('0') + 1, std::string::npos);
	if (str[str.size() - 1] == '.') {
		str.pop_back();
	}
	return str;
}

std::string ScientificNotation(float x, const int decimal_digits) {
	std::stringstream ss;
	ss.precision(decimal_digits);
	ss << std::scientific << x;
	std::string str = ss.str();

	std::string mantissa = str.substr(0, str.find('e'));
	std::string exponent = str.substr(str.find('e') + 1);

	mantissa.erase(mantissa.find_last_not_of('0') + 1, std::string::npos);
	if (mantissa[mantissa.size() - 1] == '.') {
		mantissa.pop_back();
	}

	exponent.erase(exponent.find('0'), exponent.find_last_not_of('0') - 1);

	//str.erase(str.find_last_not_of('0'), str.find_first_of('e'));
	return mantissa + "e" + exponent;
}

template<typename T>
void RunCompleteTest(std::string name, GeneticAlgorithm<T> ga, int test_size, int successful_executions_section_size, int best_fitness_values_section_size, std::string directory) {

	std::cout << "Running Success Rate Test...\n";
	TestSuite::TestResults result = TestSuite::CompleteTest(ga, test_size, successful_executions_section_size, best_fitness_values_section_size);

	directory.erase(directory.find_last_of('\\') + 1);

	std::string file_name = name + " " + std::to_string(ga.parameters_.chromosome_length_);
	if (ga.parameters_.additional_parameters_.find("range_max") != ga.parameters_.additional_parameters_.end()) {
		file_name += " " + FloatFormat(ga.parameters_.additional_parameters_["range_max"], 2);
	}

	if (ga.parameters_.target_fitness_ != 0) {
		file_name += " " + ScientificNotation(ga.parameters_.target_fitness_, 2);
	}

	file_name += std::string(" ") + (ga.parameters_.recombination_rate_ == 0 ? "Standard" : "Optimised");

	std::string finish_time = FormatTime(std::time(nullptr), "%d-%m-%y %H-%M-%S");

	std::string file_path = directory + file_name + " - GA Results " + finish_time + ".txt";

	std::ofstream result_file(file_path);

	result_file << "EXECUTION PARAMETERS:\n";
	result_file << "Name: " << name << "\n";
	result_file << "Test Size: " << std::to_string(test_size) << "\n";
	result_file << "Success Rate Section Size: " << std::to_string(successful_executions_section_size) << "\n";
	result_file << "Best Fitness Values Section Size: " << std::to_string(best_fitness_values_section_size) << "\n\n";

	result_file << "GA PARAMETERS:\n";
	result_file << ga.parameters_.ToString() << "\n\n";

	result_file << "Success Rate: " << std::to_string(std::accumulate(result.successful_executions_distribution_.begin(), result.successful_executions_distribution_.end(), 0.0f)) << "\n\n";

	result_file << "OVERALL BEST FITNESS STATISTICS:\n";
	result_file << result.overall_best_fitness_values_stats_.ToString() << "\n\n";

	result_file << "FINAL EVALUATION STATISTICS:\n";
	result_file << result.final_evaluations_stats_.ToString() << "\n\n";


	result_file << "AVERAGE BEST FITNESS PLOT:\n";

	int best_fitness_values_section_counter = best_fitness_values_section_size;
	for (int i = 0; i < result.average_best_fitness_values_plot_.size(); i++) {
		result_file << std::to_string(best_fitness_values_section_counter) << ";" << std::to_string(result.average_best_fitness_values_plot_[i]) << "\n";
		best_fitness_values_section_counter += best_fitness_values_section_size;
	}

	result_file << "\n\nSUCCESSFUL EXECUTIONS DISTRIBUTION:\n";

	int successful_executions_section_counter = successful_executions_section_size;
	for (int i = 0; i < result.successful_executions_distribution_.size(); i++) {
		result_file << std::to_string(successful_executions_section_counter) << ";" << std::to_string(result.successful_executions_distribution_[i]) << "\n";
		successful_executions_section_counter += best_fitness_values_section_size;
	}

	result_file.flush();
	result_file.close();
}


GeneticAlgorithm<float> MakeRastrigin5Standard(bool original) {

	GeneticAlgorithmParameters gap = GeneticAlgorithmParameters();

	gap.chromosome_length_ = 5;

	gap.population_size_ = original ? 100 : 200;
	gap.mutation_probability_ = 0.15f;
	gap.crossover_probability_ = 0.65f;
	gap.elitism_rate_ = original ? 0.05f : 0.1f;
	gap.recombination_rate_ = 0;

	gap.additional_parameters_["range_min"] = -5.12f;
	gap.additional_parameters_["range_max"] = 5.12f;

	gap.additional_parameters_["tournament_size"] = 2;
	gap.additional_parameters_["gene_mutation_rate"] = 0.15f;
	gap.additional_parameters_["relative_mutation_size"] = original ? 0.1f : 0.2f;
	gap.additional_parameters_["crossover_ratio"] = 0.25f;

	gap.target_fitness_ = -1e-5;
	gap.max_fitness_evaluations_ = 100000;

	GeneticAlgorithm<float> ga = GeneticAlgorithm<float>(
		gap,
		UniformInitialization,
		TournamentSelection<float>,
		IntermediateCrossover,
		RealValuedMutation,
		RealValuedRecombination,
		RastriginFitness);

	return ga;
}

void SaveParameters(std::string path, GeneticAlgorithmParameters gap, int indentation) {
	std::ofstream file(path);
	nlohmann::json j = gap;
	file << j.dump(indentation);
	file.flush();
	file.close();
}

GeneticAlgorithmParameters LoadParameters(std::string path) {
	std::ifstream file(path);
	nlohmann::json j = nlohmann::json();
	j << file;
	file.close();
	return j;
}

void ClearScreen() {
	for (int i = 0; i < 10; i++) {
		std::cout << "\n\n\n\n\n\n\n\n\n\n";
	}
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

	std::vector<GeneticAlgorithm<float>> gas = MakeRealValuedGas(RastriginFitness, false, 5, 5.12f, -1e-5f, 100000, -1);
	//to_json(j, gas[0]);

	std::string path = "C:\\Users\\Samuele\\Documents\\visual studio 2017\\Projects\\CMomentumv2\\x64\\Release\\Rastrigin5.gap";

	float target_success_rate = 0.95f;

	//PROBLEMA: Se il succrate � 0 viene -inf (quindi estremamente favorevole)
	//TODO: Controllare >/>=
	while (true) {
		ClearScreen();
		GeneticAlgorithmParameters gap = LoadParameters(path);
		GeneticAlgorithm<float> ga = GeneticAlgorithm<float>(gap, UniformInitialization, TournamentSelection<float>, IntermediateCrossover, RealValuedMutation, RealValuedRecombination, RastriginFitness);

		//RunCompleteTest("Rastrigin", ga, final_test_size, 100, 100, directory);
		TestSuite::TestResults result = TestSuite::CompleteTest(ga, 20, 100, 200);

		float success_rate = std::accumulate(result.successful_executions_distribution_.begin(), result.successful_executions_distribution_.end(), 0.0f);
		float average_evaluation = result.final_evaluations_stats_.average;

		float N = std::log(1 - target_success_rate) / std::log(1 - success_rate);
		float actual_average_evaluation = average_evaluation * N;

		std::cout << "Average Best Fitness: " << FloatFormat(result.overall_best_fitness_values_stats_.average, 8) << "\n";
		std::cout << "Success Rate: " << FloatFormat(success_rate, 5) << "\n";
		std::cout << "Average Evaluations: " << FloatFormat(average_evaluation, 2) << "\n";
		std::cout << "N: " << FloatFormat(N, 2) << "\n";
		std::cout << "Actual Average Evaluations: " << FloatFormat(actual_average_evaluation, 2) << "\n";
		system("PAUSE");
	}

	return 0;
}