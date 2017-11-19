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

#define _USE_MATH_DEFINES
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
	float mutation_size = additional_parameters["mutation_size"];

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

void RealValuedRecombination(Chromosome<float>& current, float recombination_rate, std::map<std::string, float>& additional_parameters) {
	if (recombination_rate != 0) {
		float range_min = additional_parameters["range_min"];
		float range_max = additional_parameters["range_max"];

		for (int i = 0; i < current.genes_.size(); i++) {
			float y1 = (current.fitness_ - current.parent1_.fitness_) / current.fitness_;
			float y2 = (current.fitness_ - current.parent2_.fitness_) / current.fitness_;
			float x1 = (current.genes_[i] - current.parent1_.genes_[i]) / current.genes_[i];
			float x2 = (current.genes_[i] - current.parent2_.genes_[i]) / current.genes_[i];

			float m1 = x1 == 0 ? 0 : y1 / x1;
			float m2 = x2 == 0 ? 0 : y2 / x2;

			float m = m1 + m2;

			current.genes_[i] = std::min(std::max(current.genes_[i] + m * std::abs(current.parent1_.genes_[i] - current.parent2_.genes_[i]) * recombination_rate, range_min), range_max);
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
	float sum = 0;
	for each(float gene in chromosome.genes_) {
		sum += gene * gene;
	}
	return -sum;
}

float GriewankFitness(Chromosome<float>& chromosome, std::map<std::string, float>& additional_parameters) {
	float sum = 0;
	float product = 1;

	for (int i = 0; i < chromosome.genes_.size(); i++) {
		sum += chromosome.genes_[i] * chromosome.genes_[i];
		product *= std::cos(chromosome.genes_[i] / std::sqrt(i + 1));
	}

	float total = sum / 4000 - product + 1;
	return -total;
}

float RastriginFitness(Chromosome<float>& chromosome, std::map<std::string, float>& additional_parameters) {
	float sum = 0;
	for (int i = 0; i < chromosome.genes_.size(); i++) {
		sum += chromosome.genes_[i] * chromosome.genes_[i] - 10 * std::cos(2 * std::_Pi * chromosome.genes_[i]);
	}
	float total = 10 * chromosome.genes_.size() + sum;
	return -total;
}

float RosenbrockFitness(Chromosome<float>& chromosome, std::map<std::string, float>& additional_parameters) {
	float sum = 0;

	for (int i = 0; i < chromosome.genes_.size() - 1; i++) {

		float current_gene = chromosome.genes_[i];
		float next_gene = chromosome.genes_[i + 1];

		sum += 100 * (next_gene - current_gene * current_gene) * (next_gene - current_gene * current_gene) +
			(current_gene - 1) * (current_gene - 1);
	}
	return -sum;
}

float AckleyFitness(Chromosome<float>& chromosome, std::map<std::string, float>& additional_parameters) {
	//Using double instead of float due to some nasty approximation errors
	double square_sum = 0;
	double cos_sum = 0;
	for (int i = 0; i < chromosome.genes_.size(); i++) {
		square_sum += chromosome.genes_[i] * chromosome.genes_[i];
		cos_sum += std::cos(2 * M_PI * chromosome.genes_[i]);
	}

	double square_exp = -20 * std::exp(-0.2f * std::sqrt(square_sum / (double)chromosome.genes_.size()));
	double cos_exp = -std::exp(cos_sum / (double)chromosome.genes_.size());

	double total = square_exp + cos_exp + 20 + M_E;

	return -total;
}

float SchwefelFitness(Chromosome<float>& chromosome, std::map<std::string, float>& additional_parameters) {
	float sum = 0;
	for (int i = 0; i < chromosome.genes_.size(); i++) {
		sum += chromosome.genes_[i] * std::sin(std::sqrt(std::abs(chromosome.genes_[i])));
	}

	float total = 418.9829f * (float)chromosome.genes_.size() - sum;

	return -total;
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

float ActualEvaluationsNumber(int evaluations, float target_rate, float actual_rate) {
	float N = std::log(1 - target_rate) / std::log(1 - actual_rate);
	return evaluations * N;
}

int main(int argc, char **argv)
{
	std::string directory = argv[0];
	directory.erase(directory.find_last_of('\\') + 1);

	//TODO:
	//Automatizzare
	//Pausa

	//int base_test_size = 10;
	//float test_size_increase_rate = 1;
	//float elimination_rate = 0.9f;

	//int final_test_size = 100000;

	//Chromosome<float> c = Chromosome<float>();
	//c.genes_.push_back(420.9687f);
	//c.genes_.push_back(420.9687f);
	//float a = SchwefelFitness(c, std::map<std::string, float>());

	std::string name = "";
	std::string path = "";
	int test_size = -1;
	std::string function_name = "";
	bool use_adjusted = false;

	bool proceed = false;
	do {
		ClearScreen();
		std::cout << "Insert the name of the test file\n";
		std::cin >> name;
		path = directory + "Functions\\" + name + ".gap";
		proceed = (bool)std::ifstream(path);

		if (!proceed) {
			std::cout << "File not found. Do you want to create it? (y/n)\n";
			char proceed_answer;
			std::cin >> proceed_answer;
			if (proceed_answer == 'y') {
				SaveParameters(path, GeneticAlgorithmParameters(), 2);
			}
		}

	} while (!proceed);

	ClearScreen();

	std::cout << "Insert the number of tests\n";
	std::cin >> test_size;

	std::map<std::string, std::function<float(Chromosome<float>&, std::map<std::string, float>&)>> functions = {
		{"Rastrigin", RastriginFitness},
		{"Sphere", SphereFitness},
		{"Schwefel", SchwefelFitness},
		{"Rosenbrock", RosenbrockFitness},
		{"Griewank", GriewankFitness},
		{"Ackley", AckleyFitness}
	};

	do {
		ClearScreen();
		std::cout << "Insert the function name\n";
		std::cin >> function_name;
		proceed = functions.count(function_name) != 0;

	} while (!proceed);

	std::cout << "Do you want to use the adjusted evaluations value for ranking? (y/n)";
	char adjusted_answer;
	std::cin >> adjusted_answer;
	use_adjusted = adjusted_answer == 'y';

	float target_success_rate = 0.95f;

	char menu_answer;

	float best_evaluation = FLT_MAX;
	GeneticAlgorithmParameters best_gap = GeneticAlgorithmParameters();

	do {
		ClearScreen();
		GeneticAlgorithmParameters gap = LoadParameters(path);
		GeneticAlgorithm<float> ga = GeneticAlgorithm<float>(gap, UniformInitialization, TournamentSelection<float>, IntermediateCrossover, RealValuedMutation, RealValuedRecombination, functions[function_name]);

		//RunCompleteTest("Rastrigin", ga, 5, 100, 100, directory);
		TestSuite::TestResults result = TestSuite::CompleteTest(ga, test_size, 100, 200);

		float success_rate = std::accumulate(result.successful_executions_distribution_.begin(), result.successful_executions_distribution_.end(), 0.0f);

		float average_evaluation = result.final_evaluations_stats_.average;

		float N = std::log(1 - target_success_rate) / std::log(1 - success_rate);
		float adjusted_average_evaluation = average_evaluation * N;

		if (adjusted_average_evaluation < 0) {
			adjusted_average_evaluation = -adjusted_average_evaluation;
		}

		std::cout << "Average Overall Best Fitness: " << FloatFormat(result.overall_best_fitness_values_stats_.average, 8) << "\n";
		std::cout << "Success Rate: " << FloatFormat(success_rate, 5) << "\n";
		std::cout << "Average Evaluations: " << FloatFormat(average_evaluation, 2) << "\n";
		std::cout << "N: " << FloatFormat(N, 2) << "\n";
		std::cout << "Adjusted Average Evaluations: " << FloatFormat(adjusted_average_evaluation, 2) << "\n";

		float final_evaluation = use_adjusted ? adjusted_average_evaluation : average_evaluation;

		if (final_evaluation < best_evaluation) {
			std::cout << "FOUND NEW BEST!\n\n";
			best_evaluation = final_evaluation;
			best_gap = gap;
		}

		do {
			std::cout << "1. New test\n2. Show best\n3. Save best\n0. Quit\n";
			std::cin >> menu_answer;

			if (menu_answer == '2') {
				std::cout << best_gap.ToString() << "\n";
			}
			else if (menu_answer == '3') {
				SaveParameters(path, best_gap, 2);
			}

		} while (menu_answer != '1' && menu_answer != '0');

	} while (menu_answer != '0');

	return 0;
}