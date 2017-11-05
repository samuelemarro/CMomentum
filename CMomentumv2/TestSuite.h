#pragma once
#include "stdafx.h"
#include <vector>
#include <thread>
#include <process.h>
#include <random>
#include <stdlib.h>
#include <math.h>
#include <numeric>
#include <type_traits>
#include <omp.h>
#include <limits>

static class MathUtility {
	friend class TestSuite;
	friend class DataPoint;
private:
	template<typename T,
		typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
		static float Median(std::vector<T> vector) {
		return Median(vector, false);
	}
	template<typename T,
		typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
		static float Median(std::vector<T> vector, bool is_sorted) {

		float median = 0;

		if (!is_sorted) {
			std::sort(vector.begin(), vector.end());
		}

		if (vector.size() % 2 == 0) {
			median = (float)(vector[vector.size() / 2 - 1] + vector[vector.size() / 2]) / 2;
		}
		else {
			median = vector[(vector.size() - 1) / 2];
		}

		return median;
	}

	template<typename T,
		typename = typename std::enable_if<std::is_arithmetic<T>::value, T>::type>
		static float Average(std::vector<T> vector) {
		float sum = std::accumulate(vector.begin(), vector.end(), 0.0f);
		return sum / static_cast<float>(vector.size());
	}

	static float StandardDeviation(std::vector<float> data) {
		return StandardDeviation(data, Average(data));
	}
	static float StandardDeviation(std::vector<float> data, float average) {
		float sum = 0;
		for (float value : data) {
			sum += (value - average) * (value - average);
		}
		return std::sqrt(sum / data.size());
	}
};

struct DataPoint {
public:
	float min = std::numeric_limits<float>::has_signaling_NaN ? std::numeric_limits<float>::signaling_NaN() : FLT_MAX;
	float max = std::numeric_limits<float>::has_signaling_NaN ? std::numeric_limits<float>::signaling_NaN() : FLT_MAX;
	float average = std::numeric_limits<float>::has_signaling_NaN ? std::numeric_limits<float>::signaling_NaN() : FLT_MAX;
	float standard_deviation = std::numeric_limits<float>::has_signaling_NaN ? std::numeric_limits<float>::signaling_NaN() : FLT_MAX;
	float median = std::numeric_limits<float>::has_signaling_NaN ? std::numeric_limits<float>::signaling_NaN() : FLT_MAX;
	float q1 = std::numeric_limits<float>::has_signaling_NaN ? std::numeric_limits<float>::signaling_NaN() : FLT_MAX;
	float q3 = std::numeric_limits<float>::has_signaling_NaN ? std::numeric_limits<float>::signaling_NaN() : FLT_MAX;

	DataPoint() {

	}
	DataPoint(std::vector<float> data) : DataPoint(data, false) {

	}
	DataPoint(std::vector<float> data, bool is_sorted) {
		if (data.size() < 1) {
			throw std::invalid_argument("The size of data must be at least 1.");
		}

		min = *std::min_element(data.begin(), data.end());
		max = *std::max_element(data.begin(), data.end());
		if (!is_sorted) {
			std::sort(data.begin(), data.end());
		}
		average = MathUtility::Average(data);
		standard_deviation = MathUtility::StandardDeviation(data, average);
		median = MathUtility::Median(data, true);

		if (data.size() == 1) {
			q1 = median;
			q3 = median;
		}
		else {
			//Remove the middle element to split correctly data in two parts
			if (data.size() % 2 != 0) {
				int middle_position = (data.size() + 1) / 2;
				data.erase(data.begin() + middle_position);
			}
			int half_size = data.size() / 2;

			std::vector<float> first_half = std::vector<float>(data.begin(), data.begin() + half_size);
			std::vector<float> second_half = std::vector<float>(data.begin() + half_size, data.end());

			q1 = MathUtility::Median(first_half, true);
			q3 = MathUtility::Median(second_half, true);
		}
	}

	std::string ToString() {
		std::string result = "";
		result += "Min: " + std::to_string(min) + "\n";
		result += "Max: " + std::to_string(max) + "\n";
		result += "Average: " + std::to_string(average) + "\n";
		result += "Standard Deviation: " + std::to_string(standard_deviation) + "\n";
		result += "Q1: " + std::to_string(q1) + "\n";
		result += "Median (Q2): " + std::to_string(median) + "\n";
		result += "Q3: " + std::to_string(q3) + "\n";

		return result;
	}

};

class TestSuite {
public:

	template<typename T>
	static std::vector<std::vector<std::pair<int, float>>> ParallelTest(GeneticAlgorithm<T> base_ga, int test_size) {
		std::vector<std::vector<std::pair<int, float>>> results = std::vector<std::vector<std::pair<int, float>>>();
		int executed_tests = 0;
		std::random_device seed_generator;
#pragma omp parallel
		{
#pragma omp for
			for (int i = 1; i <= test_size; i++) {
				FastRand::Seed(seed_generator());
				GeneticAlgorithm<T> ga = base_ga;
				std::vector<std::pair<int, float>> pair = ga.RunAlgorithm();
#pragma omp critical
				{
					results.push_back(pair);
				}

				executed_tests++;
				if (executed_tests % std::max(1, test_size / 100) == 0) {
					EraseWriteLine("Progress: " + std::to_string(executed_tests * 100 / test_size) + "%");
				}
			}
		}

		std::cout << "\n";
		return results;
	}

	struct TestResults {
		DataPoint final_evaluations_stats_;
		DataPoint overall_best_fitness_values_stats_;
		std::vector<float> successful_executions_distribution_;
		std::vector<float> average_best_fitness_values_plot_;

		TestResults(DataPoint final_evaluations_stats,
			DataPoint overall_best_fitness_values_stats,
			std::vector<float> successful_executions_distribution,
			std::vector<float> average_best_fitness_values_plot)
			: final_evaluations_stats_(final_evaluations_stats),
			overall_best_fitness_values_stats_(overall_best_fitness_values_stats),
			successful_executions_distribution_(successful_executions_distribution),
			average_best_fitness_values_plot_(average_best_fitness_values_plot)
		{

		}
	};

	template<typename T>
	static TestResults CompleteTest(GeneticAlgorithm<T> base_ga, int test_size, int successful_executions_section_size, int best_fitness_values_section_size) {
		
		if (best_fitness_values_section_size < base_ga.parameters_.population_size_) {
			throw std::invalid_argument("\"best_fitness_values_section_size\" must be bigger than or equal to the population size");
		}
		
		std::vector<std::vector<std::pair<int, float>>> results = ParallelTest(base_ga, test_size);

		//======Find the final evaluations for each execution======

		std::vector<float> final_evaluations = std::vector<float>();
		for (std::vector<std::pair<int, float>> execution : results) {
			final_evaluations.push_back(execution[execution.size() - 1].first);
		}

		DataPoint final_evaluations_stats = final_evaluations.size() > 0 ? DataPoint(final_evaluations) : DataPoint();

		//======Find the best fitness for each execution======

		std::vector<float> overall_best_fitness_values = std::vector<float>();
		for (std::vector<std::pair<int, float>> execution : results) {
			float overall_best_fitness = -FLT_MAX;
			//Each pair stores the best fitness (stored_data.second) after the specific amount of evaluations (stored_data.first)
			for (std::pair<int, float> stored_data : execution) {
				if (stored_data.second > overall_best_fitness) {
					overall_best_fitness = stored_data.second;
				}
			}
			overall_best_fitness_values.push_back(overall_best_fitness);
		}

		DataPoint overall_best_fitness_values_stats = overall_best_fitness_values.size() > 0 ? DataPoint(overall_best_fitness_values) : DataPoint();

		//======Find the distribution of the successful executions======
		std::vector<float> successful_executions_distribution = std::vector<float>();
		for (int i = 0; i < results.size(); i++) {
			std::vector<std::pair<int, float>> execution = results[i];

			float overall_best_fitness = overall_best_fitness_values[i];
			int total_evaluations = execution[execution.size() - 1].first;

			//If the execution was successful
			if (overall_best_fitness >= base_ga.parameters_.target_fitness_) {
				int final_fitness_evaluations = execution[execution.size() - 1].first;

				//Find the section
				int floored = total_evaluations - (total_evaluations % successful_executions_section_size);
				int section_index = floored / successful_executions_section_size;

				//Make sure that there are enough sections
				if (successful_executions_distribution.size() < section_index + 1) {
					successful_executions_distribution.resize(section_index + 1);
				}

				//Store the success in the corresponding section
				successful_executions_distribution[section_index] += 1 / (float)results.size();
			}
		}

		//======Plot the average best fitness values by evaluation number======
		std::vector<std::vector<float>> best_fitness_values = std::vector<std::vector<float>>();

		//Find the maximum number of sections
		int max_sections = 1 + (base_ga.parameters_.max_fitness_evaluations_ - (base_ga.parameters_.max_fitness_evaluations_ % best_fitness_values_section_size)) / best_fitness_values_section_size;

		//Prepare the grid
		for (int i = 0; i < max_sections; i++) {
			std::vector<float> section_fitness_values = std::vector<float>();
			
			for (int j = 0; j < results.size(); j++) {
				section_fitness_values.push_back(0);
			}

			best_fitness_values.push_back(section_fitness_values);
		}

		//Store the fitness values in the grid
		for (int i = 0; i < results.size(); i++) {
			for (int j = 0; j < max_sections; j++) {

				int min = j * best_fitness_values_section_size;
				int max = (j + 1) * best_fitness_values_section_size;
				//Find all the pairs for the current execution where the evaluation number is within the boundaries.
				//Store the best fitness for each

				std::vector<float> best_fitness_values_found = std::vector<float>();

				for (std::pair<int, float> pair : results[i]) {
					if (pair.first > min && pair.first <= max) {
						best_fitness_values_found.push_back(pair.second);
					}
					if (pair.first > max) {
						break;
					}
				}

				if (best_fitness_values_found.size() == 0) {
					//If it couldn't find values, store the previous value
					//Q: What if we're considering the first section?
					//A: That's impossible: The GA stores the initial best fitness in the first section, so at least one value will be always found

					best_fitness_values[j][i] = best_fitness_values[j - 1][i];
				}
				else if (best_fitness_values_found.size() == 1) {
					//If we could only find one value, we store it
					best_fitness_values[j][i] = best_fitness_values_found[0];
				}
				else {
					//If we found multiple values, store their average
					best_fitness_values[j][i] = MathUtility::Average(best_fitness_values_found);
				}
			}
		}

		std::vector<float> average_best_fitness_values_plot = std::vector<float>();

		//Find the average for each section
		for (std::vector<float> data : best_fitness_values) {
			average_best_fitness_values_plot.push_back(MathUtility::Average(data));
		}

		return TestResults(final_evaluations_stats, overall_best_fitness_values_stats, successful_executions_distribution, average_best_fitness_values_plot);
	}

	static void EraseWriteLine(std::string text) {
		std::cout << text;
		std::cout << std::string(text.length(), '\b');
	}
};