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

static class MathUtility {
	template<typename T>
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
			median = (vector[(vector.size() - 1) / 2] + vector[(vector.size() - 1) / 2 + 1]) / 2;
		}
		else {
			median = vector[(vector.size() - 1) / 2];
		}

		return median;
	}
	static float StandardDeviation(std::vector<float> data) {
		return StandardDeviation(data, std::accumulate(data.begin(), data.end(), 0) / static_cast<float>(data.size()));
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
	float min;
	float max;
	float average;
	float standard_deviation;
	float median;
	float q1;
	float q3;

	DataPoint() {

	}
	DataPoint(std::vector<float> data) : DataPoint(data, false) {

	}
	DataPoint(std::vector<float> data, bool is_sorted) {
		if (data.size() < 2) {
			throw std::invalid_argument("The size of data must be at least 2.");
		}

		min = *std::min_element(data.begin(), data.end());
		max = *std::max_element(data.begin(), data.end());
		if (!is_sorted) {
			std::sort(data.begin(), data.end());
		}
		average = std::accumulate(data.begin(), data.end(), 0) / static_cast<float>(data.size());
		standard_deviation = MathUtility::StandardDeviation(data, average);
		median = MathUtility::Median(data, true);

		//Remove the middle element to split correctly data in two parts
		if (data.size() % 2 != 0 && data.size() > 1) {
			int middle_position = (data.size() + 1) / 2;
			data.erase(data.begin() + middle_position);
		}
		int half_size = data.size() / 2;

		q1 = MathUtility::Median(std::vector<float>(data.begin(), data.begin() + half_size), true);
		q3 = MathUtility::Median(std::vector<float>(data.begin() + half_size, data.end()), true);
	}
};

enum TrackingType {
	Fitness,
	Diversity,
	Both
};

template<typename T>
class TestSuite {
	template<typename T>
	struct PartialTestResult {
		friend class TestSuite<T>;
		std::vector<int> evaluations_;
		GeneticAlgorithm<T> ga_;
		float median_ = INT32_MIN;

		TestSuite::PartialTestResult<T>(GeneticAlgorithm<T> ga, int evaluations_size) : ga_(ga) {
			evaluations_.reserve(evaluations_size);
		}
	};

public:
	static std::vector<TestSuite::PartialTestResult<T>> RunParallelTests(std::vector<GeneticAlgorithm<T>> gas, int test_size) {

		std::vector<TestSuite::PartialTestResult<T>> results = std::vector<TestSuite::PartialTestResult<T>>();

		for (int i = 0; i < gas.size(); i++) {
			results.push_back(PartialTestResult<T>(gas[i], test_size));
		}

		int gas_size = gas.size();

		int total_tests = gas.size() * test_size;
		int executed_tests = 0;
		std::random_device seed_generator;
#pragma omp parallel
		{
#pragma omp for
			for (int i = 0; i < gas_size * test_size; i++) {
				FastRand::Seed(seed_generator());
				int ga_index = i % gas_size;
				GeneticAlgorithm<T> ga = gas[ga_index];
				float evaluation = ga.RunAlgorithm(false, false, 0);
#pragma omp critical
				{
					results[ga_index].evaluations_.push_back(evaluation);
				}

				executed_tests++;
				if (executed_tests % std::max(1, total_tests / 100) == 0) {
					EraseWriteLine("Progress: " + std::to_string(executed_tests * 100 / total_tests) + "%");
				}
			}
		}

		for (int i = 0; i < results.size(); i++)
		{
			results[i].median_ = MathUtility::Median(results[i].evaluations_);
		}

		std::cout << "\n";

		return results;
	}
	static std::vector<int> EvaluationsTest(GeneticAlgorithm<T> base_ga, int test_size) {
		std::vector<int> evaluations = std::vector<int>();
		int executed_tests = 0;
		std::random_device seed_generator;
#pragma omp parallel
		{
#pragma omp for
			for (int i = 1; i <= test_size; i++) {
				FastRand::Seed(seed_generator());
				GeneticAlgorithm<T> ga = base_ga;
				float evaluation = ga.RunAlgorithm(false, false, 0);
#pragma omp critical
				{
					evaluations.push_back(evaluation);
				}

				executed_tests++;
				if (executed_tests % std::max(1, test_size / 100) == 0) {
					EraseWriteLine("Progress: " + std::to_string(executed_tests * 100 / test_size) + "%");
				}
			}
		}

		std::cout << "\n";
		return evaluations;
	}
	static float RunBaseTest(GeneticAlgorithm<T> base_ga, int test_size) {
		std::vector<int> evaluations = EvaluationsTest(base_ga, test_size);
		return MathUtility::Median(evaluations);
	}

	static std::vector<float> SuccessRatesTest(GeneticAlgorithm<T> base_ga, int test_size, int section_size) {
		std::vector<int> evaluations = EvaluationsTest(base_ga, test_size);
		std::vector<int> successful_tests = std::vector<int>();
		std::sort(evaluations.begin(), evaluations.end());
		int current_section = 0;
		for(int evaluation : evaluations) {
			int floored = evaluation - (evaluation % section_size);
			int index = floored / section_size;
			if (successful_tests.size() < index + 1) {
				successful_tests.resize(index + 1);
			}
			successful_tests[index]++;
		}

		std::vector<float> success_rates = std::vector<float>();

		//Compute the success rate
		success_rates.push_back(static_cast<float>(successful_tests[0]) / evaluations.size());
		for (int i = 1; i < successful_tests.size(); i++) {
			success_rates.push_back(static_cast<float>(successful_tests[i]) / evaluations.size());
		}

		return success_rates;
	}

	static std::pair<std::vector<DataPoint>, std::vector<DataPoint>> RunPopulationTest(GeneticAlgorithm<T> base_ga, int test_size, TrackingType tracking_type, int snapshot_period) {

		if (base_ga.target_fitness_ != FLT_MAX || base_ga.max_generation_ != -1) {
			std::cout << "WARNING: This type of test is designed for GAs that are executed until they reach a certain evaluation number." << std::endl;
			std::cout << "Using other types of termination might lead to unexpected results, for example the best fitness worsening after reaching the termination." << std::endl;
			std::cout << "This is caused by the fact that GAs that have reached the termination stop outputting further results." << std::endl;
		}

		bool track_fitness = tracking_type == TrackingType::Fitness || tracking_type == TrackingType::Both;
		bool track_diversity = tracking_type == TrackingType::Diversity || tracking_type == TrackingType::Both;

		//Fitness and diversity tracking. Each inner vector corresponds to a generation
		std::vector<std::vector<float>> final_fitness_values = std::vector<std::vector<float>>();
		std::vector<std::vector<float>> final_diversity_values = std::vector<std::vector<float>>();

		int executed_tests = 0;
		std::random_device seed_generator;
#pragma omp parallel
		{
#pragma omp for
			for (int i = 1; i <= test_size; i++) {
				FastRand::Seed(seed_generator());
				GeneticAlgorithm<T> ga = base_ga;
				ga.RunAlgorithm(track_fitness, track_diversity, snapshot_period);

#pragma omp critical
				{
					//If final_fitness_values has 500 generations and tracked_fitness_values has 1000, make room for more
					if (track_fitness && final_fitness_values.size() < ga.tracked_fitness_values.size()) {
						final_fitness_values.resize(ga.tracked_fitness_values.size(), std::vector<float>());
					}
					//Same goes for diversity
					if (track_diversity && final_diversity_values.size() < ga.tracked_diversity_values.size()) {
						final_diversity_values.resize(ga.tracked_diversity_values.size(), std::vector<float>());
					}

					//Add the fitness values of each generation to its corresponding slot
					if (track_fitness) {
						for (int i = 0; i < ga.tracked_fitness_values.size(); i++) {
							final_fitness_values[i].insert(final_fitness_values[i].end(), ga.tracked_fitness_values[i].begin(), ga.tracked_fitness_values[i].end());
						}
					}
					//Add the diversity value of each generation to its corresponding slot
					if(track_diversity){
						for (int i = 0; i < ga.tracked_diversity_values.size(); i++) {
							final_diversity_values[i].push_back(ga.tracked_diversity_values[i]);
						}
					}
				}

				executed_tests++;
				if (executed_tests % std::max(1, test_size / 100) == 0) {
					EraseWriteLine("Progress: " + std::to_string(executed_tests * 100 / test_size) + "%");
				}
			}
		}

		std::vector<DataPoint> fitness_datapoints = std::vector<DataPoint>();
		std::vector<DataPoint> diversity_datapoints = std::vector<DataPoint>();
		if (track_fitness) {
			for (int i = 0; i < final_fitness_values.size(); i++) {
				fitness_datapoints.push_back(DataPoint(final_fitness_values[i]));
			}
		}

		if (track_diversity) {
			for (int i = 0; i < final_diversity_values.size(); i++) {
				diversity_datapoints.push_back(DataPoint(final_diversity_values[i]));
			}
		}

		return std::make_pair(fitness_datapoints, diversity_datapoints);
	}
	static GeneticAlgorithm<T> SelectBestConfiguration(std::vector<GeneticAlgorithm<T>> gas, int base_test_size, float test_size_increase_rate, float tournament_elimination) {

		std::vector<TestSuite::PartialTestResult<T>> tournament = std::vector<TestSuite::PartialTestResult<T>>();
		int tournament_size = gas.size();
		int test_size = base_test_size;
		int tournament_number = 1;

		//tournament.resize requires an object to use in case the vector's size is increased (which is not going to happen)
		TestSuite::PartialTestResult<T> empty_test = PartialTestResult<T>(PartialTestResult<T>(gas[0], 0));

		while (tournament_size > 1) {
			std::cout << "Running tournament n." << std::to_string(tournament_number) << ": " << std::to_string(test_size) << " tests for " << std::to_string(tournament_size) << " configurations" << std::endl;

			tournament = RunParallelTests(gas, test_size);

			//Sort in ascending order by executed evaluations
			std::sort(tournament.begin(), tournament.end(),
				[](const auto& lhs, const auto& rhs) {
				return lhs.median_ < rhs.median_;
			});

			//If the tournament size doesn't change (due to rounding), the tournament size is reduced by 1. The tournament size must also be at least 1.
			tournament_size = std::max(1, std::min(static_cast<int>(static_cast<float>(tournament_size) * (1 - tournament_elimination)), tournament_size - 1));
			test_size *= (1 + test_size_increase_rate);

			//The first tournament will be resized to the same size as before
			tournament.resize(tournament_size, empty_test);

			gas.clear();
			gas.reserve(tournament_size);
			for (int i = 0; i < tournament_size; i++) {
				gas.push_back(tournament[i].ga_);
			}

			tournament_number++;
		}
		return tournament[0].ga_;
	}
	static void EraseWriteLine(std::string text) {
		std::cout << text;
		std::cout << std::string(text.length(), '\b');
	}
};