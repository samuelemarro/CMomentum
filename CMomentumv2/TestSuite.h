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
			median = (float)(vector[vector.size() / 2 - 1] + vector[vector.size() / 2]) / 2;
		}
		else {
			median = vector[(vector.size() - 1) / 2];
		}

		return median;
	}
	static float StandardDeviation(std::vector<float> data) {
		return StandardDeviation(data, std::accumulate(data.begin(), data.end(), 0.0f) / static_cast<float>(data.size()));
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
		average = std::accumulate(data.begin(), data.end(), 0.0f) / static_cast<float>(data.size());
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
				std::pair<bool, int> ga_result = ga.RunAlgorithm();
				int evaluation = ga_result.first ? ga_result.second : INT32_MAX;
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
	static std::vector<std::pair<bool, int>> EvaluationsTest(GeneticAlgorithm<T> base_ga, int test_size) {
		std::vector<std::pair<bool, int>> results = std::vector<std::pair<bool, int>>();
		int executed_tests = 0;
		std::random_device seed_generator;
#pragma omp parallel
		{
#pragma omp for
			for (int i = 1; i <= test_size; i++) {
				FastRand::Seed(seed_generator());
				GeneticAlgorithm<T> ga = base_ga;
				std::pair<bool, int> pair = ga.RunAlgorithm();
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

	//Returns a pair with 1. The success rate 2. The median number of evaluations.
	static std::pair<float, float> RunBaseTest(GeneticAlgorithm<T> base_ga, int test_size) {
		std::vector<std::pair<bool, int>> result = EvaluationsTest(base_ga, test_size);
		std::vector<int> evaluations = std::vector<int>();

		int total_successes = 0;

		for (std::pair<bool, int> pair : result) {
			if (pair.first) {
				total_successes++;
			}
			evaluations.push_back(pair.second);
		}
		return std::make_pair(MathUtility::Median(evaluations),//Median evaluations
			(float)total_successes / result.size());//Success rate
	}

	static std::pair<DataPoint, std::vector<float>> SuccessRatesTest(GeneticAlgorithm<T> base_ga, int test_size, int section_size) {
		std::vector<std::pair<bool, int>> results = EvaluationsTest(base_ga, test_size);
		std::vector<int> sections = std::vector<int>();

		std::vector<float> evaluations = std::vector<float>();

		std::sort(results.begin(), results.end(), [](auto &left, auto &right) {
			return left.second < right.second;
		});

		for (std::pair<bool, int> pair : results) {
			if (pair.first) {
				int floored = pair.second - (pair.second % section_size);
				int index = floored / section_size;
				if (sections.size() < index + 1) {
					sections.resize(index + 1);
				}
				sections[index]++;

				evaluations.push_back(pair.second);
			}
		}

		std::vector<float> success_rates = std::vector<float>();

		//Compute the success rate
		for (int i = 0; i < sections.size(); i++) {
			success_rates.push_back(static_cast<float>(sections[i]) / results.size());
		}

		DataPoint datapoint = evaluations.size() > 0 ? DataPoint(evaluations, true) : DataPoint();

		return std::make_pair(datapoint,//Datapoint of the successful executions
			success_rates);//Distribution of executions (sum = success rate)
	}

	static GeneticAlgorithm<T> SelectBestConfiguration(std::vector<GeneticAlgorithm<T>> gas, int base_test_size, float test_size_increase_rate, float tournament_elimination) {

		std::vector<TestSuite::PartialTestResult<T>> tournament = std::vector<TestSuite::PartialTestResult<T>>();
		int tournament_size = gas.size();
		int test_size = base_test_size;
		int tournament_number = 1;

		//tournament.resize requires an object to use in case the vector's size is increased (which is not going to happen)
		TestSuite::PartialTestResult<T> empty_test = PartialTestResult<T>(PartialTestResult<T>(gas[0], 0));

		if (tournament_size == 1) {
			return gas[0];
		}

		while (tournament_size > 1) {
			std::cout << "Running tournament n. " << std::to_string(tournament_number) << ": " << std::to_string(test_size) << " tests for " << std::to_string(tournament_size) << " configurations" << std::endl;

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