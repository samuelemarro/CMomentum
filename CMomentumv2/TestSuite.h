#pragma once
#include "stdafx.h"
#include <vector>
#include <thread>
#include <process.h>
#include <random>

template<typename T>
class TestSuite {

	template<typename T>
	class PartialTestResult {
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

				float evaluation = ga.RunAlgorithm();

#pragma omp critical
				{
					results[ga_index].evaluations_.push_back(evaluation);
				}

				executed_tests++;
				int executed_tests_copy = executed_tests; //Local variable used to prevent race conditions
				if (executed_tests_copy % std::max(1, total_tests / 100) == 0) {
					EraseWriteLine("Progress: " + std::to_string(executed_tests_copy * 100 / total_tests) + "%");
				}
			}
		}

		for (int i = 0; i < results.size(); i++)
		{
			results[i].median_ = Median(results[i].evaluations_);
		}

		std::cout << "\n";

		return results;
	}
	static float RunTest(GeneticAlgorithm<T> ga, int test_size) {

		std::vector<int> evaluations = std::vector<int>();

		int executed_tests = 0;
		std::random_device seed_generator;
#pragma omp parallel
		{
#pragma omp for
			for (int i = 1; i <= test_size; i++) {
				FastRand::Seed(seed_generator());
				float evaluation = ga.RunAlgorithm();
#pragma omp critical
				{
					evaluations.push_back(evaluation);
				}

				executed_tests++;
				int executed_tests_copy = executed_tests;//Local variable used to prevent race conditions
				if (executed_tests_copy % std::max(1, test_size / 100) == 0) {
					EraseWriteLine("Progress: " + std::to_string(executed_tests_copy * 100 / test_size) + "%");
				}
			}
		}

		std::cout << "\n";

		return Median(evaluations);
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
private:
	static float Median(std::vector<int> vector) {
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
	static void EraseWriteLine(std::string text) {
		std::cout << text;
		std::cout << std::string(text.length(), '\b');
	}
};