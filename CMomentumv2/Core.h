#pragma once
#include <functional>
#include <map>
#include <vector>
#include <memory>
#include <string>
#include <limits>
#include "json.hpp"

using nlohmann::json;


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

class GeneticAlgorithmParameters {
public:
	int chromosome_length_ = 64;
	int population_size_ = 100;
	float mutation_probability_ = 0.05f;
	float crossover_probability_ = 0.6f;
	float elitism_rate_ = 0.1f;
	float recombination_rate_ = 0.2f;
	std::map<std::string, float> additional_parameters_ = std::map<std::string, float>();
	int max_fitness_evaluations_ = -1;
	float max_generation_ = -1;
	float max_stagnation_ = -1;
	float target_fitness_ = FLT_MAX;

	std::string ToString() {
		std::string parameters = "";

		parameters += "Chromosome Length: " + std::to_string(chromosome_length_) + "\n";
		parameters += "Population Size: " + std::to_string(population_size_) + "\n";
		parameters += "Mutation Probability: " + std::to_string(mutation_probability_) + "\n";
		parameters += "Crossover Probability: " + std::to_string(crossover_probability_) + "\n";
		parameters += "Elitism Rate: " + std::to_string(elitism_rate_) + "\n";
		parameters += "Recombination Rate: " + std::to_string(recombination_rate_) + "\n";

		parameters += "Max Fitness Evaluations: " + std::to_string(max_fitness_evaluations_) + "\n";
		parameters += "Max Generation: " + std::to_string(max_generation_) + "\n";
		parameters += "Max Stagnation: " + std::to_string(max_stagnation_) + "\n";
		parameters += "Target Fitness: " + std::to_string(target_fitness_) + "\n";

		if (additional_parameters_.size() != 0) {
			parameters += "Additional Parameters:\n";
		}

		for (auto it = additional_parameters_.begin(); it != additional_parameters_.end(); it++) {
			parameters += "-" + it->first + ": " + std::to_string(it->second) + std::string("\n");
		}

		return parameters;
	}
};

template<typename T>
class GeneticAlgorithm {
public:

	std::function<Chromosome<T>(int, std::map<std::string, float>&)> initialization_ = nullptr;
	std::function<std::pair<Chromosome<T>*, Chromosome<T>*>(std::vector<std::unique_ptr<Chromosome<T>>>&, std::map<std::string, float>&)> selection_ = nullptr;
	std::function<void(Chromosome<T>&, Chromosome<T>&, std::map<std::string, float>&)> crossover_ = nullptr;
	std::function<void(Chromosome<T>&, std::map<std::string, float>&)> mutation_ = nullptr;
	std::function<void(Chromosome<T>&, const Parent<T>&, float, std::map<std::string, float>&)> recombination_ = nullptr;
	std::function<float(Chromosome<T>&, std::map<std::string, float>&)> fitness_function_ = nullptr;

	GeneticAlgorithmParameters parameters_;


	GeneticAlgorithm(
		GeneticAlgorithmParameters parameters,
		std::function<Chromosome<T>(int, std::map<std::string, float>&)> initialization,
		std::function<std::pair<Chromosome<T>*, Chromosome<T>*>(std::vector<std::unique_ptr<Chromosome<T>>>&, std::map<std::string, float>&)> selection,
		std::function<void(Chromosome<T>&, Chromosome<T>&, std::map<std::string, float>&)> crossover,
		std::function<void(Chromosome<T>&, std::map<std::string, float>&)> mutation,
		std::function<void(Chromosome<T>&, const Parent<T>&, float, std::map<std::string, float>&)> recombination,
		std::function<float(Chromosome<T>&, std::map<std::string, float>&)> fitness_function)
		: parameters_(parameters),
		initialization_(initialization),
		selection_(selection),
		crossover_(crossover),
		mutation_(mutation),
		recombination_(recombination),
		fitness_function_(fitness_function)
	{
		CheckFunctions();
	}

	///<summary>
	///Executes the algorithms.
	///<returns> Returns a pair with 1. Whether the execution was successful 2. The amount of evaluations. </summary>
	///</summary>
	std::vector<std::pair<int,float>> RunAlgorithm() {
		CheckFunctions();
		CheckParameters();

		int generation = 1;
		int fitness_evaluations = 0;
		float best_fitness = -FLT_MAX;
		int last_fitness_evaluations = 0;
		int stagnation = 0;

		std::vector<std::pair<int, float>> stored_best_fitness_values = std::vector<std::pair<int, float>>();

		//Create the initial population
		std::vector<std::unique_ptr<Chromosome<T>>> population = std::vector<std::unique_ptr<Chromosome<T>>>();
		population.reserve(parameters_.population_size_);

		for (int i = 0; i < parameters_.population_size_; i++) {
			Chromosome<T>* chromosome = new Chromosome<T>(initialization_(parameters_.chromosome_length_, parameters_.additional_parameters_));
			chromosome->fitness_ = fitness_function_((*chromosome), parameters_.additional_parameters_);
			chromosome->fitness_is_valid_ = true;
			population.push_back(std::unique_ptr<Chromosome<T>>(chromosome));
		}

		fitness_evaluations += parameters_.population_size_;

		//Sort the population by fitness in descending order
		sort(population.begin(), population.end(),
			[](const auto& lhs, const auto& rhs) {
			return lhs->fitness_ > rhs->fitness_;
		});

		stored_best_fitness_values.push_back(std::make_pair(fitness_evaluations, population[0]->fitness_));

		//While the stop criteria aren't satisfied (ignoring unset criteria)
		while ((best_fitness < parameters_.target_fitness_ || parameters_.target_fitness_ == FLT_MAX) &&
			(fitness_evaluations < parameters_.max_fitness_evaluations_ || parameters_.max_fitness_evaluations_ == -1) &&
			(generation < parameters_.max_generation_ || parameters_.max_generation_ == -1) &&
			(stagnation < parameters_.max_stagnation_ || parameters_.max_stagnation_ == -1)) {

			std::vector<std::unique_ptr<Chromosome<T>>> offspring = std::vector<std::unique_ptr<Chromosome<T>>>();
			offspring.reserve(parameters_.population_size_);

			//Compute the elitism size (flooring to multiples of two) 
			int elitism_size = parameters_.population_size_ * parameters_.elitism_rate_;

			if (elitism_size % 2 == 1) {
				elitism_size--;
			}

			for (int i = 0; i < parameters_.population_size_ - elitism_size; i += 2) {
				std::pair<Chromosome<T>*, Chromosome<T>*> parents = selection_(population, parameters_.additional_parameters_);

				//Copy the parents
				Chromosome<T>* child1 = new Chromosome<T>(*parents.first);
				Chromosome<T>* child2 = new Chromosome<T>(*parents.second);

				if (FastRand::RandomFloat() < parameters_.crossover_probability_) {
					crossover_(*child1, *child2, parameters_.additional_parameters_);

					child1->parent1_ = Parent<T>(parents.first->genes_, parents.first->fitness_);
					child1->parent2_ = Parent<T>(parents.second->genes_, parents.second->fitness_);
					child2->parent1_ = Parent<T>(parents.first->genes_, parents.first->fitness_);
					child2->parent2_ = Parent<T>(parents.second->genes_, parents.second->fitness_);

					child1->has_parents_ = true;
					child2->has_parents_ = true;

					child1->fitness_is_valid_ = false;
					child2->fitness_is_valid_ = false;
				}
				offspring.push_back(std::unique_ptr<Chromosome<T>>(child1));
				offspring.push_back(std::unique_ptr<Chromosome<T>>(child2));
			}

			for (auto& chromosome : offspring)
			{
				if (FastRand::RandomFloat() < parameters_.mutation_probability_) {
					mutation_(*chromosome, parameters_.additional_parameters_);
					chromosome->fitness_is_valid_ = false;
				}

				if (chromosome->has_parents_) {
					//If both parents are eligible, pick one randomly
					if (chromosome->fitness_ < chromosome->parent1_.fitness_ && chromosome->fitness_ < chromosome->parent2_.fitness_)
					{
						if (FastRand::RandomInt(2) == 0) {
							recombination_(*chromosome, chromosome->parent1_, parameters_.recombination_rate_, parameters_.additional_parameters_);
						}
						else {
							recombination_(*chromosome, chromosome->parent2_, parameters_.recombination_rate_, parameters_.additional_parameters_);
						}
					}
					else if (chromosome->fitness_ < chromosome->parent1_.fitness_) {
						recombination_(*chromosome, chromosome->parent1_, parameters_.recombination_rate_, parameters_.additional_parameters_);
					}
					else if (chromosome->fitness_ < chromosome->parent2_.fitness_) {
						recombination_(*chromosome, chromosome->parent2_, parameters_.recombination_rate_, parameters_.additional_parameters_);
					}
				}

				if (!chromosome->fitness_is_valid_) {
					chromosome->fitness_ = fitness_function_(*chromosome, parameters_.additional_parameters_);
					chromosome->fitness_is_valid_ = true;

					fitness_evaluations++;
				}
			}

			for (int i = 0; i < elitism_size; i++) {
				//The population is sorted by fitness in descending order
				Chromosome<T>* parent = new Chromosome<T>(*population[i]);
				offspring.push_back(std::unique_ptr<Chromosome<T>>(parent));
			}

			population = std::move(offspring);

			//Sort the population by fitness in descending order
			std::sort(population.begin(), population.end(),
				[](const auto& lhs, const auto& rhs) {
				return lhs->fitness_ > rhs->fitness_;
			});

			if (population[0]->fitness_ > best_fitness) {
				stagnation = 0;
			}
			else {
				stagnation += fitness_evaluations - last_fitness_evaluations;
			}

			last_fitness_evaluations = fitness_evaluations;

			best_fitness = population[0]->fitness_;

			stored_best_fitness_values.push_back(std::make_pair(fitness_evaluations, best_fitness));

			generation++;
		}

		return stored_best_fitness_values;
	}
private:
	void CheckFunctions() {
		if (initialization_ == NULL) {
			throw std::invalid_argument("The parameter \"initialization\" must be defined.");
		}
		if (selection_ == NULL) {
			throw std::invalid_argument("The parameter \"selection\" must be defined.");
		}
		if (crossover_ == NULL && parameters_.crossover_probability_ != 0) {
			throw std::invalid_argument("The parameter \"crossover\" must be defined.");
		}
		if (mutation_ == NULL && parameters_.mutation_probability_ != 0) {
			throw std::invalid_argument("The parameter \"mutation\" must be defined.");
		}
		if (recombination_ == NULL && parameters_.recombination_rate_ != 0) {
			throw std::invalid_argument("The parameter \"recombination\" must be defined.");
		}
		if (fitness_function_ == NULL) {
			throw std::invalid_argument("The parameter \"fitness_function\" must be defined.");
		}
	}
	void CheckParameters() {
		if (parameters_.population_size_ % 2 != 0 || parameters_.population_size_ < 2) {
			throw std::invalid_argument("The parameter \"population_size_\" must be an even number bigger or equal to 2.");
		}
		if (parameters_.target_fitness_ == FLT_MAX && parameters_.max_fitness_evaluations_ == -1 && parameters_.max_generation_ == -1 && parameters_.max_stagnation_ == -1) {
			throw std::invalid_argument("At least one of the following stop conditions must be set: \"target_fitness_\", \"max_fitness_evaluations_\", \"max_generation_\", \"max_stagnation_\".");
		}
	}
};

void to_json(json& j, const GeneticAlgorithmParameters& gap) {
	j = json{
		{"chromosome_length", gap.chromosome_length_},
		{"population_size", gap.population_size_},
		{"mutation_probability", gap.mutation_probability_},
		{"crossover_probability", gap.crossover_probability_},
		{"elitism_rate", gap.elitism_rate_},
		{"recombination_rate", gap.recombination_rate_},
		{"additional_parameters", gap.additional_parameters_},
		{"max_fitness_evaluations", gap.max_fitness_evaluations_},
		{"max_generation", gap.max_generation_},
		{"max_stagnation", gap.max_stagnation_},
		{"target_fitness", gap.target_fitness_}
	};
}

void from_json(const json& j, GeneticAlgorithmParameters& gap) {
	gap.chromosome_length_ = j.at("chromosome_length").get<int>();
	gap.population_size_ = j.at("population_size").get<int>();
	gap.mutation_probability_ = j.at("mutation_probability").get<float>();
	gap.crossover_probability_ = j.at("crossover_probability").get<float>();
	gap.elitism_rate_ = j.at("elitism_rate").get<float>();
	gap.recombination_rate_ = j.at("recombination_rate").get<float>();
	gap.additional_parameters_ = j.at("additional_parameters").get<std::map<std::string, float>>();
	gap.max_fitness_evaluations_ = j.at("max_fitness_evaluations").get<int>();
	gap.max_generation_ = j.at("max_generation").get<int>();
	gap.max_stagnation_ = j.at("max_stagnation").get<int>();
	gap.target_fitness_ = j.at("target_fitness").get<float>();
}
