#pragma once
#include <functional>
#include <map>
#include <vector>
#include <memory>
#include <string>

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

template<typename T>
class GeneticAlgorithm {
public:

	std::function<Chromosome<T>(int, std::map<std::string, float>&)> initialization_ = nullptr;
	std::function<std::pair<Chromosome<T>*, Chromosome<T>*>(std::vector<std::unique_ptr<Chromosome<T>>>&, std::map<std::string, float>&)> selection_ = nullptr;
	std::function<void(Chromosome<T>&, Chromosome<T>&, std::map<std::string, float>&)> crossover_ = nullptr;
	std::function<void(Chromosome<T>&, std::map<std::string, float>&)> mutation_ = nullptr;
	std::function<void(Chromosome<T>&, const Parent<T>&, float, std::map<std::string, float>&)> recombination_ = nullptr;
	std::function<float(Chromosome<T>&, std::map<std::string, float>&)> fitness_function_ = nullptr;

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
		std::function<Chromosome<T>(int, std::map<std::string, float>&)> initialization,
		std::function<std::pair<Chromosome<T>*, Chromosome<T>*>(std::vector<std::unique_ptr<Chromosome<T>>>&, std::map<std::string, float>&)> selection,
		std::function<void(Chromosome<T>&, Chromosome<T>&, std::map<std::string, float>&)> crossover,
		std::function<void(Chromosome<T>&, std::map<std::string, float>&)> mutation,
		std::function<void(Chromosome<T>&, const Parent<T>&, float, std::map<std::string, float>&)> recombination,
		std::function<float(Chromosome<T>&, std::map<std::string, float>&)> fitness_function)
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
		std::vector<std::unique_ptr<Chromosome<T>>> population = std::vector<std::unique_ptr<Chromosome<T>>>();
		population.reserve(population_size_);

		for (int i = 0; i < population_size_; i++) {
			Chromosome<T>* chromosome = new Chromosome<T>(initialization_(chromosome_length_, additional_parameters_));
			chromosome->fitness_ = fitness_function_((*chromosome), additional_parameters_);
			chromosome->fitness_is_valid_ = true;
			population.push_back(std::unique_ptr<Chromosome<T>>(chromosome));
		}

		//Sort the population by fitness in descending order
		sort(population.begin(), population.end(),
			[](const auto& lhs, const auto& rhs) {
			return lhs->fitness_ > rhs->fitness_;
		});

		while ((best_fitness < target_fitness_ || target_fitness_ == INT32_MAX) &&
			(fitness_evaluations < max_fitness_evaluations_ || max_fitness_evaluations_ == -1) &&
			(generation < max_generation_ || max_generation_ == -1)) {

			std::vector<std::unique_ptr<Chromosome<T>>> offspring = std::vector<std::unique_ptr<Chromosome<T>>>();
			offspring.reserve(population_size_);

			//Compute the elitism size (flooring to multiples of two) 
			int elitism_size = population_size_ * elitism_rate_;

			if (elitism_size % 2 == 1) {
				elitism_size--;
			}

			for (int i = 0; i < elitism_size; i++) {
				Chromosome<T>* parent = new Chromosome<T>(*population[i]);
				offspring.push_back(std::unique_ptr<Chromosome<T>>(parent));
			}

			for (int i = 0; i < population_size_ - elitism_size; i += 2) {
				std::pair<Chromosome<T>*, Chromosome<T>*> parents = selection_(population, additional_parameters_);

				//Copy the parents
				Chromosome<T>* child1 = new Chromosome<T>(*parents.first);
				Chromosome<T>* child2 = new Chromosome<T>(*parents.second);

				if (FastRand::RandomFloat() < crossover_probability_) {
					crossover_(*child1, *child2, additional_parameters_);

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
				if (FastRand::RandomFloat() < mutation_probability_) {
					mutation_(*chromosome, additional_parameters_);
					chromosome->fitness_is_valid_ = false;
				}

				if (!chromosome->fitness_is_valid_) {
					chromosome->fitness_ = fitness_function_(*chromosome, additional_parameters_);
					chromosome->fitness_is_valid_ = true;

					fitness_evaluations++;

					if (chromosome->has_parents_) {
						//If both parents are eligible, pick one randomly
						if (chromosome->fitness_ < chromosome->parent1_.fitness_ && chromosome->fitness_ < chromosome->parent2_.fitness_)
						{
							if (FastRand::RandomInt(2) == 0) {
								recombination_(*chromosome, chromosome->parent1_, recombination_rate_, additional_parameters_);
							}
							else {
								recombination_(*chromosome, chromosome->parent2_, recombination_rate_, additional_parameters_);
							}
						}
						else if (chromosome->fitness_ < chromosome->parent1_.fitness_) {
							recombination_(*chromosome, chromosome->parent1_, recombination_rate_, additional_parameters_);
						}
						else if (chromosome->fitness_ < chromosome->parent2_.fitness_) {
							recombination_(*chromosome, chromosome->parent2_, recombination_rate_, additional_parameters_);
						}
					}
				}
			}
			population = std::move(offspring);

			//Sort the population by fitness in descending order
			std::sort(population.begin(), population.end(),
				[](const auto& lhs, const auto& rhs) {
				return lhs->fitness_ > rhs->fitness_;
			});

			best_fitness = population[0]->fitness_;

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
