// CMomentum.cpp : definisce il punto di ingresso dell'applicazione console.
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>
#include <iostream>
#include <functional>
#include <map>
#include <memory>



int random_int(int min, int max) {
	return min + rand() / (RAND_MAX / (max - min) + 1);
}

float random_float() {
	return (float)rand() / (float)RAND_MAX;
}

template<typename T>
class Chromosome;

template<typename T>
struct Parent;

template<typename T>
class GeneticAlgorithm {
public:

	std::function<Chromosome<T>(int, std::map<std::string, float>*)> initialization_ = NULL;
	std::function<std::pair<Chromosome<T>, Chromosome<T>>(std::vector<Chromosome<T>>*, std::map<std::string, float>*)> selection_ = NULL;
	std::function<void(Chromosome<T>*, Chromosome<T>*, std::map<std::string, float>*)> crossover_ = NULL;
	std::function<void(Chromosome<T>*, std::map<std::string, float>*)> mutation_ = NULL;
	std::function<void(Chromosome<T>*, const Parent<T>, float, std::map<std::string, float>*)> recombination_ = NULL;
	std::function<float(Chromosome<T>, std::map<std::string, float>*)> fitness_ = NULL;

	std::map<std::string, float>* additional_parameters_ = new std::map<std::string, float>();

	float mutation_probability_ = 0.1f;
	float crossover_probability_ = 0.75f;
	float recombination_rate_ = 0.5f;

	int max_fitness_evaluations_ = INT32_MAX;
	int max_generations_ = INT32_MAX;
	int target_fitness_ = INT32_MAX;
	int chromosome_length_ = 64;
	int population_size_ = 300;

	~GeneticAlgorithm() {
		if (additional_parameters_ != NULL) {
			delete additional_parameters_;
		}
	}

	int RunAlgorithm() {
		CheckParameters();
		int generations = 1;
		int fitness_evaluations = 0;
		float best_fitness = INT32_MIN;

		//Create the initial population
		std::vector<Chromosome<T>> population = std::vector<Chromosome<T>>();

		for (int i = 0; i < population_size_; i++) {
			Chromosome<T> chromosome = initialization_(chromosome_length_, additional_parameters_);
			chromosome.fitness_ = fitness_(chromosome, additional_parameters_);
			chromosome.fitness_is_valid_ = true;
			population.push_back(chromosome);
		}

		while (best_fitness < target_fitness_ && fitness_evaluations < max_fitness_evaluations_ && generations < max_generations_) {

			std::vector<Chromosome<T>> offspring = std::vector<Chromosome<T>>();

			for (int i = 0; i < population_size_; i += 2) {
				std::pair<Chromosome<T>, Chromosome<T>> parents = selection_(&population, additional_parameters_);

				//Copy the parents
				Chromosome<T> child1 = parents.first;// std::shared_ptr<Chromosome<T>>(new Chromosome<T>(*parents.first.get()));
				Chromosome<T> child2 = parents.second;//std::shared_ptr<Chromosome<T>>(new Chromosome<T>(*parents.second.get()));

				if (random_float() < crossover_probability_) {
					crossover_(&child1, &child2, additional_parameters_);

					child1.parent1_ = Parent<T>(parents.first.genes_, parents.first.fitness_);
					child1.parent2_ = Parent<T>(parents.second.genes_, parents.second.fitness_);
					child2.parent1_ = Parent<T>(parents.first.genes_, parents.first.fitness_);
					child2.parent2_ = Parent<T>(parents.second.genes_, parents.second.fitness_);

					child1.has_parents_ = true;
					child2.has_parents_ = true;

					child1.fitness_is_valid_ = false;
					child2.fitness_is_valid_ = false;
				}

				offspring.push_back(child1);
				offspring.push_back(child2);
			}

			for (auto& chromosome : offspring)
			{
				if (random_float() < mutation_probability_) {
					mutation_(&chromosome, additional_parameters_);
					chromosome.fitness_is_valid_ = false;
				}

				if (!chromosome.fitness_is_valid_) {
					chromosome.fitness_ = fitness_(chromosome, additional_parameters_);
					chromosome.fitness_is_valid_ = true;

					fitness_evaluations++;

					if (chromosome.has_parents_) {
						//If both parents are eligible, pick one randomly
						if (chromosome.fitness_ < chromosome.parent1_.fitness_ && chromosome.fitness_ < chromosome.parent2_.fitness_)
						{
							if (random_float() < 0.5f) {
								recombination_(&chromosome, chromosome.parent1_, recombination_rate_, additional_parameters_);
							}
							else {
								recombination_(&chromosome, chromosome.parent2_, recombination_rate_, additional_parameters_);
							}
						}
						else if (chromosome.fitness_ < chromosome.parent1_.fitness_) {
							recombination_(&chromosome, chromosome.parent1_, recombination_rate_, additional_parameters_);
						}
						else if (chromosome.fitness_ < chromosome.parent2_.fitness_) {
							recombination_(&chromosome, chromosome.parent2_, recombination_rate_, additional_parameters_);
						}
					}

				}
			}

			population = offspring;//è sufficiente?

			best_fitness = INT32_MIN;

			for each (Chromosome<T> chromosome in population) {
				if (chromosome.fitness_ > best_fitness) {
					best_fitness = chromosome.fitness_;
				}
			}

			generations++;
		}

		return fitness_evaluations;
	}
private:
	void CheckParameters() {
		if (initialization_ == NULL) {
			throw std::invalid_argument("The parameter \"initialization_\" must be defined.");
		}
		if (selection_ == NULL) {
			throw std::invalid_argument("The parameter \"selection_\" must be defined.");
		}
		if (crossover_ == NULL) {
			throw std::invalid_argument("The parameter \"crossover_\" must be defined.");
		}
		if (mutation_ == NULL) {
			throw std::invalid_argument("The parameter \"mutation_\" must be defined.");
		}
		if (recombination_ == NULL) {
			throw std::invalid_argument("The parameter \"recombination_\" must be defined.");
		}
		if (fitness_ == NULL) {
			throw std::invalid_argument("The parameter \"fitness_\" must be defined.");
		}

		if (population_size_ % 2 != 0) {
			throw std::invalid_argument("The parameter \"population_size_\" must be even.");
		}
		if (target_fitness_ == INT32_MAX && max_fitness_evaluations_ == INT32_MAX && max_generations_ == INT32_MAX) {
			throw std::invalid_argument("At least one of the following parameters must not be INT32_MAX: \"target_fitness_\", \"max_fitness_evaluations_\", \"max_generations_\".");
		}
	}
};

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

void BitFlipMutation(Chromosome<bool>* chromosome, std::map<std::string, float>* additional_parameters) {
	float gene_mutation_rate = (*additional_parameters)["gene_mutation_rate"];

	for (int i = 0; i < chromosome->genes_.size(); i++) {
		if (random_float() < gene_mutation_rate) { 
			chromosome->genes_[i] = !chromosome->genes_[i];
		}
	}
}

template<typename T>
void TwoPointCrossover(Chromosome<T>* parent1, Chromosome<T>* parent2, std::map<std::string, float>* additional_parameters) {
	int pos1 = random_int(0, parent1->genes_.size());
	int pos2 = random_int(0, parent2->genes_.size());

	//Make sure that pos2 > pos1
	if (pos1 > pos2)
	{
		std::swap(pos1, pos2);
	}

	for (int i = pos1; i <= pos2; i++) {
		std::swap(parent1->genes_[i], parent2->genes_[i]);
	}

	return;
}

Chromosome<bool> BinaryInitialization(int length, std::map<std::string, float>* additional_parameters) {
	Chromosome<bool> c = Chromosome<bool>();
	for (int i = 0; i < length; i++) {
		c.genes_.push_back(random_int(0, 2));
	}
	return c;
}

void BinaryRecombinate(Chromosome<bool>* current, const Parent<bool> parent, float recombination_rate, std::map<std::string, float>* additional_parameters) {

	std::vector<int> diff = std::vector<int>();

	for (int i = 0; i < current->genes_.size(); i++) {
		if (current->genes_[i] != parent.genes_[i]) {
			diff.push_back(i);
		}
	}

	int recombinations = static_cast<int>(recombination_rate * diff.size());

	for (int i = 0; i < recombinations; i++) {
		int random_index = random_int(0, diff.size());
		int diff_index = diff[random_index];

		current->genes_[diff_index] = !current->genes_[diff_index];
	}
}

template<typename T>
std::pair<Chromosome<T>, Chromosome<T>> TournamentSelection(std::vector<Chromosome<T>>* population, std::map<std::string, float>* additional_parameters) {

	int tournament_size = static_cast<int>((*additional_parameters)["tournament_size"]);

	//Chromosome<T> first_parent = RunTournament(population, tournament_size);
	//Chromosome<T> second_parent = RunTournament(population, tournament_size);
	//return std::pair<Chromosome<T>, Chromosome<T>>(first_parent, second_parent);
	return std::pair<Chromosome<T>, Chromosome<T>>(RunTournament(population, tournament_size), RunTournament(population, tournament_size));
}

template<typename T>
Chromosome<T> RunTournament(std::vector<Chromosome<T>>* population, int tournament_size) {
	std::vector<Chromosome<T>> tournament = std::vector<Chromosome<T>>();

	for (int i = 0; i < tournament_size; i++) {
		tournament.push_back((*population)[random_int(0, population->size())]);
	}

	//Sort in descending order
	sort(tournament.begin(), tournament.end(),
		[](const auto& lhs, const auto& rhs) {
		return lhs.fitness_ > rhs.fitness_;
	});

	return tournament[0];
}

float OneMaxFitness(Chromosome<bool> chromosome, std::map<std::string, float>* additional_parameters) {
	float fitness = 0;
	for each(bool gene in chromosome.genes_) {
		if (!gene) {
			fitness--;
		}
	}
	return fitness;
}

int RunTest() {
	GeneticAlgorithm<bool> ga = GeneticAlgorithm<bool>();

	ga.initialization_ = BinaryInitialization;
	ga.selection_ = TournamentSelection<bool>;
	ga.crossover_ = TwoPointCrossover<bool>;
	ga.mutation_ = BitFlipMutation;
	ga.recombination_ = BinaryRecombinate;
	ga.fitness_ = OneMaxFitness;

	ga.target_fitness_ = 0;

	(*ga.additional_parameters_)["tournament_size"] = 3;
	(*ga.additional_parameters_)["gene_mutation_rate"] = 0.05f;

	return ga.RunAlgorithm();
}

int main()
{
	srand(time(NULL));
	
	std::vector<int> evaluations = std::vector<int>();
	for (int i = 0; i < 1000; i++) {
		evaluations.push_back(RunTest());
		std::cout << i << std::endl;
	}

	std::sort(evaluations.begin(), evaluations.end());

	float median = 0;

	if (evaluations.size() % 2 == 0) {
		median = (evaluations[(evaluations.size() - 1) / 2] + evaluations[(evaluations.size() - 1) / 2 + 1]) / 2;
	}
	else {
		median = evaluations[(evaluations.size() - 1) / 2];
	}

	std::cout << median << std::endl;
	system("PAUSE");
	return 0;
}

