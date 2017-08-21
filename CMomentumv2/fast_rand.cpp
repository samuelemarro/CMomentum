#include "stdafx.h"
#include "fast_rand.h"
#include <stdlib.h>
#include <random>

thread_local std::default_random_engine engine;
thread_local std::uniform_real_distribution<float> real_distribution = std::uniform_real_distribution<float>(0.0f, 1.0f);

//Used to seed the generator.
void FastRand::Seed(unsigned int seed)
{
	engine = std::default_random_engine(seed);
}
int FastRand::RandomInt(int max) {
	return std::uniform_int_distribution<int>(0, max - 1)(engine);
}
int FastRand::RandomInt(int min, int max) {
	return std::uniform_int_distribution<int>(min, max - 1)(engine);
}

float FastRand::RandomFloat() {
	return real_distribution(engine);
}
float FastRand::RandomFloat(float max) {
	return RandomFloat() * max;
}
float FastRand::RandomFloat(float min, float max) {
	return min + RandomFloat() * (max - min);
}
int FastRand::PolinomialInt(float probability, int max) {
	return std::binomial_distribution<int>(max - 1, probability)(engine);
}
int FastRand::PolinomialInt(float probability, int min, int max) {
	return min + std::binomial_distribution<int>(max - min - 1, probability)(engine);
}
