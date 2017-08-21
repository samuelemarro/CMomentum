#pragma once
static class FastRand {
public:
	//Used to seed the generator.
	static void Seed(unsigned int seed);
	static int RandomInt(int max);
	static int RandomInt(int min, int max);
	static float RandomFloat();
	static float RandomFloat(float max);
	static float RandomFloat(float min, float max);
	static int PolinomialInt(float probability, int max);
	static int PolinomialInt(float probability, int min, int max);
};