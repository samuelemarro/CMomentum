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
private:
	//fastrand routine returns one integer, similar output value range as C lib.
	static inline int Generate();
};