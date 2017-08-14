#include "stdafx.h"
#include "fast_rand.h"
#include <stdlib.h>

static unsigned int g_seed;

//Used to seed the generator.
void FastRand::Seed(unsigned int seed)
{
	g_seed = seed;
}
int FastRand::RandomInt(int max) {
	return Generate() % max;
}
int FastRand::RandomInt(int min, int max) {
	return min + Generate() % (max - min);
}

float FastRand::RandomFloat() {
	return (float)Generate() / RAND_MAX;
}
float FastRand::RandomFloat(float max) {
	return RandomFloat() * max;
}
float FastRand::RandomFloat(float min, float max) {
	return min + RandomFloat() * (max - min);
}

//fastrand routine returns one integer, similar output value range as C lib.
inline int FastRand::Generate()
{
	g_seed = (214013 * g_seed + 2531011);

	return (g_seed >> 16) & 0x7FFF;
}