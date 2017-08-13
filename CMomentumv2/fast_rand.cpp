#include "stdafx.h"
#include "fast_rand.h"
#include <stdlib.h>

static const float rand_max_inverse = 1 / RAND_MAX;

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
	return (float)Generate() * rand_max_inverse;
}

//fastrand routine returns one integer, similar output value range as C lib.
inline int FastRand::Generate()
{
	g_seed = (214013 * g_seed + 2531011);

	return (g_seed >> 16) & 0x7FFF;
}