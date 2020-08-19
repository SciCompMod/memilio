#ifndef EPI_ABM_RNG_H
#define EPI_ABM_RNG_H

#include <random>

namespace epi
{

std::mt19937_64& thread_local_rng();

} // namespace epi

#endif