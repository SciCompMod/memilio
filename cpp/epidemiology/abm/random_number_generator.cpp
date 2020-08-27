#include "epidemiology/abm/random_number_generator.h"

namespace epi
{

RandomNumberGenerator& thread_local_rng()
{
    static thread_local auto rng = RandomNumberGenerator{std::mt19937_64(std::random_device()())};
    return rng;
}

} // namespace epi