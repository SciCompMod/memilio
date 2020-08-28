#ifndef EPI_ABM_RANDOM_NUMBER_GENERATOR_H
#define EPI_ABM_RANDOM_NUMBER_GENERATOR_H

#include <random>
#include <functional>

namespace epi
{

/**
 * abstract random number generator.
 * models UniformRandomBitGenerator concept.
 */
struct RandomNumberGenerator {
    using result_type = std::mt19937_64::result_type;
    static constexpr result_type min()
    {
        return std::mt19937_64::min();
    }
    static constexpr result_type max()
    {
        return std::mt19937_64::max();
    }
    result_type operator()()
    {
        return generator();
    }
    std::function<result_type()> generator;
};

/**
 * get a random number generator that is static and local to this thread.
 * @return a random number generator that is static and local to this thread.
 */
RandomNumberGenerator& thread_local_rng();

} // namespace epi

#endif