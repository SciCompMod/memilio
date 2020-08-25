#ifndef EPI_ABM_RNG_H
#define EPI_ABM_RNG_H

#include <random>
#include <functional>

namespace epi
{

struct Rng {
    using result_type = std::mt19937_64::result_type;
    static result_type min()
    {
        return std::mt19937_64::min();
    }
    static result_type max()
    {
        return std::mt19937_64::max();
    }
    result_type operator()()
    {
        return generator();
    }
    std::function<result_type()> generator;
};

Rng& thread_local_rng();

} // namespace epi

#endif