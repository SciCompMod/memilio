#include "epidemiology/abm/rng.h"

namespace epi
{

Rng& thread_local_rng()
{
    static thread_local auto rng = Rng{std::mt19937_64(std::random_device()())};
    return rng;
}

}