#include "epidemiology/abm/rng.h"

namespace epi
{

std::mt19937_64& thread_local_rng()
{
    static thread_local auto rng = std::mt19937_64(std::random_device()());
    return rng;
}

}