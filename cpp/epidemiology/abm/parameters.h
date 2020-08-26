#ifndef EPI_ABM_PARAMETERS_H
#define EPI_ABM_PARAMETERS_H

#include <limits>

namespace epi
{

struct GlobalInfectionParameters {
    double incubation_time = 1;
};

struct LocalInfectionParameters {
    double susceptible_to_exposed_by_carrier  = 1;
    double susceptible_to_exposed_by_infected = 1;
    double carrier_to_infected                = 1;
    double carrier_to_recovered               = 1;
    double infected_to_recovered              = 1;
    double infected_to_dead                   = 1;
    double recovered_to_susceptible           = 1;
    double detect_infection                   = 0.5;
    double death_factor                       = 1;
    double effective_contacts                 = std::numeric_limits<double>::max();
};

} // namespace epi
#endif