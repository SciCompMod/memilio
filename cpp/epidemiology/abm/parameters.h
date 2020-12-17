#ifndef EPI_ABM_PARAMETERS_H
#define EPI_ABM_PARAMETERS_H

#include "epidemiology/abm/age.h"
#include "epidemiology/utils/eigen.h"
#include <limits>

namespace epi
{

/**
 * parameter that depends on the age group.
* EXPERIMENTAL; will be merged with new model framework soon.
 */
class AgeDependentParameter
{
public:
    AgeDependentParameter() = default;
    AgeDependentParameter(double d)
        : values(Eigen::VectorXd::Constant(6, d))
    {
    }
    const double& operator[](AbmAgeGroup a) const
    {
        return values[size_t(a)];
    }
    double& operator[](AbmAgeGroup a)
    {
        return values[size_t(a)];
    }
    template <class V, class = std::enable_if_t<std::is_assignable<Eigen::ArrayXd, V>::value, int>>
    AgeDependentParameter& operator=(V&& v)
    {
        values = v;
        return *this;
    }
    const Eigen::ArrayXd& array() const
    {
        return values;
    }
    Eigen::ArrayXd& array()
    {
        return values;
    }

private:
    Eigen::ArrayXd values = Eigen::ArrayXd(6);
};

/**
 * parameters of the infection that are the same everywhere within the world.
 */
class GlobalInfectionParameters {
public:
    AgeDependentParameter incubation_period                  = 1.;
    AgeDependentParameter susceptible_to_exposed_by_carrier  = 1.;
    AgeDependentParameter susceptible_to_exposed_by_infected = 1.;
    AgeDependentParameter carrier_to_infected                = 1.;
    AgeDependentParameter carrier_to_recovered               = 1.;
    AgeDependentParameter infected_to_recovered              = 1.;
    AgeDependentParameter infected_to_dead                   = 1.;
    AgeDependentParameter recovered_to_susceptible           = 1.;
    AgeDependentParameter detect_infection                   = 0.5;
};

/**
 * parameters of the infection that depend on the location.
 */
struct LocalInfectionParameters {
    double death_factor       = 1;
    double effective_contacts = std::numeric_limits<double>::max();
};

} // namespace epi
#endif
