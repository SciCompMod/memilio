#ifndef EPI_ABM_PARAMETERS_H
#define EPI_ABM_PARAMETERS_H

#include "epidemiology/abm/age.h"
#include "epidemiology/utils/eigen.h"
#include <limits>

namespace epi
{

/**
 * parameter that depends an enum, e.g. the age group.
* EXPERIMENTAL; will be merged with new model framework soon.
 */
template<class E>
class DependentParameter
{
public:
    DependentParameter() = default;
    DependentParameter(double d)
        : values(Eigen::VectorXd::Constant(Eigen::Index(E::Count), d))
    {
    }
    const double& operator[](E a) const
    {
        return values[size_t(a)];
    }
    double& operator[](E a)
    {
        return values[size_t(a)];
    }
    template <class V, class = std::enable_if_t<std::is_assignable<Eigen::ArrayXd, V>::value, int>>
    DependentParameter& operator=(V&& v)
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
    Eigen::ArrayXd values { Eigen::Index(E::Count) };
};

/**
 * parameters of the infection that are the same everywhere within the world.
 */
class GlobalInfectionParameters {
public:
    DependentParameter<AbmAgeGroup> incubation_period                  = 1.;
    DependentParameter<AbmAgeGroup> susceptible_to_exposed_by_carrier  = 1.;
    DependentParameter<AbmAgeGroup> susceptible_to_exposed_by_infected = 1.;
    DependentParameter<AbmAgeGroup> carrier_to_infected                = 1.;
    DependentParameter<AbmAgeGroup> carrier_to_recovered               = 1.;
    DependentParameter<AbmAgeGroup> infected_to_recovered              = 1.;
    DependentParameter<AbmAgeGroup> infected_to_dead                   = 1.;
    DependentParameter<AbmAgeGroup> recovered_to_susceptible           = 1.;
    DependentParameter<AbmAgeGroup> detect_infection                   = 0.5;
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
