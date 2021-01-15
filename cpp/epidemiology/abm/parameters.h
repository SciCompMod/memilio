#ifndef EPI_ABM_PARAMETERS_H
#define EPI_ABM_PARAMETERS_H

#include "epidemiology/abm/age.h"
#include "epidemiology/abm/time.h"
#include "epidemiology/utils/eigen.h"
#include <limits>

namespace epi
{

/**
 * parameter that depends on an enum, e.g. the age group.
 * EXPERIMENTAL; will be merged with new model framework soon.
 * @tparam E enum that is used as an index.
 */
template<class E>
class DependentParameter
{
public:
    /**
     * default constructor.
     */
    DependentParameter() = default;

    /**
     * constant constructor.
     * same value for each group.
     * @param d Value for every group.
     */
    DependentParameter(double d)
        : values(Eigen::VectorXd::Constant(Eigen::Index(E::Count), d))
    {
    }

    /**
     * array constructor.
     * @param a array expression with one value per group.
     * @tparam A array expression.
     */
    template <class A, class = std::enable_if_t<std::is_convertible<Eigen::ArrayBase<Eigen::ArrayXd>, std::decay_t<A>>::value, int>>
    DependentParameter(A&& a)
        : values(std::forward<A>(a))
    {
    }

    /**
     * get the value for one group.
     * @param e group index.
     * @return value for the specified group
     */
    const double& operator[](E e) const
    {
        return values[size_t(e)];
    }
    double& operator[](E e)
    {
        return values[size_t(e)];
    }

    /**
     * set values, constant or array.
     * @param a scalar for every group or array expression with one value for each group
     * @return reference to this object
     * @tparam A scalar or array expression
     */
    template <class A, class = std::enable_if_t<std::is_assignable<Eigen::ArrayXd, A>::value, int>>
    DependentParameter& operator=(A&& a)
    {
        values = std::forward<A>(a);
        return *this;
    }

    /**
     * as Eigen3 array for componentwise operations.
     * @return all values as array expression. 
     */
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

/**
 * parameters that govern the migration between locations
 */
struct MigrationParameters {
    TimePoint lockdown_date = epi::TimePoint(std::numeric_limits<int>::max());
};

} // namespace epi
#endif
