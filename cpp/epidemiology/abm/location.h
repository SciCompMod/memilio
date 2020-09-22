#ifndef EPI_ABM_LOCATION_H
#define EPI_ABM_LOCATION_H

#include "epidemiology/abm/parameters.h"
#include "epidemiology/abm/state.h"
#include "epidemiology/abm/location_type.h"

#include <Eigen/Core>
#include <array>
#include <random>

namespace epi
{
class Person;

/**
 * all locations in the simulated world where persons gather.
 */
class Location
{
public:
    /**
     * construct a Location of a certain type.
     * @param type the type of the location
     */
    Location(LocationType type);

    /** 
     * a person interacts with the population at this location, may change infection state.
     * @param person the person that interacts with the population
     * @param dt length of the current simulation time step
     * @param global_params global infection parameters
     * @return new infection state of the person
     */
    InfectionState interact(const Person& person, double dt, const GlobalInfectionParameters& global_params) const;

    /** 
     * add a person to the population at this location.
     * @param person the person arriving
    */
    void add_person(const Person& person);

    /** 
     * remove a person from the population of this location.
     * @param person the person leaving
     */
    void remove_person(const Person& person);

    /** 
     * notification that one person in this location changed infection state.
     * @param person the person that changed infection state
     * @param old_state the previous infection state of the person
     */
    void changed_state(const Person& person, InfectionState old_state);

    /** 
     * prepare the location for the next simulation step.
     * @param dt the duration of the simulation step
     * @param global_params global infection parameters
     */
    void begin_step(double dt, const GlobalInfectionParameters& global_params);

    /** 
     * number of persons at this location in one infection state.
     * @return number of persons at this location that are in the specified infection state
     */
    int get_subpopulation(InfectionState s) const;

    /** 
     * number of persons at this location for all infection states.
     * vector is indexed by InfectionState.
     * @return number of persons in all infection states.
     * */
    Eigen::Ref<const Eigen::VectorXi> get_subpopulations() const;

    /**
     * @return parameters of the infection that are specific to this location
     */
    LocalInfectionParameters& get_infection_parameters()
    {
        return m_parameters;
    }
    const LocalInfectionParameters& get_infection_parameters() const
    {
        return m_parameters;
    }

private:
    void change_subpopulation(InfectionState s, int delta);

private:
    LocationType m_type;
    int m_num_persons = 0;
    std::array<int, size_t(InfectionState::Count)> m_subpopulations;
    LocalInfectionParameters m_parameters;
    double m_cached_exposure_rate;
};
} // namespace epi

#endif
