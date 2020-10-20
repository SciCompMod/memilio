#ifndef EPI_ABM_SIMULATOR_H
#define EPI_ABM_SIMULATOR_H

#include "epidemiology/abm/world.h"
#include "epidemiology/utils/time_series.h"

namespace epi
{

/**
 * run the simulation in discrete steps, evolve the world and report results.
 */
class AbmSimulation
{
    using ResultVector = Eigen::Matrix<int, Eigen::Index(InfectionState::Count), 1>;

public:
    /**
     * create a simulation.
     * @param t the starting time of the simulation
     * @param world the world to simulate
     */
    AbmSimulation(double t, World&& world);

    /** 
     * run the simulation from the current time to tmax.
     * @param tmax time to stop
     */
    void advance(double tmax);

    /**
     * get the result of the simulation.
     * sum over all locations of the number of persons in an infection state.
     * @return the result of the simulation.
     */
    const TimeSeries<double>& get_result() const
    {
        return m_result;
    }

private:
    void store_result_at(double t);

    World m_world;
    TimeSeries<double> m_result;
    double m_t;
    double m_dt;
};

} // namespace epi

#endif