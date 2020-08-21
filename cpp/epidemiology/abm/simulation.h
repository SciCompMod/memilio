#ifndef EPI_ABM_SIMULATOR_H
#define EPI_ABM_SIMULATOR_H

#include "epidemiology/abm/world.h"
#include <Eigen/Core>

namespace epi
{

class AbmSimulation
{
    using ResultVector = Eigen::Matrix<int, Eigen::Index(InfectionState::Count), 1>;
public:
    AbmSimulation(double t, World&& world);
    void advance(double tmax);

    auto get_result()
    {
        return m_result;
    }

private:
    void store_result_at(double t);

    double m_t;
    double m_dt;
    World m_world;
    std::vector<std::pair<double, ResultVector>> m_result; //will be replaced with TimeSeries<double>
};

} // namespace epi

#endif