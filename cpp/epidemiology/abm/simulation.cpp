#include "epidemiology/abm/simulation.h"

namespace epi
{

AbmSimulation::AbmSimulation(double t, World&& world)
    : m_world(std::move(world))
    , m_t(t)
    , m_dt(0.1)
    , m_result(Eigen::Index(InfectionState::Count))
{
    store_result_at(t);
}

void AbmSimulation::advance(double tmax)
{
    auto t = m_t;
    while (t < tmax) { //TODO: FP unstable
        auto dt = std::min(m_dt, tmax - t);
        m_world.evolve(dt);
        store_result_at(t);
        t += m_dt;
    }
}

void AbmSimulation::store_result_at(double t)
{
    m_result.add_time_point(t);
    m_result.get_last_value().setZero();
    for (auto&& node : m_world.get_nodes())
    {
        m_result.get_last_value() += node.get_subpopulations().cast<double>();
    }
}

} // namespace epi