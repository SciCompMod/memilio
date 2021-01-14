#include "epidemiology/abm/simulation.h"

namespace epi
{

AbmSimulation::AbmSimulation(TimePoint t, World&& world)
    : m_world(std::move(world))
    , m_result(Eigen::Index(InfectionState::Count))
    , m_t(t)
    , m_dt(days(1))
{
    store_result_at(t);
}

void AbmSimulation::advance(TimePoint tmax)
{
    auto t = m_t;
    while (t < tmax) {
        auto dt = std::min(m_dt, tmax - t);
        m_world.evolve(t, dt);
        t += m_dt;
        store_result_at(t);
    }
}

void AbmSimulation::store_result_at(TimePoint t)
{
    m_result.add_time_point(t.days());
    m_result.get_last_value().setZero();
    for (auto&& location : m_world.get_locations()) {
        m_result.get_last_value() += location.get_subpopulations().cast<double>();
    }
}

} // namespace epi