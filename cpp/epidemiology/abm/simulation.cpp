#include "epidemiology/abm/simulation.h"

namespace epi
{

AbmSimulation::AbmSimulation(double t, World&& world)
    : m_world(std::move(world))
    , m_t(t)
    , m_dt(0.1)
    , m_result()
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
    m_result.emplace_back(t, ResultVector::Zero());
    for (auto&& node : m_world.get_nodes())
    {
        for (size_t i = 0; i < size_t(InfectionState::Count); i++)
        {
            m_result.back().second[i] += node.get_subpopulation(InfectionState(i));
        }
    }
}

} // namespace epi