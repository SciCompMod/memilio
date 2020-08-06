#include "epidemiology/abm/simulation.h"

namespace epi
{

AbmSimulation::AbmSimulation(double t, const World& world)
    : m_world(world)
    , m_t(t)
    , m_dt(0.1)
{
}

void AbmSimulation::advance(double tmax)
{
    auto t = m_t;
    while (t < tmax) { //TODO: FP unstable
        auto dt = std::min(m_dt, tmax - t);
        m_world.evolve(dt);
        //TODO: save world statistics for result
        t += m_dt;
    }
}

} // namespace epi