#ifndef EPI_ABM_SIMULATOR_H
#define EPI_ABM_SIMULATOR_H

#include "epidemiology/abm/world.h"

namespace epi
{

class AbmSimulation
{
public:
    AbmSimulation(double t, const World& world);
    void advance(double tmax);
    //results,...

private:
    double m_t;
    double m_dt;
    World m_world;
};

} // namespace epi

#endif