
#include "ide_endemic_secir/simulation.h"
#include "ide_endemic_secir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include <memory>

namespace mio
{
namespace endisecir
{

void Simulation::advance(ScalarType tmax)
{
    mio::log_info("Simulating IDE-END-SECIR from t0 = {} until tmax = {} with dt = {}.",
                  m_model->transitions.get_last_time(), tmax, m_dt);

    m_model->compparameters->set_transitiondistributions_support_max(m_dt);
    m_model->compparameters->set_transitiondistributions(m_dt);
    m_model->compparameters->set_transitiondistributions_derivative(m_dt);
    m_model->compparameters->set_infectivity(m_dt);
    m_model->compparameters->set_reproductionnumber_c(m_dt);
    m_model->compparameters->set_FoI_0();
    m_model->compparameters->set_InitFoI(m_dt);

    // For every time step:
    while (m_model->transitions.get_last_time() < tmax - m_dt / 2) {

        m_model->transitions.add_time_point(m_model->transitions.get_last_time() + m_dt);
        m_model->transitions_update.add_time_point(m_model->transitions_update.get_last_time() + m_dt);
        m_model->populations.add_time_point(m_model->populations.get_last_time() + m_dt);
        m_model->populations_update.add_time_point(m_model->populations_update.get_last_time() + m_dt);
        m_model->m_forceofinfection.add_time_point(m_model->m_forceofinfection.get_last_time() + m_dt);
        m_model->m_forceofinfectionupdate.add_time_point(m_model->m_forceofinfectionupdate.get_last_time() + m_dt);
        m_model->m_totalpopulation.add_time_point(m_model->m_totalpopulation.get_last_time() + m_dt);
        m_model->m_totalpopulationupdate.add_time_point(m_model->m_totalpopulationupdate.get_last_time() + m_dt);
        m_model->m_normalizedpopulations.add_time_point(m_model->m_normalizedpopulations.get_last_time() + m_dt);

        // Compute susceptibles:
        m_model->compute_susceptibles(m_dt);

        // Compute flows:
        m_model->flows_currents_timestep(m_dt);

        // Compute remaining compartments:
        m_model->update_compartments(m_dt);

        // Compute m_populationsize:
        m_model->compute_populationsize();

        // Compute normalized compartments:
        m_model->compute_normalizedcompartments();

        // Compute m_forceofinfection;
        m_model->compute_forceofinfection(m_dt);
    }
}

TimeSeries<ScalarType> simulate(ScalarType tmax, ScalarType dt, Model const& m_model)
{
    m_model.check_constraints();
    Simulation sim(m_model, dt);
    sim.advance(tmax);
    return sim.get_compartments();
}

} // namespace endisecir

} // namespace mio
