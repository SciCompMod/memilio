
#include "ide_endemic_secir/simulation.h"
#include "ide_endemic_secir/computed_parameters.h"
#include "ide_endemic_secir/model.h"
#include "ide_endemic_secir/normalized_model.h"
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
    m_model->compparameters->set_B(m_dt);
    m_model->compparameters->set_infectivity(m_dt);
    m_model->compparameters->set_FoI_0(m_dt);
    m_model->compparameters->set_InitFoI(m_dt);
    m_model->compparameters->set_reproductionnumber_c(m_dt);
    m_model->compparameters->set_T(m_dt);
    m_model->compparameters->set_W(m_dt);
    m_model->initialization_compute_forceofinfection();

    m_normmodel->compparameters->set_transitiondistributions_support_max(m_dt);
    m_normmodel->compparameters->set_transitiondistributions(m_dt);
    m_normmodel->compparameters->set_transitiondistributions_derivative(m_dt);
    m_normmodel->compparameters->set_B(m_dt);
    m_normmodel->compparameters->set_infectivity(m_dt);
    m_normmodel->compparameters->set_FoI_0(m_dt);
    m_normmodel->compparameters->set_InitFoI(m_dt);
    m_normmodel->compparameters->set_reproductionnumber_c(m_dt);
    m_normmodel->compparameters->set_T(m_dt);
    m_normmodel->compparameters->set_W(m_dt);
    m_normmodel->initialization_compute_forceofinfection();

    // For every time step:
    while (m_model->transitions.get_last_time() < tmax - m_dt / 2) {

        //standard model:
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

        // normalized model:
        m_normmodel->transitions.add_time_point(m_normmodel->transitions.get_last_time() + m_dt);
        m_normmodel->populations.add_time_point(m_normmodel->populations.get_last_time() + m_dt);
        m_normmodel->m_forceofinfection.add_time_point(m_normmodel->m_forceofinfection.get_last_time() + m_dt);

        // Compute susceptibles:
        m_normmodel->compute_susceptibles(m_dt);

        // Compute flows:
        m_normmodel->flows_currents_timestep(m_dt);

        // Compute remaining compartments:
        m_normmodel->update_compartments(m_dt);

        // Compute m_forceofinfection;
        m_normmodel->compute_forceofinfection(m_dt);
    }
}

TimeSeries<ScalarType> simulate(ScalarType tmax, ScalarType dt, CompParameters const& m_compparams,
                                Model const& m_model, NormModel const& m_normmodel)
{
    m_model.check_constraints();
    Simulation sim(m_compparams, m_model, m_normmodel, dt);
    sim.advance(tmax);
    return sim.get_compartments();
}

} // namespace endisecir

} // namespace mio
