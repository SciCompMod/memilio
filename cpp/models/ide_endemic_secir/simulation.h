#ifndef IDE_END_SECIR_SIMULATION_H
#define IDE_END_SECIR_SIMULATION_H

#include "ide_endemic_secir/model.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include <memory>
#include <vector>

namespace mio
{
namespace endisecir
{

class Simulation
{
public:
    Simulation(Model const& model, ScalarType dt = 0.1)
        : m_model(std::make_unique<Model>(model))
        , m_dt(dt)
    {
    }

    /**
     * @brief Run the simulation until time tmax.
     */
    void advance(ScalarType tmax);

    /**
     * @brief Get the result of the simulation for the compartments
     * Return the number of persons in all #InfectionState%s
     */
    TimeSeries<ScalarType> get_compartments()
    {
        return m_model->populations;
    }

    TimeSeries<ScalarType> get_normalizedcompartments()
    {
        return m_model->m_normalizedpopulations;
    }

    /**
     * @brief Get the result of the simulation for the compartments
     * Return the number of persons in all #InfectionState%s
     */
    TimeSeries<ScalarType>& get_compartments() const
    {
        return m_model->populations;
    }

    TimeSeries<ScalarType>& get_normalizedcompartments() const
    {
        return m_model->m_normalizedpopulations;
    }

    /**
     * @brief Get the transitions between the different #InfectionState%s.
     */
    TimeSeries<ScalarType> const& get_transitions()
    {
        return m_model->transitions;
    }

    TimeSeries<ScalarType> const& get_forceofinfections()
    {
        return m_model->m_forceofinfection;
    }

    TimeSeries<ScalarType> const& get_totalpopulations()
    {
        return m_model->m_totalpopulation;
    }

    ScalarType const& get_reproductionnumber_c()
    {
        return m_model->m_reproductionnumber_c;
    }

    std::vector<ScalarType> const& get_equilibriumcompartments()
    {
        return m_model->m_equilibriumnormalizedcompartments;
    }

    ScalarType const& get_equilibrium_forceofinfection()
    {
        return m_model->m_equilibriumforceofinfection;
    }

    /**
     * @brief returns the simulation model used in simulation.
     */
    const Model& get_model() const
    {
        return *m_model;
    }

    /**
     * @brief returns the simulation model used in simulation.
     */
    Model& get_model()
    {
        return *m_model;
    }

    /** 
     * @brief get the size of the time step of the simulation.
     */
    ScalarType get_stepsize()
    {
        return m_dt;
    }

private:
    std::unique_ptr<Model> m_model; ///< Unique pointer to the simulated Model.
    ScalarType m_dt; ///< Time step used for numerical computations in simulation.
};

TimeSeries<ScalarType> simulate(ScalarType tmax, ScalarType dt, Model const& model);

} // namespace endisecir
} // namespace mio

#endif //IDE_END_SECIR_SIMULATION_H