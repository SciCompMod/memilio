#ifndef IDE_END_SECIR_SIMULATION_H
#define IDE_END_SECIR_SIMULATION_H

#include "ide_endemic_secir/computed_parameters.h"
#include "ide_endemic_secir/model.h"
#include "ide_endemic_secir/normalized_model.h"
#include "memilio/config.h"
#include "memilio/utils/miompi.h"
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
    Simulation(CompParameters const& compparams, Model const& model, NormModel const& normmodel, ScalarType dt = 0.1)
        : m_compparams(std::make_unique<CompParameters>(compparams))
        , m_model(std::make_unique<Model>(model))
        , m_normmodel(std::make_unique<NormModel>(normmodel))
        , m_dt(dt)
    {
    }

    /**
     * @brief Run the simulation until time tmax.
     */
    void advance(ScalarType tmax);

    /**
     * @brief Get the result of the simulation for the compartments of m_model
     * Return the number of persons in all #InfectionState%s
     */
    TimeSeries<ScalarType> get_compartments()
    {
        return m_model->populations;
    }
    /**
     * @brief Get the result of the simulation for the compartments of m_model, where we use the update scheme.
     * Return the number of persons in all #InfectionState%s
     */
    TimeSeries<ScalarType> get_compartments_update()
    {
        return m_model->populations_update;
    }

    /**
     * @brief Get the result of the simulation for the normalized compartments, where we use m_model and the compartments
     * computed using the sum scheme.
     * Return the number of persons in all #InfectionState%s
     */
    TimeSeries<ScalarType> get_normalizedcompartments()
    {
        return m_model->m_normalizedpopulations;
    }
    /**
     * @brief Get the result of the simulation for the compartments of the normalized model.
     * Return the number of persons in all #InfectionState%s
     */
    TimeSeries<ScalarType> get_normmodel_compartments()
    {
        return m_normmodel->populations;
    }

    /**
     * @brief Get the result of the simulation for the compartments.
     * Return the number of persons in all #InfectionState%s
     */
    TimeSeries<ScalarType>& get_compartments() const
    {
        return m_model->populations;
    }

    /**
     * @brief Get the result of the simulation for the compartments of m_model, where we use the update scheme.
     * Return the number of persons in all #InfectionState%s
     */
    TimeSeries<ScalarType>& get_compartments_update() const
    {
        return m_model->populations_update;
    }

    /**
     * @brief Get the result of the simulation for the normalized compartments, where we use m_model.
     * Return the number of persons in all #InfectionState%s
     */
    TimeSeries<ScalarType>& get_normalizedcompartments() const
    {
        return m_model->m_normalizedpopulations;
    }

    /**
     * @brief Get the result of the simulation for the compartments of the normalized model.
     * Return the number of persons in all #InfectionState%s
     */
    TimeSeries<ScalarType>& get_normmodel_compartments() const
    {
        return m_normmodel->populations;
    }

    /**
     * @brief Get the transitions between the different #InfectionState%s for m_model.
     */
    TimeSeries<ScalarType> const& get_transitions()
    {
        return m_model->transitions;
    }
    /**
     * @brief Get the transitions between the different #InfectionState%s for m_model using the update scheme.
     */
    TimeSeries<ScalarType> const& get_transitions_update()
    {
        return m_model->transitions_update;
    }

    /**
     * @brief Get the transitions between the different #InfectionState%s  of m_normmodel.
     */
    TimeSeries<ScalarType> const& get_normmodel_transitions()
    {
        return m_normmodel->transitions;
    }

    /**
     * @brief Get the force of infection term of m_model.
     */
    TimeSeries<ScalarType> const& get_forceofinfections()
    {
        return m_model->m_forceofinfection;
    }

    /**
     * @brief Get the force of infection term of m_model using the update scheme.
     */
    TimeSeries<ScalarType> const& get_forceofinfections_update()
    {
        return m_model->m_forceofinfectionupdate;
    }

    /**
     * @brief Get the force of infection term of m_normmodel.
     */
    TimeSeries<ScalarType> const& get_normmodel_forceofinfections()
    {
        return m_normmodel->m_forceofinfection;
    }

    /**
     * @brief Get the total population of m_model.
     */
    TimeSeries<ScalarType> const& get_totalpopulations()
    {
        return m_model->m_totalpopulation;
    }

    /**
     * @brief Get the total population of m_model using the update scheme.
     */
    TimeSeries<ScalarType> const& get_totalpopulations_update()
    {
        return m_model->m_totalpopulationupdate;
    }

    /**
     * @brief Get the derivative of the total population of m_model.
     */
    TimeSeries<ScalarType> const& get_totalpopulations_derivative()
    {
        return m_model->m_totalpopulation_derivative;
    }

    /**
     * @brief Get the reproduction numer.
     */
    ScalarType const& get_reproductionnumber_c()
    {
        return m_model->compparameters->m_reproductionnumber_c;
    }

    /**
     * @brief Get T.
     */
    std::vector<ScalarType> const& get_T()
    {
        return m_model->compparameters->m_T;
    }

    /**
     * @brief Get V.
     */
    std::vector<ScalarType> const& get_V()
    {
        return m_model->compparameters->m_V;
    }

    /**
     * @brief Get W.
     */
    std::vector<ScalarType> const& get_W()
    {
        return m_model->compparameters->m_W;
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
     * @brief returns the simulation normmodel used in simulation.
     */
    const NormModel& get_normmodel() const
    {
        return *m_normmodel;
    }

    /**
      * @brief returns the simulation normmodel used in simulation.
      */
    NormModel& get_normmodel()
    {
        return *m_normmodel;
    }

    /** 
     * @brief get the size of the time step of the simulation.
     */
    ScalarType get_stepsize()
    {
        return m_dt;
    }

    void compute_difference_normalizations()
    {
        for (int state = 0; state < (int)InfectionState::Count - 1; state++) {
            m_difference_normalizedcompartments.get_last_value()[state] =
                std::abs(m_normmodel->populations.get_last_value()[state] -
                         m_model->m_normalizedpopulations.get_last_value()[state]);
        }
        m_difference_normalizedFoI.get_last_value()[0] = std::abs(m_normmodel->m_forceofinfection.get_last_value()[0] -
                                                                  m_model->m_forceofinfection.get_last_value()[0]);
    }

    TimeSeries<ScalarType> const& get_difference_normalizationcomp()
    {
        return m_difference_normalizedcompartments;
    }

    TimeSeries<ScalarType> const& get_difference_normalizationFoi()
    {
        return m_difference_normalizedFoI;
    }

private:
    std::unique_ptr<CompParameters> m_compparams; ///< Unique pointer to the computed Parameters.
    std::unique_ptr<Model> m_model; ///< Unique pointer to the simulated Model.
    std::unique_ptr<NormModel> m_normmodel; ///< Unique pointer to the simulated normalized Model.
    ScalarType m_dt; ///< Time step used for numerical computations in simulation.
    TimeSeries<ScalarType> m_difference_normalizedcompartments{TimeSeries<ScalarType>(
        (int)InfectionState::Count - 1)}; ///< TimeSeries containing the difference of the compartments
    // computed by NormModel and the normalized compartments computed in Model.
    TimeSeries<ScalarType> m_difference_normalizedFoI{
        TimeSeries<ScalarType>(1)}; ///< TimeSeries containing the difference of the force of infection terms
    // computed by NormModel and Model.
};

TimeSeries<ScalarType> simulate(ScalarType tmax, ScalarType dt, Model const& model);

} // namespace endisecir
} // namespace mio

#endif //IDE_END_SECIR_SIMULATION_H