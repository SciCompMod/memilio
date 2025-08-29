#ifndef IDE_END_SECIR_MODEL_H
#define IDE_END_SECIR_MODEL_H

#include "ide_endemic_secir/infection_state.h"
#include "ide_endemic_secir/parameters.h"
#include "ide_endemic_secir/computed_parameters.h"
#include "memilio/config.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/time_series.h"

#include "vector"
#include <Eigen/src/Core/util/Meta.h>
#include <memory>
#include <vector>

namespace mio
{
namespace endisecir
{
// Forward declaration of friend classes/functions of Model.
class Model;
class Simulation;

class Model
{
    using ParameterSet = Parameters;

public:
    /**
    * @brief Constructor to create an endemic IDE_SECIR model.
    *
    * @param[in] TODO!!!!!!
    */

    Model(CompParameters const& compparams);

    /**
    * @brief Checks constraints on model parameters and initial data.
    * @return Returns true if one (or more) constraint(s) are not satisfied, otherwise false.
    */
    bool check_constraints() const;

    /**
     * @brief Setter for the tolerance used to calculate the maximum support of the TransitionDistributions.
     *
     * @param[in] new_tol New tolerance.
     */

    // ---- Public parameters. ----
    std::shared_ptr<CompParameters> compparameters;
    TimeSeries<ScalarType>
        transitions; ///< TimesSeries containing points of time and the corresponding number of individuals transitioning
    // from one #InfectionState to another as defined in #Infection%s.
    TimeSeries<ScalarType>
        transitions_update; ///< TimesSeries containing points of time and the corresponding number of individuals transitioning
    // from one #InfectionState to another as defined in #Infection%s. In this case we use the update formula version.
    TimeSeries<ScalarType> populations; ///< TimeSeries containing points of time and the corresponding number of people
    // in defined #InfectionState%s. In this case we compute them by a sum.
    TimeSeries<ScalarType>
        populations_update; ///< TimeSeries containing points of time and the corresponding number of people
    // in defined #InfectionState%s. We compute them by an update formula.

private:
    // ---- Functionality for the iterations of a simulation.

    /**
    * @brief Computes number of Susceptibles for the current last time in populations.
    *
    * Number is computed by using the previous number of Susceptibles, total population of the last time point and 
    * the force of infection (also from the last time point).
    * Number is stored at the matching index in populations and populations_update.
    * 
    * @param[in] dt Time discretization step size.
    */
    void compute_susceptibles(ScalarType dt);

    /**
     * @brief Computes sizeof a flow
     * 
     * Computes size of one Transition from #InfectionTransition, specified in idx_InfectionTransitions, for the time index
     * current_time_index.
     *
     * @param[in] idx_InfectionTransitions Specifies the considered transition from #InfectionTransition
     * @param[in] idx_IncomingFlow Index of the transition in #InfectionTransition, which goes to the considered starting
     *      compartment of the transition specified in idx_InfectionTransitions. Size of considered flow is calculated via 
     *      the value of this incoming flow.
     * @param[in] idx_CurrentCompartment Index of the Compartment we flow out.
     * @param[in] current_time_index The time index the transition should be computed for.
     * @param[in] dt Time discretization step size.
     */
    void compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow,
                      Eigen::Index idx_CurrentCompartment, Eigen::Index current_time_index, ScalarType dt);

    /**
     * @brief Computes size of a flow
     * 
     * Computes size of one Transition from #InfectionTransition, specified in idx_InfectionTransitions, for the current
     * last time value in transitions.
     *
     * @param[in] idx_InfectionTransitions Specifies the considered transition from #InfectionTransition
     * @param[in] idx_IncomingFlow Index of the transition in #InfectionTransition, which goes to the considered starting
     *      compartment of the transition specified in idx_InfectionTransitions. Size of considered flow is calculated via 
     *      the value of this incoming flow.
     * @param[in] idx_CurrentCompartment Index of the Compartment we flow out.
     * @param[in] dt Time discretization step size.
     */
    void compute_flow(Eigen::Index idx_InfectionTransitions, Eigen::Index idx_IncomingFlow,
                      Eigen::Index idx_CurrentCompartment, ScalarType dt);

    /**
     * @brief Sets all required transitions for the current last timestep in transitions.
     * 
     * New values are stored in transitions. Most values are computed via the function compute_flow()
     *
     * @param[in] dt time discretization step size.
     */
    void flows_currents_timestep(ScalarType dt);

    /**
     * @brief Updates the values of one compartment, specified in infectionState, using all past transitions.
     * 
     * New value is stored in populations. The values are calculated using all past values for the incoming flows 
     * including the current time step. 
     * Therefore the flows of the current time step should be calculated before using this function.
     * @param[in] infectionState Specifies the #InfectionState we want to update.
     * @param[in] IncomingFlows 
     * @param[in] NaturalDeathispossible Boolian that determines if there is the possibility of Natural Death in infectionState.
     * @param[in] Transitionispossible Boolian that determines if it is possible to transition from the current InfectionState
     * into another
     * @param[in] dt
     */
    void update_compartment_with_sum(InfectionState infectionState,
                                     std::vector<InfectionTransition> const& IncomingFlows, bool NaturalDeathispossible,
                                     bool Transitionispossible, ScalarType dt);
    /**
     * @brief Updates the values of one compartment, specified in infectionState, using the transitions.
     * 
     * New value is stored in populations_update. The values are calculated using the compartment size in the previous 
     * time step and the related flows of the current time step. 
     * Therefore the flows of the current time step should be calculated before using this function.
     * @param[in] infectionState Specifies the #InfectionState we want to update.
     * @param[in] IncomingFlows 
     * @param[in] OutgoingFlows
     * @param[in] NaturalDeathispossible Boolian that determines if there is the possibility of Natural Death in infectionState.
     */
    void update_compartment_from_flow(InfectionState infectionState,
                                      std::vector<InfectionTransition> const& IncomingFlows,
                                      std::vector<InfectionTransition> const& OutgoingFlowsm,
                                      bool NaturalDeathispossible, ScalarType dt);

    /** 
     * @brief Updates the values of all compartments except Susceptible  
     *
     * New value is stored in populations. Values are computed via the function update_compartment_from_flow
     */
    void update_compartments(ScalarType dt);

    /**
     * @brief Compute the population size for the current last time in populations.
     * 
     * The population size is computed as the sum of all compartments. 
     * The population size is stored in the vector m_populationsize.
     */
    void compute_populationsize();

    /**
     * @brief Compute the normalized compartments for the current last time in m_normalizedpopulations.
     * 
     * The normalized compartments are computed as populations / m_populationsize.
     */
    void compute_normalizedcompartments();

    /**
     * @brief Computes the force of infection for the current last time in transitions.
     *
     * Computed value is stored in m_forceofinfection.
     *
     * @param[in] dt Time discretization step size.
     */
    void compute_forceofinfection(ScalarType dt);

    // ---- Private parameters. ----

    TimeSeries<ScalarType> m_forceofinfection{
        TimeSeries<ScalarType>(1)}; ///< TimeSeries containing the Force of infection term for every time point,
    // needed for the numerical scheme.
    TimeSeries<ScalarType> m_forceofinfectionupdate{
        TimeSeries<ScalarType>(1)}; ///< TimeSeries containing the Force of infection term for every time point,
    // needed for the numerical scheme. For the numerical scheme using the update formula for the compartments.
    TimeSeries<ScalarType> m_totalpopulation{TimeSeries<ScalarType>(
        1)}; ///< TimeSeries containing the total population size of the considered region for each time point.
    TimeSeries<ScalarType> m_totalpopulationupdate{TimeSeries<ScalarType>(
        1)}; ///< TimeSeries containing the total population size of the considered region for each time point.
    //In this case we use the compartments from populations_update.
    TimeSeries<ScalarType> m_totalpopulation_derivative{TimeSeries<ScalarType>(
        1)}; ///< TimeSeries containing the derivative of the total population size of the considered
    // region for each time point.
    TimeSeries<ScalarType> m_normalizedpopulations{
        TimeSeries<ScalarType>(Eigen::Index(InfectionState::Count) -
                               1)}; ///< TimeSeries containing points of time and the corresponding portion
    // of people in defined #IndectionState%s.

    // ---- Friend classes/functions. ----
    // In the Simulation class, the actual simulation is performed which is why it needs access to the here
    // defined (and private) functions to solve the model equations.
    friend class Simulation;
};

} // namespace endisecir
} // namespace mio

#endif //IDE_END_SECIR_MODEL_H