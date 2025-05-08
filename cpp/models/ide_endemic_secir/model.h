#ifndef IDE_END_SECIR_MODEL_H
#define IDE_END_SECIR_MODEL_H

#include "ide_endemic_secir/infection_state.h"
#include "ide_endemic_secir/parameters.h"
#include "memilio/config.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/time_series.h"

#include "vector"
#include <Eigen/src/Core/util/Meta.h>
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
    * @param[in] states_init A vector, containing the number of individuals in each compartment at time 0.
    */

    Model(TimeSeries<ScalarType>&& states_init);

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
    void set_tol_for_support_max(ScalarType new_tol)
    {
        m_tol = new_tol;
    }

    // ---- Public parameters. ----
    ParameterSet parameters; ///< ParameterSet of Model Parameters.
    TimeSeries<ScalarType>
        transitions; ///< TimesSeries containing points of time and the corresponding number of individuals transitioning
    // from one #InfectionState to another as defined in #Infection%s.
    TimeSeries<ScalarType> populations; ///< TimeSeries containing points of time and the corresponding number of people
    // in defined #InfectionState%s.

private:
    // ---- Functionality for the iterations of a simulation.

    /**
    * @brief Computes number of Susceptibles for the current last time in populations.
    *
    * Number is computed by using the previous number of Susceptibles, total population of the last time point and 
    * the force of infection (also from the last time point).
    * Number is stored at the matching index in populations.
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
     * @brief Updates the values of one compartment, specified in infectionState, using the transitions.
     * 
     * New value is stored in populations. The values are calculated using the compartment size in the previous 
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
     * @brief Compute the normalized compartments for the current last time in normalized_populations.
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

    // ---- Functionality to set vectors with necessary information regarding TransitionDistributions. ----
    /**
     * @brief Setter for the vector m_transitiondistributions_support_max that contains the support_max for all 
     * TransitionDistributions.
     *
     * This determines how many summands are required when calculating flows, the force of infection or compartments.
     *
     * @param[in] dt Time step size.
     */
    void set_transitiondistributions_support_max(ScalarType dt);

    /**
     * @brief Setter for the vector m_transitiondistributions.
     *
     * In the computation of the force of infection in the initialization function of the force of infection term,
     * we evaluate the survival functions of the TransitionDistributions InfectedNoSymptomsToInfectedSymptoms, 
     * InfectedNoSymptomsToRecovered, InfectedSymptomsToInfectedSevere and InfectedSymptomsToRecovered, weighted by 
     * the corresponding TransitionProbabilities, at the same time points.
     * Here, we compute these contributions to the force of infection term, an store the vector 
     * m_transitiondistributions_in_forceofinfection so that we can access this vector for all following computations.
     * 
     * @param[in] dt Time step size.
     */
    void set_transitiondistributions(ScalarType dt);

    void set_calctime();

    /**
     * @brief Setter for the vector m_transitiondistributions_derivative that contains the approximated derivative for 
     * all TransitionDistributions for all necessary time points.
     *
     * The derivative is approximated using a backwards difference scheme.
     * The number of necessary time points for each TransitionDistribution is determined using 
     * m_transitiondistributions_support_max.
     *
     * @param[in] dt Time step size.
     */
    void set_transitiondistributions_derivative(ScalarType dt);

    /**
     *@brief Setter for the vector m_meaninfectivity that contains the approximated value of the mean infectivity for all
     * for all necessary time points.
     * 
     * @param[in] dt Time step size.
     */
    void set_meaninfectivity(ScalarType dt);

    // // /**
    // //  *@brief Setter for the vector m_initalvaluesforceofinfection that contains the approximated value of the part of the
    // //  * force of infection term that is influenced by the initial values.
    // //  *
    // //  * @param[in] dt Time step size.
    // //  */
    // void set_initalvaluesforceofinfection(ScalarType dt);

    // ---- Functionality for the model analysis. ----

    /**
     * @brief Setter for the Reproduction number R_c.
     *
     */
    void set_reproductionnumber_c(ScalarType dt);

    //TODO Beschreibung
    void set_probability_of_transition(ScalarType dt);

    //TODO Beschreibung
    void set_meansojourntime(ScalarType dt);

    /**
     * @brief Setter for the equilibrium of the normalized model
     *
     * If the m_reproductionnumber_c < 1 we use the disease free equilibrium and if m_reproductionnumber_c > 1 we compute
     * the endemic equilibrium.
     *
    */
    void set_equilibrium();

    // ---- Private parameters. ----

    TimeSeries<ScalarType> m_forceofinfection{
        TimeSeries<ScalarType>(1)}; ///< TimeSeries containing the Force of infection term for every time point,
    // needed for the numerical scheme.
    TimeSeries<ScalarType> m_totalpopulation{TimeSeries<ScalarType>(
        1)}; ///< TimeSeries containing the total population size of the considered region for each time point.
    TimeSeries<ScalarType> m_normalizedpopulations{
        TimeSeries<ScalarType>(Eigen::Index(InfectionState::Count) -
                               1)}; ///< TimeSeries containing points of time and the corresponding portion
    // of people in defined #IndectionState%s.

    ScalarType m_reproductionnumber_c; ///< The control Reproduction number

    ScalarType m_tol{1e-10}; ///< Tolerance used to calculate the maximum support of the TransitionDistributions.
    ScalarType m_calctime{
        0}; ///< A ScalarType wit he calc time determined by the support max of the transition distributions.
    std::vector<ScalarType> m_transitiondistributions_support_max{
        std::vector<ScalarType>((int)InfectionTransition::Count, 0.)}; ///< A vector containing the support_max
    // for all TransitionDistributions.
    std::vector<std::vector<ScalarType>> m_transitiondistributions_in_forceofinfection{
        std::vector<std::vector<ScalarType>>(
            2, std::vector<ScalarType>(1, 0.))}; ///> A vectpr containing the weighted TransitionDistributions
    // needed  in the compuation of the initial functions for the compartments and the mean infectivity at all necessary
    // time points.
    std::vector<std::vector<ScalarType>> m_transitiondistributions_derivative{std::vector<std::vector<ScalarType>>(
        (int)InfectionTransition::Count, std::vector<ScalarType>(1, 0.))}; ///> A Vector containing
    // the approximated derivative for all TransitionDistributions for all necessary time points.
    std::vector<ScalarType>
        m_meaninfectivity; ///> a vector containing the approximated mean infectivity for all time points.

    std::vector<ScalarType> m_equilibriumnormalizedcompartments{
        std::vector<ScalarType>(int(InfectionState::Count), 0.)}; ///< Vector containing
    // the equilibrium for the normalized compartments.

    ScalarType m_equilibriumforceofinfection{0.}; ///< The equilibrium of the force of infection.

    std::vector<ScalarType> m_probabilityoftransition{
        std::vector<ScalarType>((int)InfectionTransition::Count, 0.)}; ///<TODO

    std::vector<ScalarType> m_meansojourntime{std::vector<ScalarType>((int)InfectionState::Count - 2, 0.)}; ///<TODO

    // std::vector<ScalarType> m_initalvaluesforceofinfection; ///> a vector containing the parts of the force of infection
    // // that are influenced by the inital values of our model, i.e phi_0 and f.

    // ---- Friend classes/functions. ----
    // In the Simulation class, the actual simulation is performed which is why it needs access to the here
    // defined (and private) functions to solve the model equations.
    friend class Simulation;
};

} // namespace endisecir
} // namespace mio

#endif //IDE_END_SECIR_MODEL_H