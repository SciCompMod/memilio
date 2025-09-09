#ifndef IDE_END_SIRD_COMPPARAMS_H
#define IDE_END_SIRD_COMPPARAMS_H

#include "ide_endemic_sird/infection_state.h"
#include "ide_endemic_sird/parameters.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"

#include "vector"
#include <Eigen/src/Core/util/Meta.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace mio
{
namespace endisird
{
class CompParameters;
class Model;
class Simulation;

class CompParameters
{
    using ParameterSet = Parameters;

public:
    CompParameters(TimeSeries<ScalarType>&& states_init)
        : parameters{Parameters()}
        , m_statesinit{std::move(states_init)}
    {
        m_totalpopulationinit = std::accumulate(m_statesinit[0].begin(), m_statesinit[0].end(), 0) -
                                m_statesinit[0][(int)InfectionState::Dead];
    }

    bool check_constraints() const
    {
        if (!(static_cast<size_t>(m_statesinit.get_num_elements()) == static_cast<size_t>(InfectionState::Count))) {
            log_error(
                " A variable given for model construction is not valid. Number of elements in vector of populations "
                "does not match the required number.");
            return true;
        }
        for (int i = 0; i < static_cast<int>(InfectionState::Count); i++) {
            if (m_statesinit[0][i] < 0) {
                log_error("Initialization failed. Initial values for populations are less than zero.");
                return true;
            }
        }
        return parameters.check_constraints();
    }

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

private:
    // ---- Functionality to set vectors with necessary information regarding TransitionDistributions. ----
    /**
     * @brief Setter for the vector m_transitiondistributions_support_max that contains the support_max for all 
     * TransitionDistributions.
     *
     * This determines how many summands are required when calculating flows, the force of infection or compartments.
     *
     * @param[in] dt Time step size.
     */
    void set_transitiondistributions_support_max(ScalarType dt)
    {
        // The transition Susceptible to Exposed is not needed in the computations.
        for (int transition = 1; transition < (int)InfectionTransition::Count; transition++) {
            m_transitiondistributions_support_max[transition] =
                parameters.get<TransitionDistributions>()[transition].get_support_max(dt, m_tol);
        }
    }

    /**
     * @brief Setter for the vector m_transitiondistributions.
     *
     * Here we compute the weighted transition distributions for every compartment. Meaning for every compartment we weight
     * the distributions of the outgoing transitions with their probability and compute the sum of the two weighed distributions.
     * 
     * @param[in] dt Time step size.
     */
    void set_transitiondistributions(ScalarType dt)
    {

        // Vector containing the outgoing transitions for the Infected state.

        std::vector<int> vector_transitions = {(int)InfectionTransition::InfectedToDead,
                                               (int)InfectionTransition::InfectedToRecovered};

        Eigen::Index calc_time_index =
            std::ceil(std::max(m_transitiondistributions_support_max[(Eigen::Index)vector_transitions[0]],
                               m_transitiondistributions_support_max[(Eigen::Index)vector_transitions[1]]) /
                      dt);

        std::vector<ScalarType> vec_contribution_to_foi_2(calc_time_index + 1, 0.);
        for (Eigen::Index i = 0; i <= calc_time_index; i++) {
            ///Compute the state_age.
            ScalarType state_age = (ScalarType)i * dt;

            vec_contribution_to_foi_2[i] +=
                parameters.get<TransitionProbabilities>()[vector_transitions[0]] *
                    parameters.get<TransitionDistributions>()[vector_transitions[0]].eval(state_age) +
                parameters.get<TransitionProbabilities>()[vector_transitions[1]] *
                    parameters.get<TransitionDistributions>()[vector_transitions[1]].eval(state_age);
            m_transitiondistributions = vec_contribution_to_foi_2;
        }
    }

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
    void set_transitiondistributions_derivative(ScalarType dt)
    {
        // The transition SusceptibleToExposed is not needed in the computations.
        for (int transition = 1; transition < (int)InfectionTransition::Count; transition++) {
            Eigen::Index support_max_index =
                (Eigen::Index)std::ceil(m_transitiondistributions_support_max[transition] / dt);

            // Create vec_derivative that stores the value of the approximated derivative for all necessary time points.
            std::vector<ScalarType> vec_derivative(support_max_index + 1, 0.);

            for (Eigen::Index i = 0; i <= support_max_index; i++) {
                // Compute state_age.
                ScalarType state_age = (ScalarType)i * dt;
                //Compute the apprximate derivative by a backwards difference scheme.
                vec_derivative[i] = (parameters.get<TransitionDistributions>()[transition].eval(state_age) -
                                     parameters.get<TransitionDistributions>()[transition].eval(state_age - dt)) /
                                    dt;
            }
            m_transitiondistributions_derivative[transition] = vec_derivative;
        }
    }

    /**
     *@brief Setter for the vector m_meaninfectivity that contains the approximated value of the mean infectivity for all
     * for all necessary time points.
     * 
     * This values is needed the compute the reproduction numer and the force of infection term.
     *
     * @param[in] dt Time step size.
     */
    void set_infectivity(ScalarType dt)
    {
        // Compute the calc_time_indedx corresponding to a_1:
        ScalarType calc_time_index = m_transitiondistributions.size();

        m_infectivity = std::vector<ScalarType>(calc_time_index, 0.);

        for (Eigen::Index time_point_index = 1; time_point_index < calc_time_index; time_point_index++) {
            ScalarType time_point_age       = (ScalarType)time_point_index * dt;
            m_infectivity[time_point_index] = parameters.get<RiskOfInfectionFromSymptomatic>().eval(time_point_age) *
                                              m_transitiondistributions[time_point_index] *
                                              std::exp(-parameters.get<NaturalDeathRate>() * time_point_age);
        }
    }

    /**
     *@brief Setter for the vectors m_FoI_0 and m_NormFoI_0 that contain the approximated values of the function FoI_0,
     * that is a part of the force of infection term for all necessary time points.
     */
    void set_FoI_0(ScalarType dt)
    {
        Eigen::Index calc_time_index = m_transitiondistributions.size();
        m_FoI_0                      = std::vector<ScalarType>(calc_time_index, 0.);
        m_NormFoI_0                  = std::vector<ScalarType>(calc_time_index, 0.);
        for (Eigen::Index time_point_index = 0; time_point_index < calc_time_index; time_point_index++) {
            ScalarType time_point_age = (ScalarType)time_point_index * dt;
            m_FoI_0[time_point_index] =
                (parameters.get<TransmissionProbabilityOnContact>().eval(time_point_age) *
                 parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(time_point_age)(0, 0)) *
                m_statesinit[0][(int)InfectionState::Infected] * m_infectivity[time_point_index];
            m_NormFoI_0[time_point_index] =
                (parameters.get<TransmissionProbabilityOnContact>().eval(time_point_age) *
                 parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(time_point_age)(0, 0)) *
                (m_statesinit[0][(int)InfectionState::Infected] / m_totalpopulationinit) *
                m_infectivity[time_point_index];
        }
    }

    // ---- Parameters needed for the analysis of the model ----
    /**
     * @brief Setter for the Reproduction number m_reproductionnumber_c.
     *
     */
    void set_reproductionnumber_c(ScalarType dt)
    {
        // Determine the corresponding time index.

        Eigen::Index calc_time_index = m_infectivity.size() - 1;
        ScalarType sum               = 0;
        for (int i = 0; i <= calc_time_index; i++) {
            sum += m_infectivity[i];
        }
        m_reproductionnumber_c = parameters.get<TransmissionProbabilityOnContact>().eval(0) *
                                 parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(0)(0, 0) * dt *
                                 sum;
    }

    /**
     * @brief Setter for m_T.
     *
    * @param[in] dt step size.
     */
    void set_T(ScalarType dt)
    {
        // The value T_z1^z2 is not defined for the transition SusceptiblesToInfected.
        for (int transition = 1; transition < (int)InfectionTransition::Count; transition++) {
            Eigen::Index support_max_index =
                (Eigen::Index)std::ceil(m_transitiondistributions_support_max[transition] / dt);
            ScalarType sum = 0;

            for (Eigen::Index i = 0; i <= support_max_index; i++) {
                // Compute state_age.
                ScalarType state_age = (ScalarType)i * dt;
                //Compute the apprximate derivative by a backwards difference scheme.
                sum -= std::exp(-parameters.get<NaturalDeathRate>() * state_age) *
                       m_transitiondistributions_derivative[transition][i];
            }
            m_T[transition] = parameters.get<TransitionProbabilities>()[transition] * sum;
        }
    }

    /**
     * @brief Setter for m_W.
     *
     * m_W contains the values W_i 
     *
     * @param[in] dt step size.
     */
    void set_W(ScalarType dt)
    {

        Eigen::Index calc_time_index = m_transitiondistributions.size();
        ScalarType sum               = 0;

        for (Eigen::Index i = 0; i < calc_time_index; i++) {
            // Compute state_age.
            ScalarType state_age = (ScalarType)i * dt;
            //Compute the apprximate derivative by a backwards difference scheme.
            sum += std::exp(-parameters.get<NaturalDeathRate>() * state_age) * m_transitiondistributions[i];
        }
        m_W = sum;
    }

    void set_equilibria()
    {
        // We start by setting the equilibrium for the infected compartments.
        ScalarType p = (m_W / m_T[(int)InfectionTransition::InfectedToDead]) *
                           (parameters.get<NaturalDeathRate>() - 1 / m_W - parameters.get<NaturalBirthRate>()) -
                       parameters.get<NaturalBirthRate>() /
                           (m_reproductionnumber_c - m_T[(int)InfectionTransition::InfectedToDead]);
        ScalarType q = parameters.get<NaturalBirthRate>() * m_W /
                       ((m_reproductionnumber_c - m_T[(int)InfectionTransition::InfectedToDead]) *
                        m_T[(int)InfectionTransition::InfectedToDead]) *
                       ((1 / m_W) - parameters.get<NaturalDeathRate>() + parameters.get<NaturalBirthRate>() -
                        m_reproductionnumber_c);
        m_compartments_equilibrium[0][(int)InfectionState::Infected] = -p / 2 + std::sqrt((p / 2) * (p / 2) - q);
        m_compartments_equilibrium[1][(int)InfectionState::Infected] = -p / 2 - std::sqrt((p / 2) * (p / 2) - q);

        // From this we set the other equilibria points.
        // Susceptibles:
        m_compartments_equilibrium[0][(int)InfectionState::Susceptible] =
            parameters.get<NaturalBirthRate>() /
            (m_compartments_equilibrium[0][(int)InfectionState::Infected] *
                 (m_reproductionnumber_c / m_W - m_T[(int)InfectionTransition::InfectedToDead] / m_W) +
             parameters.get<NaturalBirthRate>());
        m_compartments_equilibrium[1][(int)InfectionState::Susceptible] =
            parameters.get<NaturalBirthRate>() /
            (m_compartments_equilibrium[1][(int)InfectionState::Infected] *
                 (m_reproductionnumber_c / m_W - m_T[(int)InfectionTransition::InfectedToDead] / m_W) +
             parameters.get<NaturalBirthRate>());
        // Force of Infection:
        m_FoI_equilibrium[0] =
            m_compartments_equilibrium[0][(int)InfectionState::Infected] * (m_reproductionnumber_c / m_W);
        m_FoI_equilibrium[1] =
            m_compartments_equilibrium[1][(int)InfectionState::Infected] * (m_reproductionnumber_c / m_W);
        // InfectedToDead:
        m_transitions_equilibrium[0][(int)InfectionTransition::InfectedToDead] =
            m_compartments_equilibrium[0][(int)InfectionState::Infected] *
            (m_T[(int)InfectionTransition::InfectedToDead] / m_W);
        m_transitions_equilibrium[1][(int)InfectionTransition::InfectedToDead] =
            m_compartments_equilibrium[1][(int)InfectionState::Infected] *
            (m_T[(int)InfectionTransition::InfectedToDead] / m_W);
        // InfectedToRecovered:
        m_transitions_equilibrium[0][(int)InfectionTransition::InfectedToRecovered] =
            (m_compartments_equilibrium[0][(int)InfectionState::Susceptible] * m_FoI_equilibrium[0] +
             (parameters.get<NaturalDeathRate>() +
              m_transitions_equilibrium[0][(int)InfectionTransition::InfectedToDead] -
              parameters.get<NaturalBirthRate>()) *
                 m_compartments_equilibrium[0][(int)InfectionState::Infected]) *
            m_T[(int)InfectionTransition::InfectedToDead];
        m_transitions_equilibrium[1][(int)InfectionTransition::InfectedToRecovered] =
            (m_compartments_equilibrium[1][(int)InfectionState::Susceptible] * m_FoI_equilibrium[0] +
             (parameters.get<NaturalDeathRate>() +
              m_transitions_equilibrium[1][(int)InfectionTransition::InfectedToDead] -
              parameters.get<NaturalBirthRate>()) *
                 m_compartments_equilibrium[1][(int)InfectionState::Infected]) *
            m_T[(int)InfectionTransition::InfectedToDead];
        // Recovered:
        m_compartments_equilibrium[0][(int)InfectionState::Recovered] =
            m_transitions_equilibrium[0][(int)InfectionTransition::InfectedToRecovered] /
            (parameters.get<NaturalBirthRate>() -
             m_transitions_equilibrium[0][(int)InfectionTransition::InfectedToDead]);
        m_compartments_equilibrium[1][(int)InfectionState::Recovered] =
            m_transitions_equilibrium[1][(int)InfectionTransition::InfectedToRecovered] /
            (parameters.get<NaturalBirthRate>() -
             m_transitions_equilibrium[1][(int)InfectionTransition::InfectedToDead]);
    }

    // ---- Private parameters. ----
    TimeSeries<ScalarType> m_statesinit; ///< TimeSeries containing the initial values for the compartments.
    ScalarType m_totalpopulationinit;
    ScalarType m_tol{1e-10}; ///< Tolerance used to calculate the maximum support of the TransitionDistributions.
    std::vector<ScalarType> m_transitiondistributions_support_max{
        std::vector<ScalarType>((int)InfectionTransition::Count, 0.)}; ///< A vector containing the support_max
    // for all TransitionDistributions.
    std::vector<ScalarType> m_transitiondistributions{
        std::vector<ScalarType>(1, 0.)}; ///< A vector containing the weighted TransitionDistributions
    // for the Infected state.
    std::vector<std::vector<ScalarType>> m_transitiondistributions_derivative{std::vector<std::vector<ScalarType>>(
        (int)InfectionTransition::Count, std::vector<ScalarType>(1, 0.))}; ///< A Vector containing
    // the approximated derivative for all TransitionDistributions for all necessary time points.
    std::vector<ScalarType> m_infectivity{
        std::vector<ScalarType>(1, 0.)}; ///< A vector containing the approximated mean infectivity for all time points.
    ScalarType m_reproductionnumber_c; ///< The control Reproduction number
    std::vector<ScalarType> m_FoI_0{std::vector<ScalarType>(1, 0.)}; ///< A vector containing the approcimated
    // value of the function FoI_0 used for the computation of the force of infection in the standard model
    std::vector<ScalarType> m_NormFoI_0{std::vector<ScalarType>(1, 0.)}; ///< A vector containing the approcimated
    // value of the function FoI_0 used for the computation of the force of infection in the normalized model
    std::vector<ScalarType> m_T{std::vector<ScalarType>((int)InfectionTransition::Count, 0.)}; ///< A vector
    // containing the approximated value for T_z1^z2 for every Flow z1 to z2.
    ScalarType m_W{0}; ///< ScalarType of the value W_i.
    std::vector<std::vector<ScalarType>> m_compartments_equilibrium{std::vector<std::vector<ScalarType>>(
        2, std::vector<ScalarType>((int)InfectionState::Count - 1, 0.))}; ///< Vector containing the
    // two computed equilibria points for the compartments of NormModel.
    std::vector<std::vector<ScalarType>> m_transitions_equilibrium{std::vector<std::vector<ScalarType>>(
        2, std::vector<ScalarType>((int)InfectionState::Count - 1, 0.))}; ///< Vector containing the
    // two computed equilibria points for the transitions of NormModel.
    std::vector<ScalarType> m_FoI_equilibrium{std::vector<ScalarType>(2, 0.)}; ///< A Vector containing the two
    // computed equilibria points for the force of infection of NormModel.
    // ---- Friend classes/functions. ----
    friend class Model;
    friend class NormModel;
    friend class Simulation;
};

} // namespace endisird

} // namespace mio

#endif //IDE_END_SIRD_COMPPARAMS_H