#ifndef IDE_END_SECIR_COMPPARAMS_H
#define IDE_END_SECIR_COMPPARAMS_H

#include "ide_endemic_secir/infection_state.h"
#include "ide_endemic_secir/parameters.h"
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
namespace endisecir
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
        //Exposed state:
        Eigen::Index support_max_index = (Eigen::Index)std::ceil(
            m_transitiondistributions_support_max[(int)InfectionTransition::ExposedToInfectedNoSymptoms] / dt);

        std::vector<ScalarType> vec_contribution_to_foi_1(support_max_index + 1, 0.);
        for (Eigen::Index i = 0; i <= support_max_index; i++) {
            ///Compute the state_age.
            ScalarType state_age = (ScalarType)i * dt;
            vec_contribution_to_foi_1[i] =
                parameters.get<TransitionDistributions>()[(int)InfectionTransition::ExposedToInfectedNoSymptoms].eval(
                    state_age);
        }
        m_transitiondistributions[(int)InfectionState::Exposed] = vec_contribution_to_foi_1;

        // Vector containing the transitions for the states with two outgoing flows.
        // We will compute the Transition Distribution for InfectedNoSymptoms and InfectedSymptoms until m_caltime, as we need
        // that many values to compute the mean_infectivity.

        std::vector<std::vector<int>> vector_transitions = {
            {(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms,
             (int)InfectionTransition::InfectedNoSymptomsToRecovered},
            {(int)InfectionTransition::InfectedSymptomsToInfectedSevere,
             (int)InfectionTransition::InfectedSymptomsToRecovered},
            {(int)InfectionTransition::InfectedSevereToInfectedCritical,
             (int)InfectionTransition::InfectedSevereToRecovered},
            {(int)InfectionTransition::InfectedCriticalToDead, (int)InfectionTransition::InfectedCriticalToRecovered}};

        for (int state = 2; state < (int)InfectionState::Count - 2; state++) {
            Eigen::Index calc_time_index = std::ceil(
                std::max(m_transitiondistributions_support_max[(Eigen::Index)vector_transitions[state - 2][0]],
                         m_transitiondistributions_support_max[(Eigen::Index)vector_transitions[state - 2][1]]) /
                dt);

            std::vector<ScalarType> vec_contribution_to_foi_2(calc_time_index + 1, 0.);
            for (Eigen::Index i = 0; i <= calc_time_index; i++) {
                ///Compute the state_age.
                ScalarType state_age = (ScalarType)i * dt;

                vec_contribution_to_foi_2[i] +=
                    parameters.get<TransitionProbabilities>()[vector_transitions[state - 2][0]] *
                        parameters.get<TransitionDistributions>()[vector_transitions[state - 2][0]].eval(state_age) +
                    parameters.get<TransitionProbabilities>()[vector_transitions[state - 2][1]] *
                        parameters.get<TransitionDistributions>()[vector_transitions[state - 2][1]].eval(state_age);
            }
            m_transitiondistributions[state] = vec_contribution_to_foi_2;
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
     *@brief Setter for the vector m_B that contains the approximated value of B for all for all necessary time points.
     * 
     * The value m_B is needed for the compuation of the force of infection term and also for the compuation of 
     * m_infectivity.
     *
     * @param[in] dt Time step size.
     */
    void set_B(ScalarType dt)
    {
        // Determine the calc_time_index that is the maximum of the support max of the TransitionDistributions that are
        // used in this computation.
        ScalarType calc_time =
            std::max(m_transitiondistributions_support_max[(int)InfectionTransition::InfectedSymptomsToInfectedSevere],
                     m_transitiondistributions_support_max[(int)InfectionTransition::InfectedSymptomsToRecovered]);
        Eigen::Index calc_time_index = (Eigen::Index)std::ceil(calc_time / dt) - 1;

        m_B = std::vector<ScalarType>(calc_time_index + 1, 0.);

        // We start the for loop for time_point_index = 1, as the next for loop where we compute the sum, goes up to
        // time_point_index-1. Therefore, for time_point_index = 0, the for loop would start at 0 and go up to -1, which
        // is not possible. The value m_B(0) is just set 0 zero.
        for (Eigen::Index time_point_index = 1; time_point_index <= calc_time_index; time_point_index++) {

            ScalarType sum                 = 0;
            Eigen::Index max_support_index = std::ceil(
                m_transitiondistributions_support_max[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] /
                dt);
            Eigen::Index starting_point = std::max(0, (int)time_point_index - (int)max_support_index);

            for (int i = starting_point; i <= time_point_index - 1; i++) {
                ScalarType state_age_i = static_cast<ScalarType>(i) * dt;

                sum +=
                    parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age_i) *
                    m_transitiondistributions[static_cast<int>(InfectionState::InfectedSymptoms)][i] *
                    m_transitiondistributions_derivative[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
                                                        [time_point_index - i];
            }
            m_B[time_point_index] =
                parameters
                    .get<TransitionProbabilities>()[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
                dt * sum;
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
        ScalarType calc_time_1 = std::max(
            m_transitiondistributions_support_max[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms],
            m_transitiondistributions_support_max[(int)InfectionTransition::InfectedNoSymptomsToRecovered]);
        Eigen::Index calc_time_index_1 = (Eigen::Index)std::ceil(calc_time_1 / dt) - 1;

        // Compute the calc_time_indedx corresponding to a_2:
        Eigen::Index calc_time_index_2 = m_B.size() - 1;

        // All in all calc_time_index for the infectivity:
        Eigen::Index calc_time_index = std::max(calc_time_index_1 + 1, calc_time_index_2 + 1);

        m_infectivity = std::vector<ScalarType>(calc_time_index, 0.);

        for (Eigen::Index time_point_index = 1; time_point_index < calc_time_index; time_point_index++) {
            ScalarType a_1 = 0;
            ScalarType a_2 = 0;

            ScalarType time_point_age = (ScalarType)time_point_index * dt;
            //We compute a_1 and a_2 at time_point.
            Eigen::Index max_support_index = std::ceil(
                m_transitiondistributions_support_max[(int)InfectionTransition::ExposedToInfectedNoSymptoms] / dt);
            Eigen::Index starting_point = std::max(0, (int)time_point_index - (int)max_support_index);

            //compute a_1:

            for (int i = starting_point; i < std::min(calc_time_index_1, time_point_index); i++) {

                ScalarType state_age_i = static_cast<ScalarType>(i) * dt;
                a_1 += parameters.get<RelativeTransmissionNoSymptoms>().eval(state_age_i) *
                       m_transitiondistributions[static_cast<int>(InfectionState::InfectedNoSymptoms)][i] *
                       m_transitiondistributions_derivative[(int)InfectionTransition::ExposedToInfectedNoSymptoms]
                                                           [time_point_index - i];
            }
            //compute a_2:

            for (int i = starting_point; i < std::min(calc_time_index_2, time_point_index); i++) {

                a_2 += m_transitiondistributions_derivative[(int)InfectionTransition::ExposedToInfectedNoSymptoms]
                                                           [time_point_index - i] *
                       m_B[i];
            }
            m_infectivity[time_point_index] =
                dt * std::exp(-parameters.get<NaturalDeathRate>() * time_point_age) * (a_2 - a_1);
        }
    }

    /**
     *@brief Setter for the vectors m_FoI_0 and m_NormFoI_0 that contain the approximated values of the function FoI_0,
     * that is a part of the force of infection term for all necessary time points.
     */
    void set_FoI_0(ScalarType dt)
    {
        Eigen::Index calc_time_index = m_transitiondistributions[(int)InfectionState::InfectedNoSymptoms].size();
        m_FoI_0                      = std::vector<ScalarType>(calc_time_index, 0.);
        m_NormFoI_0                  = std::vector<ScalarType>(calc_time_index, 0.);
        for (Eigen::Index time_point_index = 0; time_point_index < calc_time_index; time_point_index++) {
            ScalarType time_point_age = (ScalarType)time_point_index * dt;
            m_FoI_0[time_point_index] =
                (parameters.get<TransmissionProbabilityOnContact>().eval(time_point_age) *
                 parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(time_point_age)(0, 0)) *
                (m_statesinit[0][(int)InfectionState::InfectedNoSymptoms] *
                     m_transitiondistributions[(int)InfectionState::InfectedNoSymptoms][time_point_index] *
                     std::exp(-parameters.get<NaturalDeathRate>() * time_point_age) *
                     parameters.get<RelativeTransmissionNoSymptoms>().eval(time_point_age) +
                 m_statesinit[0][(int)InfectionState::InfectedSymptoms] *
                     m_transitiondistributions[(int)InfectionState::InfectedSymptoms][time_point_index] *
                     std::exp(-parameters.get<NaturalDeathRate>() * time_point_age) *
                     parameters.get<RiskOfInfectionFromSymptomatic>().eval(time_point_age));
            m_NormFoI_0[time_point_index] =
                (parameters.get<TransmissionProbabilityOnContact>().eval(time_point_age) *
                 parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(time_point_age)(0, 0)) *
                ((m_statesinit[0][(int)InfectionState::InfectedNoSymptoms] *
                      m_transitiondistributions[(int)InfectionState::InfectedNoSymptoms][time_point_index] *
                      std::exp(-parameters.get<NaturalDeathRate>() * time_point_age) *
                      parameters.get<RelativeTransmissionNoSymptoms>().eval(time_point_age) +
                  m_statesinit[0][(int)InfectionState::InfectedSymptoms] *
                      m_transitiondistributions[(int)InfectionState::InfectedSymptoms][time_point_index] *
                      std::exp(-parameters.get<NaturalDeathRate>() * time_point_age) *
                      parameters.get<RiskOfInfectionFromSymptomatic>().eval(time_point_age)) /
                 m_totalpopulationinit);
        }
    }

    /**
     * @brief Setter for the vectors m_InitFoI and m_NormInitFoI that contain the approximated values of the functions 
     * f(standard model) and g(normalized model) that is a part of the force of infection term for all necessary time 
     * points.
     *
     * @param[in] dt Time step size.
     */
    void set_InitFoI(ScalarType dt)
    {
        Eigen::Index calc_time_index = std::max(m_infectivity.size(), m_B.size());
        m_InitFoI.resize(calc_time_index, 0.);
        m_NormInitFoI.resize(calc_time_index, 0.);
        for (Eigen::Index time_point_index = 0; time_point_index < calc_time_index; time_point_index++) {
            ScalarType sum            = 0;
            ScalarType sum_norm       = 0;
            ScalarType time_point_age = (ScalarType)time_point_index * dt;
            if (time_point_index < (int)m_B.size()) {
                sum -= m_statesinit[0][(int)InfectionState::InfectedNoSymptoms] *
                       std::exp(-parameters.get<NaturalBirthRate>() * time_point_age) * m_B[time_point_index];
                sum_norm -= m_statesinit[0][(int)InfectionState::InfectedNoSymptoms] / m_totalpopulationinit *
                            std::exp(-parameters.get<NaturalBirthRate>() * time_point_age) * m_B[time_point_index];
            }
            if (time_point_index < (int)m_infectivity.size()) {
                sum += m_statesinit[0][(int)InfectionState::Exposed] * m_infectivity[time_point_index];
                sum_norm += m_statesinit[0][(int)InfectionState::Exposed] / m_totalpopulationinit *
                            m_infectivity[time_point_index];
            }
            m_InitFoI[time_point_index] =
                parameters.get<TransmissionProbabilityOnContact>().eval(time_point_age) *
                parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(time_point_age)(0, 0) * sum;
            m_NormInitFoI[time_point_index] =
                parameters.get<TransmissionProbabilityOnContact>().eval(time_point_age) *
                parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(time_point_age)(0, 0) * sum_norm;
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
        // The value T_z1^z2 is not defined for the transition SusceptiblesToExposed.
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
     * @brief Setter for m_V.
     *
    * @param[in] dt step size.
     */
    void set_V()
    {
        // The value V^z is not defined for the compartments Susceptibles and Exposed
        // InfectedNoSymptoms:
        m_V[(int)InfectionState::InfectedNoSymptoms] = m_T[(int)InfectionTransition::ExposedToInfectedNoSymptoms];
        // InfectedSymptoms:
        m_V[(int)InfectionState::InfectedSymptoms] =
            m_T[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
            m_V[(int)InfectionState::InfectedNoSymptoms];
        // InfectedSevere:
        m_V[(int)InfectionState::InfectedSevere] = m_T[(int)InfectionTransition::InfectedSymptomsToInfectedSevere] *
                                                   m_V[(int)InfectionState::InfectedSymptoms];
        // InfectedCritical:
        m_V[(int)InfectionState::InfectedCritical] =
            m_T[(int)InfectionTransition::InfectedSevereToInfectedCritical] * m_V[(int)InfectionState::InfectedSevere];
        // Dead:
        m_V[(int)InfectionState::Dead] =
            m_T[(int)InfectionTransition::InfectedCriticalToDead] * m_V[(int)InfectionState::InfectedCritical];
        // Recovered:
        m_V[(int)InfectionState::Recovered] =
            m_T[(int)InfectionTransition::InfectedNoSymptomsToRecovered] *
                m_V[(int)InfectionState::InfectedNoSymptoms] +
            m_T[(int)InfectionTransition::InfectedSymptomsToRecovered] * m_V[(int)InfectionState::InfectedSymptoms] +
            m_T[(int)InfectionTransition::InfectedSevereToRecovered] * m_V[(int)InfectionState::InfectedSevere] +
            m_T[(int)InfectionTransition::InfectedCriticalToRecovered] * m_V[(int)InfectionState::InfectedCritical];
    }

    /**
     * @brief Setter for m_W.
     *
     * m_W contains the values W_z for compartments z. As W_z is only defined for the compartments Exposed, InfectedNoSymptoms,
     * InfectedSymptoms, InfectedSevere and InfectedCrititcal we will set W_z for the other compartments to zero.
     * 
     * @param[in] dt step size.
     */
    void set_W(ScalarType dt)
    {
        for (int state = 1; state < (int)InfectionState::Count - 2; state++) {
            Eigen::Index calc_time_index = m_transitiondistributions[state].size();
            ScalarType sum               = 0;

            for (Eigen::Index i = 0; i < calc_time_index; i++) {
                // Compute state_age.
                ScalarType state_age = (ScalarType)i * dt;
                //Compute the apprximate derivative by a backwards difference scheme.
                sum += std::exp(-parameters.get<NaturalDeathRate>() * state_age) * m_transitiondistributions[state][i];
            }
            m_W[state] = sum;
        }
    }

    // ---- Private parameters. ----
    TimeSeries<ScalarType> m_statesinit; ///< TimeSeries containing the initial values for the compartments.
    ScalarType m_totalpopulationinit;
    ScalarType m_tol{1e-10}; ///< Tolerance used to calculate the maximum support of the TransitionDistributions.
    std::vector<ScalarType> m_transitiondistributions_support_max{
        std::vector<ScalarType>((int)InfectionTransition::Count, 0.)}; ///< A vector containing the support_max
    // for all TransitionDistributions.
    std::vector<std::vector<ScalarType>> m_transitiondistributions{std::vector<std::vector<ScalarType>>(
        (int)InfectionState::Count,
        std::vector<ScalarType>(1, 0.))}; ///< A vector containing the weighted TransitionDistributions.
    std::vector<std::vector<ScalarType>> m_transitiondistributions_derivative{std::vector<std::vector<ScalarType>>(
        (int)InfectionTransition::Count, std::vector<ScalarType>(1, 0.))}; ///< A Vector containing
    // the approximated derivative for all TransitionDistributions for all necessary time points.
    std::vector<ScalarType> m_B{std::vector<ScalarType>(1, 0.)}; ///< A Vector contaiing the appriximated
    // values for B.
    std::vector<ScalarType> m_infectivity{
        std::vector<ScalarType>(1, 0.)}; ///< A vector containing the approximated mean infectivity for all time points.
    ScalarType m_reproductionnumber_c; ///< The control Reproduction number
    std::vector<ScalarType> m_FoI_0{std::vector<ScalarType>(1, 0.)}; ///< A vector containing the approcimated
    // value of the function FoI_0 used for the computation of the force of infection in the standard model
    std::vector<ScalarType> m_NormFoI_0{std::vector<ScalarType>(1, 0.)}; ///< A vector containing the approcimated
    // value of the function FoI_0 used for the computation of the force of infection in the normalized model
    std::vector<ScalarType> m_InitFoI = std::vector<ScalarType>(1, 0.); ///< A vector containing the approcimated
    // value of the function FoI_0 used for the computation of the force of infection in the standard model
    std::vector<ScalarType> m_NormInitFoI{std::vector<ScalarType>(1, 0.)}; ///< A vector containing the approcimated
    // value of the function FoI_0 used for the computation of the force of infection in the normalized model
    std::vector<ScalarType> m_T{std::vector<ScalarType>((int)InfectionTransition::Count, 0.)}; ///< A vector
    // containing the approximated value for T_z1^z2 for every Flow z1 to z2.
    std::vector<ScalarType> m_V{std::vector<ScalarType>((int)InfectionTransition::Count, 0.)}; ///< A vector
    // containing the approximated value for V^z for every compartment z.
    std::vector<ScalarType> m_W{std::vector<ScalarType>((int)InfectionState::Count, 0.)}; ///< A vector containing+
    // the approximated value for W_z for every compartment z.
    // ---- Friend classes/functions. ----
    friend class Model;
    friend class NormModel;
    friend class Simulation;
}; // namespace endisecir

} // namespace endisecir

} // namespace mio

#endif //IDE_END_SECIR_COMPPARAMS_H