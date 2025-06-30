/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Lena Ploetzke
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#ifndef MIO_GLCT_SECIR_MODEL_H
#define MIO_GLCT_SECIR_MODEL_H

#include "glct_secir/parameters.h"
#include "glct_secir/infection_state.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/compartments/compartmental_model.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/config.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"

namespace mio
{
namespace glsecir
{

/**
 * @brief Class that defines an GLCT-SECIR model.
 *
 * @tparam NumExposed The number of subcompartments used for the Exposed compartment.
 * @tparam NumInfectedNoSymptoms The number of subcompartments used for the InfectedNoSymptoms compartment. 
 * @tparam NumInfectedSymptoms The number of subcompartments used for the InfectedSymptoms compartment.
 * @tparam NumInfectedSevere The number of subcompartments used for the InfectedSevere compartment.
 * @tparam NumInfectedCritical The number of subcompartments used for the InfectedCritical compartment.
 */
template <size_t NumExposed, size_t NumInfectedNoSymptoms, size_t NumInfectedSymptoms, size_t NumInfectedSevere,
          size_t NumInfectedCritical>
class Model
    : public CompartmentalModel<
          ScalarType,
          LctInfectionState<InfectionState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms,
                            NumInfectedSevere, NumInfectedCritical, 1, 1>,
          mio::Populations<ScalarType,
                           LctInfectionState<InfectionState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms,
                                             NumInfectedSevere, NumInfectedCritical, 1, 1>>,
          Parameters>
{

public:
    using LctState =
        LctInfectionState<InfectionState, 1, NumExposed, NumInfectedNoSymptoms, NumInfectedSymptoms, NumInfectedSevere,
                          NumInfectedCritical, 1, 1>; ///< This class specifies the number of subcompartments.
    using Base = CompartmentalModel<ScalarType, LctState, mio::Populations<ScalarType, LctState>, Parameters>;
    using typename Base::ParameterSet;
    using typename Base::Populations;

    /// @brief Default constructor.
    Model()
        : Base(Populations({Index<LctState>(LctState::Count)}, 0.), ParameterSet())
    {
    }

    /**
     * @brief Checks that the model satisfies all constraints (e.g. parameter or population constraints), and 
     *  logs an error if constraints are not satisfied.
     *
     * @return Returns true if one or more constraints are not satisfied, false otherwise.
     */
    bool check_constraints() const
    {
        auto params = this->parameters;

        // --- Check that the dimensions are consistent. ---
        if ((Eigen::Index)LctState::template get_num_subcompartments<InfectionState::Exposed>() !=
            params.template get<StartingProbabilitiesExposed>().rows()) {
            log_error("Constraint check: Dimension of the parameters does not match the number of subcompartments for "
                      "the Exposed "
                      "compartment.");
            return true;
        }
        if ((Eigen::Index)LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>() !=
            params.template get<StartingProbabilitiesInfectedNoSymptoms>().rows()) {
            log_error(
                "Constraint check: Dimension of the parameters does not match the number of subcompartments for the "
                "InfectedNoSymptoms compartment.");
            return true;
        }
        if ((Eigen::Index)LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>() !=
            params.template get<StartingProbabilitiesInfectedSymptoms>().rows()) {
            log_error(
                "Constraint check: Dimension of the parameters does not match the number of subcompartments for the "
                "InfectedSymptoms compartment.");
            return true;
        }
        if ((Eigen::Index)LctState::template get_num_subcompartments<InfectionState::InfectedSevere>() !=
            params.template get<StartingProbabilitiesInfectedSevere>().rows()) {
            log_error("Constraint check: Dimension of the parameters does not match the number of subcompartments for "
                      "the InfectedSevere "
                      "compartment.");
            return true;
        }
        if ((Eigen::Index)LctState::template get_num_subcompartments<InfectionState::InfectedCritical>() !=
            params.template get<StartingProbabilitiesInfectedCritical>().rows()) {
            log_error(
                "Constraint check: Dimension of the parameters does not match the number of subcompartments for the "
                "InfectedCritical compartment.");
            return true;
        }

        return (params.check_constraints() || this->populations.check_constraints());
    }

    /**
     * @brief Evaluates the right-hand-side f of the GLCT dydt = f(y, t).
     *
     * The GLCT-SECIR model is defined through ordinary differential equations of the form dydt = f(y, t). 
     * y is a vector containing number of individuals for each (sub-) compartment.
     * This function evaluates the right-hand-side f of the ODE and can be used in an ODE solver.
     *
     * @param[in] pop The current state of the population in the geographic unit we are considering.
     * @param[in] y The current state of the model (or a subpopulation) as a flat array.
     * @param[in] t The current time.
     * @param[out] dydt A reference to the calculated output.
     */
    void get_derivatives(Eigen::Ref<const Eigen::VectorX<ScalarType>> pop,
                         Eigen::Ref<const Eigen::VectorX<ScalarType>> y, ScalarType t,
                         Eigen::Ref<Eigen::VectorX<ScalarType>> dydt) const override
    {
        dydt.setZero();

        auto params           = this->parameters;
        auto total_population = pop.sum() - pop[LctState::template get_first_index<InfectionState::Dead>()];

        // Calculate sum of all subcompartments for InfectedNoSymptoms.
        ScalarType InfectedNoSymptoms_sum =
            pop.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                        LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>())
                .sum();
        // Calculate sum of all subcompartments for InfectedSymptoms.
        ScalarType InfectedSymptoms_sum =
            pop.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                        LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>())
                .sum();

        // --- Susceptibles. ---
        ScalarType season_val = 1 + params.template get<Seasonality>() *
                                        sin(3.141592653589793 * ((params.template get<StartDay>() + t) / 182.5 + 0.5));
        dydt[0] = -y[0] / total_population * season_val * params.template get<TransmissionProbabilityOnContact>() *
                  params.template get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(t)(0, 0) *
                  (params.template get<RelativeTransmissionNoSymptoms>() * InfectedNoSymptoms_sum +
                   params.template get<RiskOfInfectionFromSymptomatic>() * InfectedSymptoms_sum);

        // --- Exposed. ---
        dydt.segment(LctState::template get_first_index<InfectionState::Exposed>(),
                     LctState::template get_num_subcompartments<InfectionState::Exposed>()) -=
            dydt[0] * params.template get<StartingProbabilitiesExposed>();
        dydt.segment(LctState::template get_first_index<InfectionState::Exposed>(),
                     LctState::template get_num_subcompartments<InfectionState::Exposed>()) +=
            params.template get<TransitionMatrixExposedToInfectedNoSymptoms>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::Exposed>(),
                      LctState::template get_num_subcompartments<InfectionState::Exposed>());

        // --- InfectedNoSymptoms. ---
        // Flow from Exposed To InfectedNoSymptoms.
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedNoSymptoms>()) =
            -(params.template get<TransitionMatrixExposedToInfectedNoSymptoms>() *
              Eigen::VectorX<ScalarType>::Ones(LctState::template get_num_subcompartments<InfectionState::Exposed>()))
                 .transpose() *
            y.segment(LctState::template get_first_index<InfectionState::Exposed>(),
                      LctState::template get_num_subcompartments<InfectionState::Exposed>()) *
            params.template get<StartingProbabilitiesInfectedNoSymptoms>();
        // Flow from InfectedNoSymptoms To InfectedSymptoms.
        size_t dimensionInfectedNoSymptomsToInfectedSymptoms =
            params.template get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>().rows();
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                     dimensionInfectedNoSymptomsToInfectedSymptoms) +=
            params.template get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                      dimensionInfectedNoSymptomsToInfectedSymptoms);
        // Flow from InfectedNoSymptoms To Recovered.
        size_t dimensionInfectedNoSymptomsToRecovered =
            params.template get<TransitionMatrixInfectedNoSymptomsToRecovered>().rows();
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() +
                         dimensionInfectedNoSymptomsToInfectedSymptoms,
                     dimensionInfectedNoSymptomsToRecovered) +=
            params.template get<TransitionMatrixInfectedNoSymptomsToRecovered>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() +
                          dimensionInfectedNoSymptomsToInfectedSymptoms,
                      dimensionInfectedNoSymptomsToRecovered);
        // Add flow directly to Recovered compartment.
        dydt[LctState::template get_first_index<InfectionState::Recovered>()] +=
            -(params.template get<TransitionMatrixInfectedNoSymptomsToRecovered>() *
              Eigen::VectorX<ScalarType>::Ones(dimensionInfectedNoSymptomsToRecovered))
                 .transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>() +
                          dimensionInfectedNoSymptomsToInfectedSymptoms,
                      dimensionInfectedNoSymptomsToRecovered);

        // --- InfectedSymptoms. ---
        // Flow from InfectedNoSymptoms To InfectedSymptoms.
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedSymptoms>()) =
            -params.template get<StartingProbabilitiesInfectedSymptoms>() *
            (params.template get<TransitionMatrixInfectedNoSymptomsToInfectedSymptoms>() *
             Eigen::VectorX<ScalarType>::Ones(dimensionInfectedNoSymptomsToInfectedSymptoms))
                .transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedNoSymptoms>(),
                      dimensionInfectedNoSymptomsToInfectedSymptoms);
        // Flow from InfectedSymptoms To InfectedSevere.
        size_t dimensionInfectedSymptomsToInfectedSevere =
            params.template get<TransitionMatrixInfectedSymptomsToInfectedSevere>().rows();
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                     dimensionInfectedSymptomsToInfectedSevere) +=
            params.template get<TransitionMatrixInfectedSymptomsToInfectedSevere>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                      dimensionInfectedSymptomsToInfectedSevere);
        // Flow from InfectedSymptoms To Recovered.
        size_t dimensionInfectedSymptomsToRecovered =
            params.template get<TransitionMatrixInfectedSymptomsToRecovered>().rows();
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>() +
                         dimensionInfectedSymptomsToInfectedSevere,
                     dimensionInfectedSymptomsToRecovered) +=
            params.template get<TransitionMatrixInfectedSymptomsToRecovered>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>() +
                          dimensionInfectedSymptomsToInfectedSevere,
                      dimensionInfectedSymptomsToRecovered);
        // Add flow directly to Recovered compartment.
        dydt[LctState::template get_first_index<InfectionState::Recovered>()] +=
            -(params.template get<TransitionMatrixInfectedSymptomsToRecovered>() *
              Eigen::VectorX<ScalarType>::Ones(dimensionInfectedSymptomsToRecovered))
                 .transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>() +
                          dimensionInfectedSymptomsToInfectedSevere,
                      dimensionInfectedSymptomsToRecovered);

        // --- InfectedSevere. ---
        // Flow from InfectedSymptoms To InfectedSevere.
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedSevere>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedSevere>()) =
            -params.template get<StartingProbabilitiesInfectedSevere>() *
            (params.template get<TransitionMatrixInfectedSymptomsToInfectedSevere>() *
             Eigen::VectorX<ScalarType>::Ones(dimensionInfectedSymptomsToInfectedSevere))
                .transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedSymptoms>(),
                      dimensionInfectedSymptomsToInfectedSevere);
        // Flow from InfectedSevere To InfectedCritical.
        size_t dimensionInfectedSevereToInfectedCritical =
            params.template get<TransitionMatrixInfectedSevereToInfectedCritical>().rows();
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedSevere>(),
                     dimensionInfectedSevereToInfectedCritical) +=
            params.template get<TransitionMatrixInfectedSevereToInfectedCritical>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedSevere>(),
                      dimensionInfectedSevereToInfectedCritical);
        // Flow from InfectedSevere To Recovered.
        size_t dimensionInfectedSevereToRecovered =
            params.template get<TransitionMatrixInfectedSevereToRecovered>().rows();
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedSevere>() +
                         dimensionInfectedSevereToInfectedCritical,
                     dimensionInfectedSevereToRecovered) +=
            params.template get<TransitionMatrixInfectedSevereToRecovered>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedSevere>() +
                          dimensionInfectedSevereToInfectedCritical,
                      dimensionInfectedSevereToRecovered);
        // Add flow directly to Recovered compartment.
        dydt[LctState::template get_first_index<InfectionState::Recovered>()] +=
            -(params.template get<TransitionMatrixInfectedSevereToRecovered>() *
              Eigen::VectorX<ScalarType>::Ones(dimensionInfectedSevereToRecovered))
                 .transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedSevere>() +
                          dimensionInfectedSevereToInfectedCritical,
                      dimensionInfectedSevereToRecovered);

        // --- InfectedCritical. ---
        // Flow from InfectedSevere To InfectedCritical.
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedCritical>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedCritical>()) =
            -params.template get<StartingProbabilitiesInfectedCritical>() *
            (params.template get<TransitionMatrixInfectedSevereToInfectedCritical>() *
             Eigen::VectorX<ScalarType>::Ones(dimensionInfectedSevereToInfectedCritical))
                .transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedSevere>(),
                      dimensionInfectedSevereToInfectedCritical);
        // Flow from InfectedCritical To Dead.
        size_t dimensionInfectedCriticalToDead = params.template get<TransitionMatrixInfectedCriticalToDead>().rows();
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedCritical>(),
                     dimensionInfectedCriticalToDead) +=
            params.template get<TransitionMatrixInfectedCriticalToDead>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedCritical>(),
                      dimensionInfectedCriticalToDead);
        // Flow from InfectedCritical To Recovered.
        size_t dimensionInfectedCriticalToRecovered =
            params.template get<TransitionMatrixInfectedCriticalToRecovered>().rows();
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedCritical>() +
                         dimensionInfectedCriticalToDead,
                     dimensionInfectedCriticalToRecovered) +=
            params.template get<TransitionMatrixInfectedCriticalToRecovered>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedCritical>() +
                          dimensionInfectedCriticalToDead,
                      dimensionInfectedCriticalToRecovered);
        // Add flow directly to Recovered compartment.
        dydt[LctState::template get_first_index<InfectionState::Recovered>()] +=
            -(params.template get<TransitionMatrixInfectedCriticalToRecovered>() *
              Eigen::VectorX<ScalarType>::Ones(dimensionInfectedCriticalToRecovered))
                 .transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedCritical>() +
                          dimensionInfectedCriticalToDead,
                      dimensionInfectedCriticalToRecovered);

        // --- Dead. ---
        dydt[LctState::template get_first_index<InfectionState::Dead>()] =
            -(params.template get<TransitionMatrixInfectedCriticalToDead>() *
              Eigen::VectorX<ScalarType>::Ones(dimensionInfectedCriticalToDead))
                 .transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedCritical>(),
                      dimensionInfectedCriticalToDead);
    }

    /**
     * @brief Cumulates a simulation result with subcompartments to produce a result that divides the population only
     *   into the infection states defined in InfectionState.
     *
     * If the model is used for simulation, we will get a result in form of a TimeSeries with infection states divided 
     * in subcompartments.
     * The function calculates a TimeSeries without subcompartments from another TimeSeries with subcompartments. 
     * This is done by summing up the corresponding subcompartments.
     * @param[in] subcompartments_ts Result of a simulation with the model.
     * @return Result of the simulation divided in infection states without subcompartments. 
     *  Returns TimeSeries with values -1 if calculation is not possible.
     */
    TimeSeries<ScalarType> calculate_compartments(const TimeSeries<ScalarType>& subcompartments_ts) const
    {
        TimeSeries<ScalarType> compartments_ts((Eigen::Index)InfectionState::Count);
        if (!(LctState::Count == subcompartments_ts.get_num_elements())) {
            log_error("Result does not match InfectionStates of the model.");
            // Return a TimeSeries with values -1.
            Eigen::VectorX<ScalarType> error_output =
                Eigen::VectorX<ScalarType>::Constant((Eigen::Index)InfectionState::Count, -1);
            compartments_ts.add_time_point(-1, error_output);
            return compartments_ts;
        }
        Eigen::VectorX<ScalarType> compartments((Eigen::Index)InfectionState::Count);
        for (Eigen::Index i = 0; i < subcompartments_ts.get_num_time_points(); ++i) {
            compartments_ts.add_time_point(subcompartments_ts.get_time(i),
                                           LctState::calculate_compartments(subcompartments_ts[i]));
        }
        return compartments_ts;
    }
};

} // namespace glsecir
} // namespace mio

#endif // MIO_GLCT_SECIR_MODEL_H
