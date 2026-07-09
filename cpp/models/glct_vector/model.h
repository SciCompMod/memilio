/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Annika Jungklaus
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

#ifndef MIO_GLCT_VECTOR_MODEL_H
#define MIO_GLCT_VECTOR_MODEL_H

#include "glct_vector/parameters.h"
#include "glct_vector/infection_state.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/compartments/compartmental_model.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/config.h"
#include "memilio/utils/index.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"

#include <cstdlib>
#include <numbers>

namespace mio
{
namespace glvector
{

/**
 * @brief Class that defines an GLCT-VECTOR model.
 *
 * @tparam NumExposedHuman The number of subcompartments used for the ExposedHuman compartment.
 * @tparam NumInfectedHuman The number of subcompartments used for the InfectedHuman compartment.
 * @tparam NumExposedVector The number of subcompartments used for the ExposedVector compartment.
 */
template <typename FP, size_t NumExposedHuman, size_t NumInfectedHuman, size_t NumExposedVector>
class Model
    : public CompartmentalModel<FP, LctInfectionState<FP, InfectionState, 1, 1, NumExposedHuman, NumInfectedHuman, 1, 1, 1, 1, 1, 1,
                            NumExposedVector, 1, 1>,
                            mio::Populations<FP, LctInfectionState<FP, InfectionState, 1, 1, NumExposedHuman, NumInfectedHuman, 1, 1, 1, 1, 1, 1,
                            NumExposedVector, 1, 1>>,
                            Parameters<FP>>
{

public:
    using LctState = LctInfectionState<FP, InfectionState, 1, 1, NumExposedHuman, NumInfectedHuman, 1, 1, 1, 1, 1, 1,
                            NumExposedVector, 1, 1>; ///< This class specifies the number of subcompartments.
    using Base     = CompartmentalModel<FP, LctState, mio::Populations<FP, LctState>, Parameters<FP>>;
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
        if ((Eigen::Index)LctState::template get_num_subcompartments<InfectionState::ExposedHuman>() !=
            params.template get<StartingProbabilitiesExposedHuman<FP>>().rows()) {
            log_error("Constraint check: Dimension of the parameters does not match the number of subcompartments for "
                      "the Exposed_Human "
                      "compartment.");
            return true;
        }
                if ((Eigen::Index)LctState::template get_num_subcompartments<InfectionState::ExposedVector>() !=
            params.template get<StartingProbabilitiesExposedVector<FP>>().rows()) {
            log_error("Constraint check: Dimension of the parameters does not match the number of subcompartments for "
                      "the Exposed_Vector "
                      "compartment.");
            return true;
        }
        if ((Eigen::Index)LctState::template get_num_subcompartments<InfectionState::InfectedHuman>() !=
            params.template get<StartingProbabilitiesInfectedHuman<FP>>().rows()) {
            log_error(
                "Constraint check: Dimension of the parameters does not match the number of subcompartments for the "
                "Infected_Human compartment.");
            return true;
        }
        return (params.check_constraints() || this->populations.check_constraints());
    }

    /**
     * @brief Evaluates the right-hand-side f of the GLCT dydt = f(y, t).
     *
     * The GLCT-SEIR-VECTOR model is defined through ordinary differential equations of the form dydt = f(y, t).
     * y is a vector containing number of individuals for each (sub-) compartment.
     * This function evaluates the right-hand-side f of the ODE and can be used in an ODE solver.
     *
     * @param[in] pop The current state of the population (human and vector) in the geographic unit we are considering.
     * @param[in] y The current state of the model (or a subpopulation) as a flat array.
     * @param[in] t The current time.
     * @param[out] dydt A reference to the calculated output.
     */
    void get_derivatives(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                         Eigen::Ref<Eigen::VectorX<FP>> dydt) const override
    {
        using std::sin;

        dydt.setZero();

        auto params           = this->parameters;

        // Total current human population not counting the dead.
        FP total_population_human = pop[LctState::template get_first_index<InfectionState::SusceptibleHuman>()] 
                                + pop[LctState::template get_first_index<InfectionState::RecoveredHuman>()] 
                                + pop.segment(LctState::template get_first_index<InfectionState::ExposedHuman>(),
                                  LctState::template get_num_subcompartments<InfectionState::ExposedHuman>()).sum()
                                + pop.segment(LctState::template get_first_index<InfectionState::InfectedHuman>(),
                                  LctState::template get_num_subcompartments<InfectionState::InfectedHuman>()).sum();
        // Initial human population. Here we assume 0 dead and unborn people at the start of the simulation.
        FP initial_population_human = pop[LctState::template get_first_index<InfectionState::UnbornHuman>()]
                                + pop[LctState::template get_first_index<InfectionState::SusceptibleHuman>()] 
                                + pop[LctState::template get_first_index<InfectionState::RecoveredHuman>()] 
                                + pop[LctState::template get_first_index<InfectionState::DeadHuman>()] 
                                + pop[LctState::template get_first_index<InfectionState::NaturalDeathHuman>()] 
                                + pop.segment(LctState::template get_first_index<InfectionState::ExposedHuman>(),
                                  LctState::template get_num_subcompartments<InfectionState::ExposedHuman>()).sum()
                                + pop.segment(LctState::template get_first_index<InfectionState::InfectedHuman>(),
                                  LctState::template get_num_subcompartments<InfectionState::InfectedHuman>()).sum();

        // Calculate sum of all subcompartments for InfectedHumans.
        FP InfectedHuman_sum =
            pop.segment(LctState::template get_first_index<InfectionState::InfectedHuman>(),
                        LctState::template get_num_subcompartments<InfectionState::InfectedHuman>())
                .sum();
        // Calculate sum of all subcompartments for TransmitterVector.
        FP TransmitterVector_sum =
            pop.segment(LctState::template get_first_index<InfectionState::TransmitterVector>(),
                        LctState::template get_num_subcompartments<InfectionState::TransmitterVector>())
                .sum();

        // --- Vector compartments. ---
        // --- AquaticVector. ---
        // We only model females in the vector population (with the parameter SexRatioVector scaling the birth of new individuals)
        FP AdultFemalesVector = pop[LctState::template get_first_index<InfectionState::SusceptibleVector>()] 
                                + pop.segment(LctState::template get_first_index<InfectionState::ExposedVector>(),
                                  LctState::template get_num_subcompartments<InfectionState::ExposedVector>()).sum()
                                + pop.segment(LctState::template get_first_index<InfectionState::TransmitterVector>(),
                                  LctState::template get_num_subcompartments<InfectionState::TransmitterVector>()).sum();
        // Birth, only considering females.
        // If the system is overpopulated, i.e. the vector population exceeds the carrying capacity, the amount born can be negative (death of larvae).
        FP seasonal_variation = initial_population_human * params.template get<RatioVectorToHuman<FP>>() * (1 + 
                                params.template get<SeasonalVariationImpact<FP>>() * sin(2 * std::numbers::pi_v<ScalarType> *
                                    t / params.template get<SeasonalCycleDuration<FP>>()));
        dydt[LctState::template get_first_index<InfectionState::AquaticVector>()] = params.template get<OvipositionRateVector<FP>>()
                                * (1 - AdultFemalesVector / seasonal_variation) * params.template get<SexRatioVector<FP>>()
                                * AdultFemalesVector;
        dydt[LctState::template get_first_index<InfectionState::UnbornVector>()] = - dydt[LctState::template get_first_index<InfectionState::AquaticVector>()];
        // Death.
        dydt[LctState::template get_first_index<InfectionState::AquaticVector>()] -= 
                                params.template get<MortalityRateAquaticVector<FP>>()
                                * y[LctState::template get_first_index<InfectionState::AquaticVector>()];
        dydt[LctState::template get_first_index<InfectionState::DeadVector>()] = - dydt[LctState::template get_first_index<InfectionState::AquaticVector>()]
                                - dydt[LctState::template get_first_index<InfectionState::UnbornVector>()];
        // Aging into adult stage.
        dydt[LctState::template get_first_index<InfectionState::AquaticVector>()] -=
                                params.template get<TransitionRateAquaticToAdultVector<FP>>()
                                * y[LctState::template get_first_index<InfectionState::AquaticVector>()];

        // --- SusceptibleVector. ---
        // Death.
        dydt[LctState::template get_first_index<InfectionState::SusceptibleVector>()] = - params.template get<MortalityRateAdultVector<FP>>()
                                * y[LctState::template get_first_index<InfectionState::SusceptibleVector>()];
        dydt[LctState::template get_first_index<InfectionState::DeadVector>()] -= dydt[LctState::template get_first_index<InfectionState::SusceptibleVector>()];
        // Inflow from the aquatic stage
        dydt[LctState::template get_first_index<InfectionState::SusceptibleVector>()] += params.template get<TransitionRateAquaticToAdultVector<FP>>() 
                                * y[LctState::template get_first_index<InfectionState::AquaticVector>()];
        // Outflow to ExposedVector.
        dydt[LctState::template get_first_index<InfectionState::SusceptibleVector>()] -= params.template get<BitingRateVector<FP>>()  *
                                params.template get<TransmissionProbabilityOnContactHumanToVector<FP>>() * InfectedHuman_sum / total_population_human
                                * y[LctState::template get_first_index<InfectionState::SusceptibleVector>()];

        // --- ExposedVector. ---
        // Death.
        dydt.segment(LctState::template get_first_index<InfectionState::ExposedVector>(),
                     LctState::template get_num_subcompartments<InfectionState::ExposedVector>()) = - params.template get<MortalityRateAdultVector<FP>>()
                            * y.segment(LctState::template get_first_index<InfectionState::ExposedVector>(), 
                            LctState::template get_num_subcompartments<InfectionState::ExposedVector>());
        dydt[LctState::template get_first_index<InfectionState::DeadVector>()] += - dydt.segment(LctState::template get_first_index<InfectionState::ExposedVector>(),
                     LctState::template get_num_subcompartments<InfectionState::ExposedVector>()).sum();
        // Inflow from SusceptibleVector.
        dydt.segment(LctState::template get_first_index<InfectionState::ExposedVector>(),
                     LctState::template get_num_subcompartments<InfectionState::ExposedVector>()) += params.template get<BitingRateVector<FP>>()
                 * params.template get<TransmissionProbabilityOnContactHumanToVector<FP>>() * InfectedHuman_sum / total_population_human
                * y[LctState::template get_first_index<InfectionState::SusceptibleVector>()] * params.template get<StartingProbabilitiesExposedVector<FP>>();
        // Outflow to TransmitterVector.
        dydt.segment(LctState::template get_first_index<InfectionState::ExposedVector>(),
                     LctState::template get_num_subcompartments<InfectionState::ExposedVector>()) +=
            params.template get<TransitionMatrixExposedToTransmitterVector<FP>>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::ExposedVector>(),
                      LctState::template get_num_subcompartments<InfectionState::ExposedVector>());

        // --- TransmitterVector. ---
        // Death.
        dydt[LctState::template get_first_index<InfectionState::TransmitterVector>()] = - params.template get<MortalityRateAdultVector<FP>>()
                            * y[LctState::template get_first_index<InfectionState::TransmitterVector>()];
        dydt[LctState::template get_first_index<InfectionState::DeadVector>()] -= dydt.segment(LctState::template get_first_index<InfectionState::TransmitterVector>(),
                     LctState::template get_num_subcompartments<InfectionState::TransmitterVector>()).sum();
        // Inflow from ExposedVector.
        dydt[LctState::template get_first_index<InfectionState::TransmitterVector>()] -=
            (params.template get<TransitionMatrixExposedToTransmitterVector<FP>>() *
              Eigen::VectorX<FP>::Ones(LctState::template get_num_subcompartments<InfectionState::ExposedVector>()))
                 .transpose() *
            y.segment(LctState::template get_first_index<InfectionState::ExposedVector>(),
                      LctState::template get_num_subcompartments<InfectionState::ExposedVector>());


        // --- Human compartments. ---
        // --- SusceptibleHuman. ---
        // Birth.
        dydt[LctState::template get_first_index<InfectionState::SusceptibleHuman>()] = params.template get<BirthRateHuman<FP>>() 
                            * total_population_human;
        dydt[LctState::template get_first_index<InfectionState::UnbornHuman>()] = - dydt[LctState::template get_first_index<InfectionState::SusceptibleHuman>()];
        // Natural death.
        dydt[LctState::template get_first_index<InfectionState::SusceptibleHuman>()] -= params.template get<MortalityRateHuman<FP>>() 
                            * y[LctState::template get_first_index<InfectionState::SusceptibleHuman>()];
        dydt[LctState::template get_first_index<InfectionState::NaturalDeathHuman>()] = params.template get<MortalityRateHuman<FP>>() 
                            * y[LctState::template get_first_index<InfectionState::SusceptibleHuman>()];       
        // Outflow to ExposedHuman.
        dydt[LctState::template get_first_index<InfectionState::SusceptibleHuman>()] -=
            y[LctState::template get_first_index<InfectionState::SusceptibleHuman>()] / total_population_human * params.template get<BitingRateVector<FP>>() * 
            params.template get<TransmissionProbabilityOnContactVectorToHuman<FP>>() * TransmitterVector_sum;

        // --- ExposedHuman. ---
        // Natural death.
        dydt.segment(LctState::template get_first_index<InfectionState::ExposedHuman>(),
                     LctState::template get_num_subcompartments<InfectionState::ExposedHuman>()) = - params.template get<MortalityRateHuman<FP>>()
                            * y.segment(LctState::template get_first_index<InfectionState::ExposedHuman>(), 
                            LctState::template get_num_subcompartments<InfectionState::ExposedHuman>());
        dydt[LctState::template get_first_index<InfectionState::NaturalDeathHuman>()] -= dydt.segment(LctState::template get_first_index<InfectionState::ExposedHuman>(),
                     LctState::template get_num_subcompartments<InfectionState::ExposedHuman>()).sum();
        // Inflow from SusceptibleHuman.
        dydt.segment(LctState::template get_first_index<InfectionState::ExposedHuman>(),
                     LctState::template get_num_subcompartments<InfectionState::ExposedHuman>()) +=
                y[LctState::template get_first_index<InfectionState::SusceptibleHuman>()] / total_population_human * params.template get<BitingRateVector<FP>>() * 
            params.template get<TransmissionProbabilityOnContactVectorToHuman<FP>>() * TransmitterVector_sum * params.template get<StartingProbabilitiesExposedHuman<FP>>();
        // Outflow to InfectedHuman.
        dydt.segment(LctState::template get_first_index<InfectionState::ExposedHuman>(),
                     LctState::template get_num_subcompartments<InfectionState::ExposedHuman>()) +=
            params.template get<TransitionMatrixExposedToInfectedHuman<FP>>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::ExposedHuman>(),
                      LctState::template get_num_subcompartments<InfectionState::ExposedHuman>());

        // --- InfectedHuman. ---
        // Natural death.
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedHuman>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedHuman>()) = - params.template get<MortalityRateHuman<FP>>()
                            * y.segment(LctState::template get_first_index<InfectionState::InfectedHuman>(), 
                            LctState::template get_num_subcompartments<InfectionState::InfectedHuman>());
        dydt[LctState::template get_first_index<InfectionState::NaturalDeathHuman>()] -= dydt.segment(LctState::template get_first_index<InfectionState::InfectedHuman>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedHuman>()).sum();
        // Inflow from ExpsoedHuman.
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedHuman>(),
                     LctState::template get_num_subcompartments<InfectionState::InfectedHuman>()) -=
            (params.template get<TransitionMatrixExposedToInfectedHuman<FP>>() *
              Eigen::VectorX<FP>::Ones(LctState::template get_num_subcompartments<InfectionState::ExposedHuman>())).transpose()
            * y.segment(LctState::template get_first_index<InfectionState::ExposedHuman>(),
                      LctState::template get_num_subcompartments<InfectionState::ExposedHuman>()) *
            params.template get<StartingProbabilitiesInfectedHuman<FP>>();
        // The outflow has to be split between RecoveredHuman and DeadHuman
        size_t dimensionInfectedToRecovered =
            params.template get<TransitionMatrixInfectedToRecoveredHuman<FP>>().rows();
        size_t dimensionInfectedToDead =
            params.template get<TransitionMatrixInfectedToDeadHuman<FP>>().rows();
        // Outflow to RecoveredHuman.
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedHuman>(),
                     dimensionInfectedToRecovered) +=
            params.template get<TransitionMatrixInfectedToRecoveredHuman<FP>>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedHuman>(),
                      dimensionInfectedToRecovered);
        // Outflow to DeadHuman.
        dydt.segment(LctState::template get_first_index<InfectionState::InfectedHuman>() + dimensionInfectedToRecovered,
                     dimensionInfectedToDead) +=
            params.template get<TransitionMatrixInfectedToDeadHuman<FP>>().transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedHuman>() + dimensionInfectedToRecovered,
                      dimensionInfectedToDead);

        // --- RecoveredHuman. ---
        // Natural death.
        dydt[LctState::template get_first_index<InfectionState::RecoveredHuman>()] = - params.template get<MortalityRateHuman<FP>>() 
                            * y[LctState::template get_first_index<InfectionState::RecoveredHuman>()];
        dydt[LctState::template get_first_index<InfectionState::NaturalDeathHuman>()] -= dydt[LctState::template get_first_index<InfectionState::RecoveredHuman>()];
        // Inflow from InfectedHuman.
        dydt[LctState::template get_first_index<InfectionState::RecoveredHuman>()] -=
            (params.template get<TransitionMatrixInfectedToRecoveredHuman<FP>>() *
              Eigen::VectorX<FP>::Ones(dimensionInfectedToRecovered)).transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedHuman>(),
                      dimensionInfectedToRecovered);

        // --- DeadHuman. ---
        // Inflow from InfectedHuman.
        dydt[LctState::template get_first_index<InfectionState::DeadHuman>()] -=
            (params.template get<TransitionMatrixInfectedToDeadHuman<FP>>() *
              Eigen::VectorX<FP>::Ones(dimensionInfectedToDead)).transpose() *
            y.segment(LctState::template get_first_index<InfectionState::InfectedHuman>() +
                          dimensionInfectedToRecovered,
                      dimensionInfectedToDead);
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
    TimeSeries<FP> calculate_compartments(const TimeSeries<FP>& subcompartments_ts) const
    {
        TimeSeries<FP> compartments_ts((Eigen::Index)InfectionState::Count);
        if (!(LctState::Count == subcompartments_ts.get_num_elements())) {
            log_error("Result does not match InfectionStates of the model.");
            // Return a TimeSeries with values -1.
            Eigen::VectorX<FP> error_output = Eigen::VectorX<FP>::Constant((Eigen::Index)InfectionState::Count, -1);
            compartments_ts.add_time_point(-1, error_output);
            return compartments_ts;
        }
        Eigen::VectorX<FP> compartments((Eigen::Index)InfectionState::Count);
        for (Eigen::Index i = 0; i < subcompartments_ts.get_num_time_points(); ++i) {
            compartments_ts.add_time_point(subcompartments_ts.get_time(i),
                                           LctState::calculate_compartments(subcompartments_ts[i]));
        }
        return compartments_ts;
    }
};

} // namespace glvector
} // namespace mio

#endif // MIO_GLCT_VECTOR_MODEL_H
