/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Wadim Koslow, Daniel Abele, Martin J. Kühn
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
#ifndef MIO_ODE_SECIRVVS_MODEL_H
#define MIO_ODE_SECIRVVS_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/epidemiology/populations.h"
#include "ode_secirvvs/analyze_result.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/parameters.h"
#include "memilio/math/smoother.h"
#include "memilio/math/eigen_util.h"
#include "memilio/math/interpolation.h"

#include <numbers>
#include <complex>

namespace mio
{
namespace osecirvvs
{
// clang-format off
using Flows = TypeList<
    //naive
    Flow<InfectionState::SusceptibleNaive,                            InfectionState::ExposedNaive>,
    Flow<InfectionState::ExposedNaive,                                InfectionState::InfectedNoSymptomsNaive>,
    Flow<InfectionState::InfectedNoSymptomsNaive,                     InfectionState::InfectedSymptomsNaive>,
    Flow<InfectionState::InfectedNoSymptomsNaive,                     InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedNoSymptomsNaiveConfirmed,            InfectionState::InfectedSymptomsNaiveConfirmed>,
    Flow<InfectionState::InfectedNoSymptomsNaiveConfirmed,            InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsNaive,                       InfectionState::InfectedSevereNaive>,
    Flow<InfectionState::InfectedSymptomsNaive,                       InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsNaiveConfirmed,              InfectionState::InfectedSevereNaive>,
    Flow<InfectionState::InfectedSymptomsNaiveConfirmed,              InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedSevereNaive,                         InfectionState::InfectedCriticalNaive>,
    Flow<InfectionState::InfectedSevereNaive,                         InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedSevereNaive,                         InfectionState::DeadNaive>,
    Flow<InfectionState::InfectedCriticalNaive,                       InfectionState::DeadNaive>,
    Flow<InfectionState::InfectedCriticalNaive,                       InfectionState::SusceptibleImprovedImmunity>,
    //partial immunity
    Flow<InfectionState::SusceptiblePartialImmunity,                  InfectionState::ExposedPartialImmunity>,
    Flow<InfectionState::ExposedPartialImmunity,                      InfectionState::InfectedNoSymptomsPartialImmunity>,
    Flow<InfectionState::InfectedNoSymptomsPartialImmunity,           InfectionState::InfectedSymptomsPartialImmunity>,
    Flow<InfectionState::InfectedNoSymptomsPartialImmunity,           InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,  InfectionState::InfectedSymptomsPartialImmunityConfirmed>,
    Flow<InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,  InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsPartialImmunity,             InfectionState::InfectedSeverePartialImmunity>,
    Flow<InfectionState::InfectedSymptomsPartialImmunity,             InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsPartialImmunityConfirmed,    InfectionState::InfectedSeverePartialImmunity>,
    Flow<InfectionState::InfectedSymptomsPartialImmunityConfirmed,    InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedSeverePartialImmunity,               InfectionState::InfectedCriticalPartialImmunity>,
    Flow<InfectionState::InfectedSeverePartialImmunity,               InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedSeverePartialImmunity,               InfectionState::DeadPartialImmunity>,
    Flow<InfectionState::InfectedCriticalPartialImmunity,             InfectionState::DeadPartialImmunity>,
    Flow<InfectionState::InfectedCriticalPartialImmunity,             InfectionState::SusceptibleImprovedImmunity>,
    //improved immunity
    Flow<InfectionState::SusceptibleImprovedImmunity,                 InfectionState::ExposedImprovedImmunity>,
    Flow<InfectionState::ExposedImprovedImmunity,                     InfectionState::InfectedNoSymptomsImprovedImmunity>,
    Flow<InfectionState::InfectedNoSymptomsImprovedImmunity,          InfectionState::InfectedSymptomsImprovedImmunity>,
    Flow<InfectionState::InfectedNoSymptomsImprovedImmunity,          InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed, InfectionState::InfectedSymptomsImprovedImmunityConfirmed>,
    Flow<InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed, InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsImprovedImmunity,            InfectionState::InfectedSevereImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsImprovedImmunity,            InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsImprovedImmunityConfirmed,   InfectionState::InfectedSevereImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsImprovedImmunityConfirmed,   InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedSevereImprovedImmunity,              InfectionState::InfectedCriticalImprovedImmunity>,
    Flow<InfectionState::InfectedSevereImprovedImmunity,              InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::InfectedSevereImprovedImmunity,              InfectionState::DeadImprovedImmunity>,
    Flow<InfectionState::InfectedCriticalImprovedImmunity,            InfectionState::DeadImprovedImmunity>,
    Flow<InfectionState::InfectedCriticalImprovedImmunity,            InfectionState::SusceptibleImprovedImmunity>>;
// clang-format on

template <typename FP>
class Model
    : public FlowModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>
{
    using Base = FlowModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>;

public:
    using typename Base::ParameterSet;
    using typename Base::Populations;

    Model(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }

    Model(int num_agegroups)
        : Model(Populations({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        auto const& params   = this->parameters;
        AgeGroup n_agegroups = params.get_num_groups();

        ContactMatrixGroup<FP> const& contact_matrix = params.template get<ContactPatterns<FP>>();

        FP icu_occupancy           = 0.0;
        FP test_and_trace_required = 0.0;
        for (auto i = AgeGroup(0); i < n_agegroups; ++i) {
            test_and_trace_required +=
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                params.template get<TimeInfectedNoSymptoms<FP>>()[i] *
                (this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptomsNaive}) +
                 this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptomsPartialImmunity}) +
                 this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptomsImprovedImmunity}) +
                 this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptomsNaiveConfirmed}) +
                 this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}) +
                 this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}));
            icu_occupancy += this->populations.get_from(pop, {i, InfectionState::InfectedCriticalNaive}) +
                             this->populations.get_from(pop, {i, InfectionState::InfectedCriticalPartialImmunity}) +
                             this->populations.get_from(pop, {i, InfectionState::InfectedCriticalImprovedImmunity});
        }

        for (auto i = AgeGroup(0); i < n_agegroups; i++) {

            size_t SNi    = this->populations.get_flat_index({i, InfectionState::SusceptibleNaive});
            size_t ENi    = this->populations.get_flat_index({i, InfectionState::ExposedNaive});
            size_t INSNi  = this->populations.get_flat_index({i, InfectionState::InfectedNoSymptomsNaive});
            size_t ISyNi  = this->populations.get_flat_index({i, InfectionState::InfectedSymptomsNaive});
            size_t ISevNi = this->populations.get_flat_index({i, InfectionState::InfectedSevereNaive});
            size_t ICrNi  = this->populations.get_flat_index({i, InfectionState::InfectedCriticalNaive});

            size_t INSNCi = this->populations.get_flat_index({i, InfectionState::InfectedNoSymptomsNaiveConfirmed});
            size_t ISyNCi = this->populations.get_flat_index({i, InfectionState::InfectedSymptomsNaiveConfirmed});

            size_t SPIi    = this->populations.get_flat_index({i, InfectionState::SusceptiblePartialImmunity});
            size_t EPIi    = this->populations.get_flat_index({i, InfectionState::ExposedPartialImmunity});
            size_t INSPIi  = this->populations.get_flat_index({i, InfectionState::InfectedNoSymptomsPartialImmunity});
            size_t ISyPIi  = this->populations.get_flat_index({i, InfectionState::InfectedSymptomsPartialImmunity});
            size_t ISevPIi = this->populations.get_flat_index({i, InfectionState::InfectedSeverePartialImmunity});
            size_t ICrPIi  = this->populations.get_flat_index({i, InfectionState::InfectedCriticalPartialImmunity});

            size_t INSPICi =
                this->populations.get_flat_index({i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed});
            size_t ISyPICi =
                this->populations.get_flat_index({i, InfectionState::InfectedSymptomsPartialImmunityConfirmed});

            size_t EIIi    = this->populations.get_flat_index({i, InfectionState::ExposedImprovedImmunity});
            size_t INSIIi  = this->populations.get_flat_index({i, InfectionState::InfectedNoSymptomsImprovedImmunity});
            size_t ISyIIi  = this->populations.get_flat_index({i, InfectionState::InfectedSymptomsImprovedImmunity});
            size_t ISevIIi = this->populations.get_flat_index({i, InfectionState::InfectedSevereImprovedImmunity});
            size_t ICrIIi  = this->populations.get_flat_index({i, InfectionState::InfectedCriticalImprovedImmunity});

            size_t INSIICi =
                this->populations.get_flat_index({i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed});
            size_t ISyIICi =
                this->populations.get_flat_index({i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed});

            size_t SIIi = this->populations.get_flat_index({i, InfectionState::SusceptibleImprovedImmunity});

            FP reducExposedPartialImmunity  = params.template get<ReducExposedPartialImmunity<FP>>()[i];
            FP reducExposedImprovedImmunity = params.template get<ReducExposedImprovedImmunity<FP>>()[i];
            FP reducInfectedSymptomsPartialImmunity =
                params.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[i];
            FP reducInfectedSymptomsImprovedImmunity =
                params.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[i];
            FP reducInfectedSevereCriticalDeadPartialImmunity =
                params.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[i];
            FP reducInfectedSevereCriticalDeadImprovedImmunity =
                params.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[i];
            FP reducTimeInfectedMild = params.template get<ReducTimeInfectedMild<FP>>()[i];

            //symptomatic are less well quarantined when testing and tracing is overwhelmed so they infect more people
            auto riskFromInfectedSymptomatic =
                smoother_cosine<FP>(test_and_trace_required, params.template get<TestAndTraceCapacity<FP>>(),
                                    params.template get<TestAndTraceCapacity<FP>>() *
                                        params.template get<TestAndTraceCapacityMaxRiskSymptoms<FP>>(),
                                    params.template get<RiskOfInfectionFromSymptomatic<FP>>()[i],
                                    params.template get<MaxRiskOfInfectionFromSymptomatic<FP>>()[i]);

            auto riskFromInfectedNoSymptoms =
                smoother_cosine<FP>(test_and_trace_required, params.template get<TestAndTraceCapacity<FP>>(),
                                    params.template get<TestAndTraceCapacity<FP>>() *
                                        params.template get<TestAndTraceCapacityMaxRiskNoSymptoms<FP>>(),
                                    params.template get<RelativeTransmissionNoSymptoms<FP>>()[i], 1.0);

            for (auto j = AgeGroup(0); j < n_agegroups; j++) {
                size_t SNj    = this->populations.get_flat_index({j, InfectionState::SusceptibleNaive});
                size_t ENj    = this->populations.get_flat_index({j, InfectionState::ExposedNaive});
                size_t INSNj  = this->populations.get_flat_index({j, InfectionState::InfectedNoSymptomsNaive});
                size_t ISyNj  = this->populations.get_flat_index({j, InfectionState::InfectedSymptomsNaive});
                size_t ISevNj = this->populations.get_flat_index({j, InfectionState::InfectedSevereNaive});
                size_t ICrNj  = this->populations.get_flat_index({j, InfectionState::InfectedCriticalNaive});
                size_t SIIj   = this->populations.get_flat_index({j, InfectionState::SusceptibleImprovedImmunity});

                size_t INSNCj = this->populations.get_flat_index({j, InfectionState::InfectedNoSymptomsNaiveConfirmed});
                size_t ISyNCj = this->populations.get_flat_index({j, InfectionState::InfectedSymptomsNaiveConfirmed});

                size_t SPIj = this->populations.get_flat_index({j, InfectionState::SusceptiblePartialImmunity});
                size_t EPIj = this->populations.get_flat_index({j, InfectionState::ExposedPartialImmunity});
                size_t INSPIj =
                    this->populations.get_flat_index({j, InfectionState::InfectedNoSymptomsPartialImmunity});
                size_t ISyPIj  = this->populations.get_flat_index({j, InfectionState::InfectedSymptomsPartialImmunity});
                size_t ISevPIj = this->populations.get_flat_index({j, InfectionState::InfectedSeverePartialImmunity});
                size_t ICrPIj  = this->populations.get_flat_index({j, InfectionState::InfectedCriticalPartialImmunity});

                size_t INSPICj =
                    this->populations.get_flat_index({j, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed});
                size_t ISyPICj =
                    this->populations.get_flat_index({j, InfectionState::InfectedSymptomsPartialImmunityConfirmed});

                size_t EIIj = this->populations.get_flat_index({j, InfectionState::ExposedImprovedImmunity});
                size_t INSIIj =
                    this->populations.get_flat_index({j, InfectionState::InfectedNoSymptomsImprovedImmunity});
                size_t ISyIIj = this->populations.get_flat_index({j, InfectionState::InfectedSymptomsImprovedImmunity});
                size_t ISevIIj = this->populations.get_flat_index({j, InfectionState::InfectedSevereImprovedImmunity});
                size_t ICrIIj = this->populations.get_flat_index({j, InfectionState::InfectedCriticalImprovedImmunity});

                size_t INSIICj =
                    this->populations.get_flat_index({j, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed});
                size_t ISyIICj =
                    this->populations.get_flat_index({j, InfectionState::InfectedSymptomsImprovedImmunityConfirmed});

                // effective contact rate by contact rate between groups i and j and damping j
                FP season_val = (1 + params.template get<Seasonality<FP>>() *
                                         sin(std::numbers::pi_v<ScalarType> *
                                             ((params.template get<StartDay<FP>>() + t) / 182.5 + 0.5)));
                FP cont_freq_eff =
                    season_val * contact_matrix.get_matrix_at(SimulationTime<FP>(t))(
                                     static_cast<Eigen::Index>((size_t)i), static_cast<Eigen::Index>((size_t)j));
                // without died people
                FP Nj = pop[SNj] + pop[ENj] + pop[INSNj] + pop[ISyNj] + pop[ISevNj] + pop[ICrNj] + pop[INSNCj] +
                        pop[ISyNCj] + pop[SPIj] + pop[EPIj] + pop[INSPIj] + pop[ISyPIj] + pop[ISevPIj] + pop[ICrPIj] +
                        pop[INSPICj] + pop[ISyPICj] + pop[SIIj] + pop[EIIj] + pop[INSIIj] + pop[ISyIIj] + pop[ISevIIj] +
                        pop[ICrIIj] + pop[INSIICj] + pop[ISyIICj];
                const FP divNj = (Nj < Limits<FP>::zero_tolerance()) ? FP(0.0) : FP(1.0 / Nj);

                FP ext_inf_force_dummy = cont_freq_eff * divNj *
                                         params.template get<TransmissionProbabilityOnContact<FP>>()[(AgeGroup)i] *
                                         (riskFromInfectedNoSymptoms * (pop[INSNj] + pop[INSPIj] + pop[INSIIj]) +
                                          riskFromInfectedSymptomatic * (pop[ISyNj] + pop[ISyPIj] + pop[ISyIIj]));

                FP dummy_SN = y[SNi] * ext_inf_force_dummy;

                FP dummy_SPI = y[SPIi] * reducExposedPartialImmunity * ext_inf_force_dummy;

                FP dummy_SII = y[SIIi] * reducExposedImprovedImmunity * ext_inf_force_dummy;

                flows[this->template get_flat_flow_index<InfectionState::SusceptibleNaive,
                                                         InfectionState::ExposedNaive>({i})] += dummy_SN;
                flows[this->template get_flat_flow_index<InfectionState::SusceptiblePartialImmunity,
                                                         InfectionState::ExposedPartialImmunity>({i})] += dummy_SPI;
                flows[this->template get_flat_flow_index<InfectionState::SusceptibleImprovedImmunity,
                                                         InfectionState::ExposedImprovedImmunity>({i})] += dummy_SII;
            }

            // ICU capacity shortage is close
            // TODO: if this is used with vaccination model, it has to be adapted if CriticalPerSevere
            // is different for different vaccination status. This is not the case here and in addition, ICUCapacity
            // is set to infinity and this functionality is deactivated, so this is OK for the moment.
            FP criticalPerSevereAdjusted = smoother_cosine<FP>(
                icu_occupancy, 0.90 * params.template get<ICUCapacity<FP>>(), params.template get<ICUCapacity<FP>>(),
                params.template get<CriticalPerSevere<FP>>()[i], 0);

            FP deathsPerSevereAdjusted = params.template get<CriticalPerSevere<FP>>()[i] - criticalPerSevereAdjusted;

            /**** path of immune-naive ***/
            // Exposed
            flows[this->template get_flat_flow_index<InfectionState::ExposedNaive,
                                                     InfectionState::InfectedNoSymptomsNaive>({i})] +=
                y[ENi] / params.template get<TimeExposed<FP>>()[i];

            // InfectedNoSymptoms
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsNaive,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i] /
                params.template get<TimeInfectedNoSymptoms<FP>>()[i] * y[INSNi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsNaive,
                                                     InfectionState::InfectedSymptomsNaive>({i})] =
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                params.template get<TimeInfectedNoSymptoms<FP>>()[i] * y[INSNi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsNaiveConfirmed,
                                                     InfectionState::InfectedSymptomsNaiveConfirmed>({i})] =
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                params.template get<TimeInfectedNoSymptoms<FP>>()[i] * y[INSNCi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsNaiveConfirmed,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i] /
                params.template get<TimeInfectedNoSymptoms<FP>>()[i] * y[INSNCi];

            // InfectedSymptoms
            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsNaive,
                                                     InfectionState::InfectedSevereNaive>({i})] =
                params.template get<SeverePerInfectedSymptoms<FP>>()[i] /
                params.template get<TimeInfectedSymptoms<FP>>()[i] * y[ISyNi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsNaive,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                params.template get<TimeInfectedSymptoms<FP>>()[i] * y[ISyNi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsNaiveConfirmed,
                                                     InfectionState::InfectedSevereNaive>({i})] =
                params.template get<SeverePerInfectedSymptoms<FP>>()[i] /
                params.template get<TimeInfectedSymptoms<FP>>()[i] * y[ISyNCi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsNaiveConfirmed,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                params.template get<TimeInfectedSymptoms<FP>>()[i] * y[ISyNCi];

            // InfectedSevere
            flows[this->template get_flat_flow_index<InfectionState::InfectedSevereNaive,
                                                     InfectionState::InfectedCriticalNaive>({i})] =
                criticalPerSevereAdjusted / params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevNi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSevereNaive,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - params.template get<CriticalPerSevere<FP>>()[i]) /
                params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevNi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSevereNaive, InfectionState::DeadNaive>(
                {i})] = deathsPerSevereAdjusted / params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevNi];

            // InfectedCritical
            flows[this->template get_flat_flow_index<InfectionState::InfectedCriticalNaive, InfectionState::DeadNaive>(
                {i})] = params.template get<DeathsPerCritical<FP>>()[i] /
                        params.template get<TimeInfectedCritical<FP>>()[i] * y[ICrNi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedCriticalNaive,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - params.template get<DeathsPerCritical<FP>>()[i]) /
                params.template get<TimeInfectedCritical<FP>>()[i] * y[ICrNi];

            // /**** path of partially immune (e.g., one dose of vaccination) ***/
            // Exposed
            flows[this->template get_flat_flow_index<InfectionState::ExposedPartialImmunity,
                                                     InfectionState::InfectedNoSymptomsPartialImmunity>({i})] +=
                y[EPIi] / params.template get<TimeExposed<FP>>()[i];

            // InfectedNoSymptoms
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsPartialImmunity,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                         (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i])) /
                (params.template get<TimeInfectedNoSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[INSPIi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsPartialImmunity,
                                                     InfectionState::InfectedSymptomsPartialImmunity>({i})] =
                (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedNoSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[INSPIi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,
                                                     InfectionState::InfectedSymptomsPartialImmunityConfirmed>({i})] =
                (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedNoSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[INSPICi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                         (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i])) /
                (params.template get<TimeInfectedNoSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[INSPICi];

            // InfectedSymptoms
            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsPartialImmunity,
                                                     InfectionState::InfectedSeverePartialImmunity>({i})] =
                reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity *
                params.template get<SeverePerInfectedSymptoms<FP>>()[i] /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyPIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsPartialImmunity,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity) *
                         params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyPIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsPartialImmunityConfirmed,
                                                     InfectionState::InfectedSeverePartialImmunity>({i})] =
                reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity *
                params.template get<SeverePerInfectedSymptoms<FP>>()[i] /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyPICi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsPartialImmunityConfirmed,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity) *
                         params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyPICi];

            // InfectedSevere
            flows[this->template get_flat_flow_index<InfectionState::InfectedSeverePartialImmunity,
                                                     InfectionState::InfectedCriticalPartialImmunity>({i})] =
                reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity *
                criticalPerSevereAdjusted / params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevPIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSeverePartialImmunity,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                         params.template get<CriticalPerSevere<FP>>()[i]) /
                params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevPIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSeverePartialImmunity,
                                                     InfectionState::DeadPartialImmunity>({i})] =
                (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                deathsPerSevereAdjusted / params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevPIi];

            // InfectedCritical
            flows[this->template get_flat_flow_index<InfectionState::InfectedCriticalPartialImmunity,
                                                     InfectionState::DeadPartialImmunity>({i})] =
                (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                params.template get<DeathsPerCritical<FP>>()[i] / params.template get<TimeInfectedCritical<FP>>()[i] *
                y[ICrPIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedCriticalPartialImmunity,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                         params.template get<DeathsPerCritical<FP>>()[i]) /
                params.template get<TimeInfectedCritical<FP>>()[i] * y[ICrPIi];

            // /**** path of twice vaccinated, here called immune although reinfection is possible now ***/
            // Exposed
            flows[this->template get_flat_flow_index<InfectionState::ExposedImprovedImmunity,
                                                     InfectionState::InfectedNoSymptomsImprovedImmunity>({i})] +=
                y[EIIi] / params.template get<TimeExposed<FP>>()[i];

            // InfectedNoSymptoms
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsImprovedImmunity,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                         (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i])) /
                (params.template get<TimeInfectedNoSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[INSIIi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsImprovedImmunity,
                                                     InfectionState::InfectedSymptomsImprovedImmunity>({i})] =
                (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedNoSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[INSIIi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed,
                                                     InfectionState::InfectedSymptomsImprovedImmunityConfirmed>({i})] =
                (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedNoSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[INSIICi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                         (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i])) /
                (params.template get<TimeInfectedNoSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[INSIICi];

            // InfectedSymptoms
            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsImprovedImmunity,
                                                     InfectionState::InfectedSevereImprovedImmunity>({i})] =
                reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity *
                params.template get<SeverePerInfectedSymptoms<FP>>()[i] /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyIIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsImprovedImmunity,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity) *
                         params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyIIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsImprovedImmunityConfirmed,
                                                     InfectionState::InfectedSevereImprovedImmunity>({i})] =
                reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity *
                params.template get<SeverePerInfectedSymptoms<FP>>()[i] /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyIICi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsImprovedImmunityConfirmed,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 - (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity) *
                         params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyIICi];

            // InfectedSevere
            flows[this->template get_flat_flow_index<InfectionState::InfectedSevereImprovedImmunity,
                                                     InfectionState::InfectedCriticalImprovedImmunity>({i})] =
                reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity *
                criticalPerSevereAdjusted / params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevIIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSevereImprovedImmunity,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 -
                 (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                     params.template get<CriticalPerSevere<FP>>()[i]) /
                params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevIIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSevereImprovedImmunity,
                                                     InfectionState::DeadImprovedImmunity>({i})] =
                (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                deathsPerSevereAdjusted / params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevIIi];

            // InfectedCritical
            flows[this->template get_flat_flow_index<InfectionState::InfectedCriticalImprovedImmunity,
                                                     InfectionState::DeadImprovedImmunity>({i})] =
                (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                params.template get<DeathsPerCritical<FP>>()[i] / params.template get<TimeInfectedCritical<FP>>()[i] *
                y[ICrIIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedCriticalImprovedImmunity,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                (1 -
                 (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                     params.template get<DeathsPerCritical<FP>>()[i]) /
                params.template get<TimeInfectedCritical<FP>>()[i] * y[ICrIIi];
        }
    }

    /**
    * serialize this.
    * @see mio::serialize
    */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("Model");
        obj.add_element("Parameters", this->parameters);
        obj.add_element("Populations", this->populations);
    }

    /**
    * deserialize an object of this class.
    * @see mio::deserialize
    */
    template <class IOContext>
    static IOResult<Model> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Model");
        auto par = obj.expect_element("Parameters", Tag<ParameterSet>{});
        auto pop = obj.expect_element("Populations", Tag<Populations>{});
        return apply(
            io,
            [](auto&& par_, auto&& pop_) {
                return Model{pop_, par_};
            },
            par, pop);
    }
};

//forward declaration, see below.
template <typename FP, class Base = mio::Simulation<FP, Model<FP>>>
class Simulation;

/**
* get percentage of infections per total population.
* @param model the compartment model with initial values.
* @param t current simulation time.
* @param y current value of compartments.
* @tparam Base simulation type that uses a secir compartment model. see Simulation.
*/
template <typename FP, class Base = mio::Simulation<FP, Model<FP>>>
FP get_infections_relative(const Simulation<FP, Base>& model, FP t, const Eigen::Ref<const Eigen::VectorX<FP>>& y);

/**
 * specialization of compartment model simulation for the SECIRVVS model.
 * @tparam FP floating point type, e.g., double.
 * @tparam BaseT simulation type, default mio::Simulation. For testing purposes only!
 */
template <typename FP, class BaseT>
class Simulation : public BaseT
{
public:
    /**
     * construct a simulation.
     * @param model the model to simulate.
     * @param t0 start time
     * @param dt time steps
     */
    Simulation(mio::osecirvvs::Model<FP> const& model, FP t0 = 0., FP dt = 0.1)
        : BaseT(model, t0, dt)
        , m_t_last_npi_check(t0)
    {
    }

    /**
    * @brief Applies the effect of a new variant of a disease to the transmission probability of the model.
    *
    * This function adjusts the transmission probability of the disease for each age group based on the share of the new variant.
    * The share of the new variant is calculated based on the time `t` and the start day of the new variant.
    * The transmission probability is then updated for each age group in the model.
    *
    * Based on Equation (35) and (36) in doi.org/10.1371/journal.pcbi.1010054
    *
    * @param [in] t The current time.
    * @param [in] base_infectiousness The base infectiousness of the old variant for each age group.
    */
    void apply_variant(const FP t, const CustomIndexArray<UncertainValue<FP>, AgeGroup> base_infectiousness)
    {
        using std::min;
        using std::pow;

        auto start_day             = this->get_model().parameters.template get<StartDay<FP>>();
        auto start_day_new_variant = this->get_model().parameters.template get<StartDayNewVariant<FP>>();

        if (start_day + t >= start_day_new_variant - 1e-10) {
            const FP days_variant      = t - (start_day_new_variant - start_day);
            const FP share_new_variant = min<FP>(1.0, 0.01 * pow(2, (1. / 7) * days_variant));
            const auto num_groups      = this->get_model().parameters.get_num_groups();
            for (auto i = AgeGroup(0); i < num_groups; ++i) {
                FP new_transmission = (1 - share_new_variant) * base_infectiousness[i] +
                                      share_new_variant * base_infectiousness[i] *
                                          this->get_model().parameters.template get<InfectiousnessNewVariant<FP>>()[i];
                this->get_model().parameters.template get<TransmissionProbabilityOnContact<FP>>()[i] = new_transmission;
            }
        }
    }

    void apply_vaccination(FP t)
    {
        using std::floor;
        auto t_idx        = SimulationDay(size_t(floor(t)));
        auto& params      = this->get_model().parameters;
        size_t num_groups = (size_t)params.get_num_groups();
        auto last_value   = this->get_result().get_last_value();

        auto count = (size_t)InfectionState::Count;
        auto S     = (size_t)InfectionState::SusceptibleNaive;
        auto SV    = (size_t)InfectionState::SusceptiblePartialImmunity;
        auto R     = (size_t)InfectionState::SusceptibleImprovedImmunity;

        for (size_t i = 0; i < num_groups; ++i) {

            FP first_vacc;
            FP full_vacc;
            if (t_idx == SimulationDay(0)) {
                first_vacc = params.template get<DailyPartialVaccinations<FP>>()[{(AgeGroup)i, t_idx}];
                full_vacc  = params.template get<DailyFullVaccinations<FP>>()[{(AgeGroup)i, t_idx}];
            }
            else {
                first_vacc =
                    params.template get<DailyPartialVaccinations<FP>>()[{(AgeGroup)i, t_idx}] -
                    params.template get<DailyPartialVaccinations<FP>>()[{(AgeGroup)i, t_idx - SimulationDay(1)}];
                full_vacc = params.template get<DailyFullVaccinations<FP>>()[{(AgeGroup)i, t_idx}] -
                            params.template get<DailyFullVaccinations<FP>>()[{(AgeGroup)i, t_idx - SimulationDay(1)}];
            }

            if (last_value(count * i + S) - first_vacc < 0) {
                FP corrected = 0.99 * last_value(count * i + S);
                log_warning("too many first vaccinated at time {}: setting first_vacc from {} to {}", t, first_vacc,
                            corrected);
                first_vacc = corrected;
            }

            last_value(count * i + S) -= first_vacc;
            last_value(count * i + SV) += first_vacc;

            if (last_value(count * i + SV) - full_vacc < 0) {
                FP corrected = 0.99 * last_value(count * i + SV);
                log_warning("too many fully vaccinated at time {}: setting full_vacc from {} to {}", t, full_vacc,
                            corrected);
                full_vacc = corrected;
            }

            last_value(count * i + SV) -= full_vacc;
            last_value(count * i + R) += full_vacc;
        }
    }

    /**
     * @brief advance simulation to tmax.
     * Overwrites Simulation::advance and includes a check for dynamic NPIs in regular intervals.
     * @see Simulation::advance
     * @param tmax next stopping point of simulation
     * @return value at tmax
     */
    Eigen::Ref<Eigen::VectorX<FP>> advance(FP tmax)
    {
        using std::floor;
        using std::min;

        auto& t_end_dyn_npis   = this->get_model().parameters.get_end_dynamic_npis();
        auto& dyn_npis         = this->get_model().parameters.template get<DynamicNPIsInfectedSymptoms<FP>>();
        auto& contact_patterns = this->get_model().parameters.template get<ContactPatterns<FP>>();
        // const size_t num_groups = (size_t)this->get_model().parameters.get_num_groups();

        // in the apply_variant function, we adjust the TransmissionProbabilityOnContact parameter. We need to store
        // the base value to use it in the apply_variant function and also to reset the parameter after the simulation.
        auto base_infectiousness = this->get_model().parameters.template get<TransmissionProbabilityOnContact<FP>>();

        FP delay_npi_implementation;
        FP t        = BaseT::get_result().get_last_time();
        const FP dt = dyn_npis.get_thresholds().size() > 0 ? dyn_npis.get_interval().get() : tmax;
        while (t < tmax) {

            FP dt_eff = min<FP>(dt, tmax - t);
            dt_eff    = min<FP>(dt_eff, m_t_last_npi_check + dt - t);

            if (dt_eff >= 1.0) {
                dt_eff = 1.0;
            }

            if (t == 0) {
                //this->apply_vaccination(t); // done in init now?
                this->apply_variant(t, base_infectiousness);
            }
            BaseT::advance(t + dt_eff);
            if (t + 0.5 + dt_eff - floor(t + 0.5) >= 1) {
                this->apply_vaccination(t + 0.5 + dt_eff);
                this->apply_variant(t, base_infectiousness);
            }

            if (t > 0) {
                delay_npi_implementation =
                    this->get_model().parameters.template get<DynamicNPIsImplementationDelay<FP>>();
            }
            else { // DynamicNPIs for t=0 are 'misused' to be from-start NPIs. I.e., do not enforce delay.
                delay_npi_implementation = 0;
            }
            t = t + dt_eff;

            if (dyn_npis.get_thresholds().size() > 0) {
                if (floating_point_greater_equal<FP>(t, m_t_last_npi_check + dt)) {
                    if (t < t_end_dyn_npis) {
                        auto inf_rel = get_infections_relative<FP>(*this, t, this->get_result().get_last_value()) *
                                       dyn_npis.get_base_value();
                        auto exceeded_threshold = dyn_npis.get_max_exceeded_threshold(inf_rel);
                        if (exceeded_threshold != dyn_npis.get_thresholds().end() &&
                            (exceeded_threshold->first > m_dynamic_npi.first ||
                             t > FP(m_dynamic_npi.second))) { //old npi was weaker or is expired

                            auto t_start = SimulationTime<FP>(t + delay_npi_implementation);
                            auto t_end   = t_start + SimulationTime<FP>(dyn_npis.get_duration());
                            this->get_model().parameters.get_start_commuter_detection() = t_start.get();
                            this->get_model().parameters.get_end_commuter_detection()   = t_end.get();
                            m_dynamic_npi = std::make_pair(exceeded_threshold->first, t_end);
                            implement_dynamic_npis(contact_patterns.get_cont_freq_mat(), exceeded_threshold->second,
                                                   t_start, t_end, [](auto& g) {
                                                       return make_contact_damping_matrix(g);
                                                   });
                        }
                    }
                    m_t_last_npi_check = t;
                }
            }
            else {
                m_t_last_npi_check = t;
            }
        }
        // reset TransmissionProbabilityOnContact. This is important for the graph simulation where the advance
        // function is called multiple times for the same model.
        this->get_model().parameters.template get<TransmissionProbabilityOnContact<FP>>() = base_infectiousness;

        return this->get_result().get_last_value();
    }

private:
    FP m_t_last_npi_check;
    std::pair<FP, SimulationTime<FP>> m_dynamic_npi = {-std::numeric_limits<FP>::max(), SimulationTime<FP>(0)};
};

/**
 * @brief Specialization of simulate for SECIRVVS models using Simulation.
 *
 * @tparam FP floating point type, e.g., double.
 * @param[in] t0 start time.
 * @param[in] tmax end time.
 * @param[in] dt time step.
 * @param[in] model SECIRVVS model to simulate.
 * @param[in] integrator_core optional IntegratorCore, uses rk45 if nullptr.
 * 
 * @return Returns the result of the simulation.
 */
template <typename FP>
inline auto simulate(FP t0, FP tmax, FP dt, const Model<FP>& model,
                     std::unique_ptr<OdeIntegratorCore<FP>>&& integrator_core = nullptr)
{
    return mio::simulate<FP, Model<FP>, Simulation<FP>>(t0, tmax, dt, model, std::move(integrator_core));
}

/**
 * @brief Specialization of simulate for SECIRVVS models using the FlowSimulation.
 *
 * @tparam FP floating point type, e.g., double.
 * @param[in] t0 start time.
 * @param[in] tmax end time.
 * @param[in] dt time step.
 * @param[in] model SECIRVVS model to simulate.
 * @param[in] integrator_core optional IntegratorCore, uses rk45 if nullptr.
 * 
 * @return Returns the result of the Flowsimulation.
  */
template <typename FP>
inline auto simulate_flows(FP t0, FP tmax, FP dt, const Model<FP>& model,
                           std::unique_ptr<OdeIntegratorCore<FP>>&& integrator_core = nullptr)
{
    return mio::simulate_flows<FP, Model<FP>, Simulation<FP, mio::FlowSimulation<FP, Model<FP>>>>(
        t0, tmax, dt, model, std::move(integrator_core));
}

//see declaration above.
template <typename FP, class Base>
FP get_infections_relative(const Simulation<FP, Base>& sim, FP /*t*/, const Eigen::Ref<const Eigen::VectorX<FP>>& y)
{
    FP sum_inf = 0;
    for (auto i = AgeGroup(0); i < sim.get_model().parameters.get_num_groups(); ++i) {
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedSymptomsNaive});
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedSymptomsNaiveConfirmed});
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedSymptomsPartialImmunity});
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedSymptomsImprovedImmunity});
        sum_inf +=
            sim.get_model().populations.get_from(y, {i, InfectionState::InfectedSymptomsPartialImmunityConfirmed});
        sum_inf +=
            sim.get_model().populations.get_from(y, {i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed});
    }
    auto inf_rel = sum_inf / sim.get_model().populations.get_total();

    return inf_rel;
}

/**
 * @brief Computes the reproduction number at a given index time of the Model
 * output obtained by the Simulation. Uses the Next-Generation Matrix method
 * with 15 infected compartments per age group:
 * (3 vaccination states × 5 disease stages: E, INS, ISy, ISev, ICr).
 * @param t_idx The index time at which the reproduction number is computed.
 * @param sim The simulation holding the SECIRVVS model.
 * @tparam FP floating point type, e.g., double.
 * @tparam Base simulation type that uses a SECIRVVS compartment model. see Simulation.
 * @returns The computed reproduction number at the provided index time.
 */
template <typename FP, class Base>
IOResult<FP> get_reproduction_number(size_t t_idx, const Simulation<FP, Base>& sim)
{
    using std::sin;

    if (!(t_idx < static_cast<size_t>(sim.get_result().get_num_time_points()))) {
        return mio::failure(mio::StatusCode::OutOfRange, "t_idx is not a valid index for the TimeSeries");
    }

    auto const& params = sim.get_model().parameters;
    auto const& pop    = sim.get_model().populations;
    auto const& y      = sim.get_result().get_value(t_idx);

    const size_t num_agegroups = (size_t)params.get_num_groups();
    // 3 immunity states (Naive, PartialImmunity, ImprovedImmunity) × 5 disease states
    // (Exposed, InfectedNoSymptoms, InfectedSymptoms, InfectedSevere, InfectedCritical)
    const size_t num_infected_states_per_group = 15;
    const size_t total                         = num_infected_states_per_group * num_agegroups;
    const FP pi                                = std::numbers::pi_v<ScalarType>;

    Eigen::MatrixX<FP> F = Eigen::MatrixX<FP>::Zero(total, total);
    Eigen::MatrixX<FP> V = Eigen::MatrixX<FP>::Zero(total, total);

    // Compute independent, aggregated quantities
    FP test_and_trace_required = 0.0;
    FP icu_occupancy           = 0.0;
    for (auto i = AgeGroup(0); i < (AgeGroup)num_agegroups; ++i) {
        test_and_trace_required +=
            (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
            params.template get<TimeInfectedNoSymptoms<FP>>()[i] *
            (y[pop.get_flat_index({i, InfectionState::InfectedNoSymptomsNaive})] +
             y[pop.get_flat_index({i, InfectionState::InfectedNoSymptomsPartialImmunity})] +
             y[pop.get_flat_index({i, InfectionState::InfectedNoSymptomsImprovedImmunity})] +
             y[pop.get_flat_index({i, InfectionState::InfectedNoSymptomsNaiveConfirmed})] +
             y[pop.get_flat_index({i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed})] +
             y[pop.get_flat_index({i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed})]);
        icu_occupancy += y[pop.get_flat_index({i, InfectionState::InfectedCriticalNaive})] +
                         y[pop.get_flat_index({i, InfectionState::InfectedCriticalPartialImmunity})] +
                         y[pop.get_flat_index({i, InfectionState::InfectedCriticalImprovedImmunity})];
    }

    FP season_val =
        (1 + params.template get<Seasonality<FP>>() *
                 sin(pi * ((params.template get<StartDay<FP>>() + sim.get_result().get_time(t_idx)) / 182.5 + 0.5)));

    ContactMatrixGroup<FP> const& contact_matrix = params.template get<ContactPatterns<FP>>();

    // Precompute per-age-group quantities
    Eigen::VectorX<FP> divN(num_agegroups);
    Eigen::VectorX<FP> riskFromInfectedSymptomatic(num_agegroups);
    Eigen::VectorX<FP> riskFromInfectedNoSymptoms(num_agegroups);
    Eigen::MatrixX<FP> cont_freq_eff(num_agegroups, num_agegroups);
    Eigen::MatrixX<FP> dRiskSy_dINS(num_agegroups, num_agegroups);
    Eigen::MatrixX<FP> dRiskNS_dINS(num_agegroups, num_agegroups);

    FP TTC             = params.template get<TestAndTraceCapacity<FP>>();
    FP maxRiskSyFactor = params.template get<TestAndTraceCapacityMaxRiskSymptoms<FP>>();
    FP maxRiskNSFactor = params.template get<TestAndTraceCapacityMaxRiskNoSymptoms<FP>>();
    FP rangeSy         = TTC * (maxRiskSyFactor - 1);
    FP rangeNS         = TTC * (maxRiskNSFactor - 1);

    for (AgeGroup k(0); k < (AgeGroup)num_agegroups; k++) {
        // N_k = total alive population in age group k
        FP Nk = y[pop.get_flat_index({k, InfectionState::SusceptibleNaive})] +
                y[pop.get_flat_index({k, InfectionState::SusceptiblePartialImmunity})] +
                y[pop.get_flat_index({k, InfectionState::SusceptibleImprovedImmunity})] +
                y[pop.get_flat_index({k, InfectionState::ExposedNaive})] +
                y[pop.get_flat_index({k, InfectionState::ExposedPartialImmunity})] +
                y[pop.get_flat_index({k, InfectionState::ExposedImprovedImmunity})] +
                y[pop.get_flat_index({k, InfectionState::InfectedNoSymptomsNaive})] +
                y[pop.get_flat_index({k, InfectionState::InfectedNoSymptomsPartialImmunity})] +
                y[pop.get_flat_index({k, InfectionState::InfectedNoSymptomsImprovedImmunity})] +
                y[pop.get_flat_index({k, InfectionState::InfectedNoSymptomsNaiveConfirmed})] +
                y[pop.get_flat_index({k, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed})] +
                y[pop.get_flat_index({k, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed})] +
                y[pop.get_flat_index({k, InfectionState::InfectedSymptomsNaive})] +
                y[pop.get_flat_index({k, InfectionState::InfectedSymptomsPartialImmunity})] +
                y[pop.get_flat_index({k, InfectionState::InfectedSymptomsImprovedImmunity})] +
                y[pop.get_flat_index({k, InfectionState::InfectedSymptomsNaiveConfirmed})] +
                y[pop.get_flat_index({k, InfectionState::InfectedSymptomsPartialImmunityConfirmed})] +
                y[pop.get_flat_index({k, InfectionState::InfectedSymptomsImprovedImmunityConfirmed})] +
                y[pop.get_flat_index({k, InfectionState::InfectedSevereNaive})] +
                y[pop.get_flat_index({k, InfectionState::InfectedSeverePartialImmunity})] +
                y[pop.get_flat_index({k, InfectionState::InfectedSevereImprovedImmunity})] +
                y[pop.get_flat_index({k, InfectionState::InfectedCriticalNaive})] +
                y[pop.get_flat_index({k, InfectionState::InfectedCriticalPartialImmunity})] +
                y[pop.get_flat_index({k, InfectionState::InfectedCriticalImprovedImmunity})];

        divN[(size_t)k] = (Nk < Limits<FP>::zero_tolerance()) ? FP(0.0) : FP(1.0 / Nk);

        riskFromInfectedSymptomatic[(size_t)k] =
            smoother_cosine<FP>(test_and_trace_required, TTC, TTC * maxRiskSyFactor,
                                params.template get<RiskOfInfectionFromSymptomatic<FP>>()[k],
                                params.template get<MaxRiskOfInfectionFromSymptomatic<FP>>()[k]);

        riskFromInfectedNoSymptoms[(size_t)k] =
            smoother_cosine<FP>(test_and_trace_required, TTC, TTC * maxRiskNSFactor,
                                params.template get<RelativeTransmissionNoSymptoms<FP>>()[k], FP(1.0));

        for (AgeGroup l(0); l < (AgeGroup)num_agegroups; l++) {
            FP dPsi_dINS = (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[l]) /
                           params.template get<TimeInfectedNoSymptoms<FP>>()[l];

            // Derivative of riskFromInfectedSymptomatic[k]
            if (test_and_trace_required < TTC || test_and_trace_required > TTC * maxRiskSyFactor) {
                dRiskSy_dINS((size_t)k, (size_t)l) = 0;
            }
            else {
                dRiskSy_dINS((size_t)k, (size_t)l) = -0.5 * pi *
                                                     (params.template get<MaxRiskOfInfectionFromSymptomatic<FP>>()[k] -
                                                      params.template get<RiskOfInfectionFromSymptomatic<FP>>()[k]) /
                                                     rangeSy * dPsi_dINS *
                                                     sin(pi / rangeSy * (test_and_trace_required - TTC));
            }

            // Derivative of riskFromInfectedNoSymptoms[k]
            if (test_and_trace_required < TTC || test_and_trace_required > TTC * maxRiskNSFactor) {
                dRiskNS_dINS((size_t)k, (size_t)l) = 0;
            }
            else {
                dRiskNS_dINS((size_t)k, (size_t)l) =
                    -0.5 * pi * (FP(1.0) - params.template get<RelativeTransmissionNoSymptoms<FP>>()[k]) / rangeNS *
                    dPsi_dINS * sin(pi / rangeNS * (test_and_trace_required - TTC));
            }
        }

        for (Eigen::Index l = 0; l < (Eigen::Index)num_agegroups; l++) {
            cont_freq_eff(l, (size_t)k) =
                season_val * contact_matrix.get_matrix_at(SimulationTime<FP>(sim.get_result().get_time(t_idx)))(
                                 l, static_cast<Eigen::Index>((size_t)k));
        }
    }

    // Check invertibility via ICU capacity subblock (J matrix for ICr)
    // The 3n×3n J subblock contains the ICr diagonal and ICU coupling
    Eigen::MatrixX<FP> J = Eigen::MatrixX<FP>::Zero(3 * num_agegroups, 3 * num_agegroups);
    for (size_t i = 0; i < num_agegroups; i++) {
        AgeGroup ag(i);
        FP critPerSev = params.template get<CriticalPerSevere<FP>>()[ag];
        // Diagonal entries: 1/T_ICr
        J(i, i)                                         = 1.0 / params.template get<TimeInfectedCritical<FP>>()[ag];
        J(num_agegroups + i, num_agegroups + i)         = 1.0 / params.template get<TimeInfectedCritical<FP>>()[ag];
        J(2 * num_agegroups + i, 2 * num_agegroups + i) = 1.0 / params.template get<TimeInfectedCritical<FP>>()[ag];

        // ICU coupling (if in smoother range)
        if (!(icu_occupancy < 0.9 * params.template get<ICUCapacity<FP>>() ||
              icu_occupancy > params.template get<ICUCapacity<FP>>())) {
            FP dCritPerSev = 5.0 * critPerSev * pi / params.template get<ICUCapacity<FP>>() *
                             sin(pi / (0.1 * params.template get<ICUCapacity<FP>>()) *
                                 (icu_occupancy - 0.9 * params.template get<ICUCapacity<FP>>()));

            FP ISevN  = y[pop.get_flat_index({ag, InfectionState::InfectedSevereNaive})];
            FP ISevPI = y[pop.get_flat_index({ag, InfectionState::InfectedSeverePartialImmunity})];
            FP ISevII = y[pop.get_flat_index({ag, InfectionState::InfectedSevereImprovedImmunity})];
            FP T_ISev = params.template get<TimeInfectedSevere<FP>>()[ag];

            for (size_t j = 0; j < num_agegroups; j++) {
                // ICr_N rows: coupling to all ICr columns
                J(i, j) -= ISevN / T_ISev * dCritPerSev;
                J(i, num_agegroups + j) -= ISevN / T_ISev * dCritPerSev;
                J(i, 2 * num_agegroups + j) -= ISevN / T_ISev * dCritPerSev;
                // ICr_PI rows
                J(num_agegroups + i, j) -= ISevPI / T_ISev * dCritPerSev;
                J(num_agegroups + i, num_agegroups + j) -= ISevPI / T_ISev * dCritPerSev;
                J(num_agegroups + i, 2 * num_agegroups + j) -= ISevPI / T_ISev * dCritPerSev;
                // ICr_II rows
                J(2 * num_agegroups + i, j) -= ISevII / T_ISev * dCritPerSev;
                J(2 * num_agegroups + i, num_agegroups + j) -= ISevII / T_ISev * dCritPerSev;
                J(2 * num_agegroups + i, 2 * num_agegroups + j) -= ISevII / T_ISev * dCritPerSev;
            }
        }
    }

    if (J.determinant() == 0) {
        return mio::failure(mio::StatusCode::UnknownError, "Matrix V is not invertible");
    }

    // ---- Build F matrix ----
    // Layout per immunity state (v=0: Naive, v=1: PI, v=2: II):
    //   block v: [E_v, INS_v, ISy_v, ISev_v, ICr_v], each of size num_agegroups
    //   offset = v * 5 * num_agegroups
    // Only E rows (offset + 0..num_agegroups-1) have nonzero entries in F.
    // Only INS columns (offset + num_agegroups..2n-1) and ISy columns (offset + 2n..3n-1) are nonzero.
    for (size_t i = 0; i < num_agegroups; i++) {
        AgeGroup ag_i(i);
        FP S_N  = y[pop.get_flat_index({ag_i, InfectionState::SusceptibleNaive})];
        FP S_PI = y[pop.get_flat_index({ag_i, InfectionState::SusceptiblePartialImmunity})];
        FP S_II = y[pop.get_flat_index({ag_i, InfectionState::SusceptibleImprovedImmunity})];

        FP rho_i = params.template get<TransmissionProbabilityOnContact<FP>>()[ag_i];
        // coefficients: S_v * rho_v for each immunity state
        FP coeff_N   = S_N * rho_i;
        FP coeff_PI  = S_PI * rho_i * params.template get<ReducExposedPartialImmunity<FP>>()[ag_i];
        FP coeff_II  = S_II * rho_i * params.template get<ReducExposedImprovedImmunity<FP>>()[ag_i];
        FP coeffs[3] = {coeff_N, coeff_PI, coeff_II};

        for (size_t j = 0; j < num_agegroups; j++) {
            // Compute derivative contribution (temp) for B_{i,j}
            // temp = \sum_k \phi_{ik}/N_k * (ISy_k * dRiskSy(k,j) + INS_k * dRiskNS(k,j))
            FP temp = 0;
            for (size_t k = 0; k < num_agegroups; k++) {
                FP INS_k = y[pop.get_flat_index({AgeGroup(k), InfectionState::InfectedNoSymptomsNaive})] +
                           y[pop.get_flat_index({AgeGroup(k), InfectionState::InfectedNoSymptomsPartialImmunity})] +
                           y[pop.get_flat_index({AgeGroup(k), InfectionState::InfectedNoSymptomsImprovedImmunity})];
                FP ISy_k = y[pop.get_flat_index({AgeGroup(k), InfectionState::InfectedSymptomsNaive})] +
                           y[pop.get_flat_index({AgeGroup(k), InfectionState::InfectedSymptomsPartialImmunity})] +
                           y[pop.get_flat_index({AgeGroup(k), InfectionState::InfectedSymptomsImprovedImmunity})];
                temp += cont_freq_eff(i, k) * divN[k] * (ISy_k * dRiskSy_dINS(k, j) + INS_k * dRiskNS_dINS(k, j));
            }

            FP B_ij = cont_freq_eff(i, j) * riskFromInfectedNoSymptoms[j] * divN[j] + temp;
            FP C_ij = cont_freq_eff(i, j) * riskFromInfectedSymptomatic[j] * divN[j];

            // Fill F for all 3 immunity states (rows) × all 3 immunity states (columns)
            for (size_t v = 0; v < 3; v++) {
                size_t row = v * 5 * num_agegroups + i; // E row for immunity state v
                for (size_t w = 0; w < 3; w++) {
                    size_t col_base                          = w * 5 * num_agegroups;
                    F(row, col_base + num_agegroups + j)     = coeffs[v] * B_ij; // INS column
                    F(row, col_base + 2 * num_agegroups + j) = coeffs[v] * C_ij; // ISy column
                }
            }
        }
    }

    // ---- Build V matrix ----
    for (Eigen::Index i = 0; i < (Eigen::Index)num_agegroups; i++) {
        AgeGroup ag(i);

        FP T_E    = params.template get<TimeExposed<FP>>()[ag];
        FP T_INS  = params.template get<TimeInfectedNoSymptoms<FP>>()[ag];
        FP T_ISy  = params.template get<TimeInfectedSymptoms<FP>>()[ag];
        FP T_ISev = params.template get<TimeInfectedSevere<FP>>()[ag];

        FP recPerINS = params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[ag];
        FP sevPerISy = params.template get<SeverePerInfectedSymptoms<FP>>()[ag];

        FP critPerSevAdj = smoother_cosine<FP>(icu_occupancy, 0.90 * params.template get<ICUCapacity<FP>>(),
                                               params.template get<ICUCapacity<FP>>(),
                                               params.template get<CriticalPerSevere<FP>>()[ag], FP(0));

        FP reducExpPI = params.template get<ReducExposedPartialImmunity<FP>>()[ag];
        FP reducExpII = params.template get<ReducExposedImprovedImmunity<FP>>()[ag];
        FP reducSymPI = params.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[ag];
        FP reducSymII = params.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[ag];
        FP reducSevPI = params.template get<ReducInfectedSevereCriticalDeadPartialImmunity<FP>>()[ag];
        FP reducSevII = params.template get<ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>()[ag];
        FP kappa      = params.template get<ReducTimeInfectedMild<FP>>()[ag];

        // Naive block (offset 0)
        {
            size_t off = 0;
            // E_N
            V(off + i, off + i) = 1.0 / T_E;
            // INS_N
            V(off + num_agegroups + i, off + i)                 = -1.0 / T_E;
            V(off + num_agegroups + i, off + num_agegroups + i) = 1.0 / T_INS;
            // ISy_N
            V(off + 2 * num_agegroups + i, off + num_agegroups + i)     = -(1 - recPerINS) / T_INS;
            V(off + 2 * num_agegroups + i, off + 2 * num_agegroups + i) = 1.0 / T_ISy;
            // ISev_N
            V(off + 3 * num_agegroups + i, off + 2 * num_agegroups + i) = -sevPerISy / T_ISy;
            V(off + 3 * num_agegroups + i, off + 3 * num_agegroups + i) = 1.0 / T_ISev;
            // ICr_N -> ISev_N coupling
            V(off + 4 * num_agegroups + i, off + 3 * num_agegroups + i) = -critPerSevAdj / T_ISev;
        }

        // Partial Immunity block (offset 5n)
        {
            size_t off = 5 * num_agegroups;
            // E_PI
            V(off + i, off + i) = 1.0 / T_E;
            // INS_PI (kappa = reducTimeInfectedMild)
            V(off + num_agegroups + i, off + i)                 = -1.0 / T_E;
            V(off + num_agegroups + i, off + num_agegroups + i) = 1.0 / (T_INS * kappa);
            // ISy_PI: Omega4 = (reducSymPI/reducExpPI) * (1 - recPerINS)
            FP omega4                                                   = (reducSymPI / reducExpPI) * (1 - recPerINS);
            V(off + 2 * num_agegroups + i, off + num_agegroups + i)     = -omega4 / (T_INS * kappa);
            V(off + 2 * num_agegroups + i, off + 2 * num_agegroups + i) = 1.0 / (T_ISy * kappa);
            // ISev_PI: Omega5 = (reducSevPI/reducSymPI) * sevPerISy
            FP omega5                                                   = (reducSevPI / reducSymPI) * sevPerISy;
            V(off + 3 * num_agegroups + i, off + 2 * num_agegroups + i) = -omega5 / (T_ISy * kappa);
            V(off + 3 * num_agegroups + i, off + 3 * num_agegroups + i) = 1.0 / T_ISev;
            // ICr_PI
            V(off + 4 * num_agegroups + i, off + 3 * num_agegroups + i) = -critPerSevAdj / T_ISev;
        }

        // Improved Immunity block (offset 10n)
        {
            size_t off = 10 * num_agegroups;
            // E_II
            V(off + i, off + i) = 1.0 / T_E;
            // INS_II
            V(off + num_agegroups + i, off + i)                 = -1.0 / T_E;
            V(off + num_agegroups + i, off + num_agegroups + i) = 1.0 / (T_INS * kappa);
            // ISy_II: Omega7 = (reducSymII/reducExpII) * (1 - recPerINS)
            FP omega7                                                   = (reducSymII / reducExpII) * (1 - recPerINS);
            V(off + 2 * num_agegroups + i, off + num_agegroups + i)     = -omega7 / (T_INS * kappa);
            V(off + 2 * num_agegroups + i, off + 2 * num_agegroups + i) = 1.0 / (T_ISy * kappa);
            // ISev_II: Omega8 = (reducSevII/reducSymII) * sevPerISy
            FP omega8                                                   = (reducSevII / reducSymII) * sevPerISy;
            V(off + 3 * num_agegroups + i, off + 2 * num_agegroups + i) = -omega8 / (T_ISy * kappa);
            V(off + 3 * num_agegroups + i, off + 3 * num_agegroups + i) = 1.0 / T_ISev;
            // ICr_II
            V(off + 4 * num_agegroups + i, off + 3 * num_agegroups + i) = -critPerSevAdj / T_ISev;
        }

        // ICr diagonal (J submatrix entries)
        // Fill from the precomputed J matrix: maps local J indices to global V indices
        // J layout: [ICr_N(0..num_agegroups-1), ICr_PI(num_agegroups..2n-1), ICr_II(2n..3n-1)]
        // V layout: ICr_N at 4n, ICr_PI at 9n, ICr_II at 14n
        size_t icr_offsets[3] = {4 * num_agegroups, 9 * num_agegroups, 14 * num_agegroups};
        for (size_t v = 0; v < 3; v++) {
            for (size_t w = 0; w < 3; w++) {
                for (size_t j = 0; j < num_agegroups; j++) {
                    V(icr_offsets[v] + i, icr_offsets[w] + j) = J(v * num_agegroups + i, w * num_agegroups + j);
                }
            }
        }
    }

    // ---- Compute Next-Generation Matrix and spectral radius ----
    V = V.inverse();

    Eigen::MatrixX<FP> NextGenMatrix(total, total);
    NextGenMatrix.noalias() = F * V;

    Eigen::ComplexEigenSolver<Eigen::MatrixX<FP>> ces;
    ces.compute(NextGenMatrix);
    FP rho = ces.eigenvalues().cwiseAbs().maxCoeff();

    return mio::success(rho);
}

/**
 * @brief Computes the reproduction number for all time points of the Model output obtained by the Simulation.
 * @param sim The Model Simulation.
 * @tparam FP floating point type, e.g., double.
 * @tparam Base simulation type that uses a SECIRVVS compartment model. see Simulation.
 * @returns Eigen::Vector containing all reproduction numbers.
 */
template <typename FP, class Base>
Eigen::VectorX<FP> get_reproduction_numbers(const Simulation<FP, Base>& sim)
{
    Eigen::VectorX<FP> temp(sim.get_result().get_num_time_points());
    for (int i = 0; i < sim.get_result().get_num_time_points(); i++) {
        temp[i] = get_reproduction_number<FP>((size_t)i, sim).value();
    }
    return temp;
}

/**
 * @brief Computes the reproduction number at a given time point of the Simulation.
 * If the particular time point is not part of the output, a linearly interpolated value is returned.
 * @param t_value The time point at which the reproduction number should be computed.
 * @param sim The Model Simulation.
 * @tparam FP floating point type, e.g., double.
 * @tparam Base simulation type that uses a SECIRVVS compartment model. see Simulation.
 * @returns The computed reproduction number at the provided time point, potentially using linear interpolation.
 */
template <typename FP, class Base>
IOResult<FP> get_reproduction_number(FP t_value, const Simulation<FP, Base>& sim)
{
    if (t_value < sim.get_result().get_time(0) || t_value > sim.get_result().get_last_time()) {
        return mio::failure(mio::StatusCode::OutOfRange,
                            "Cannot interpolate reproduction number outside computed horizon of the TimeSeries");
    }

    if (t_value == sim.get_result().get_time(0)) {
        return mio::success(get_reproduction_number<FP>((size_t)0, sim).value());
    }

    const auto& times = sim.get_result().get_times();
    auto time_late    = std::distance(times.begin(), std::lower_bound(times.begin(), times.end(), t_value));

    FP y1 = get_reproduction_number<FP>(static_cast<size_t>(time_late - 1), sim).value();
    FP y2 = get_reproduction_number<FP>(static_cast<size_t>(time_late), sim).value();

    auto result = linear_interpolation(t_value, sim.get_result().get_time(time_late - 1),
                                       sim.get_result().get_time(time_late), y1, y2);
    return mio::success(static_cast<FP>(result));
}

/**
 * Get mobility factors.
 * Used by mobility graph simulation.
 * Like infection risk, mobility of infected individuals is reduced if they are well isolated.
 * @param model the compartment model with initial values.
 * @param t current simulation time.
 * @param y current value of compartments.
 * @return vector expression, same size as y, with mobility factors per compartment.
 * @tparam Base simulation type that uses a secir compartment model. see Simulation.
 */
template <typename FP, class Base = mio::Simulation<FP, Model<FP>>>
auto get_mobility_factors(const Simulation<Base>& sim, FP /*t*/, const Eigen::Ref<const Eigen::VectorX<FP>>& y)

{
    auto& params = sim.get_model().parameters;
    //parameters as arrays
    auto& p_asymp   = params.template get<RecoveredPerInfectedNoSymptoms<FP>>().array().template cast<FP>();
    auto& p_inf     = params.template get<RiskOfInfectionFromSymptomatic<FP>>().array().template cast<FP>();
    auto& p_inf_max = params.template get<MaxRiskOfInfectionFromSymptomatic<FP>>().array().template cast<FP>();
    //slice of InfectedNoSymptoms
    auto y_INS = slice(y, {Eigen::Index(InfectionState::InfectedNoSymptomsNaive),
                           Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)}) +
                 slice(y, {Eigen::Index(InfectionState::InfectedNoSymptomsPartialImmunity),
                           Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)}) +
                 slice(y, {Eigen::Index(InfectionState::InfectedNoSymptomsImprovedImmunity),
                           Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)});

    //compute isolation, same as infection risk from main model
    auto test_and_trace_required =
        ((1 - p_asymp) / params.template get<TimeInfectedNoSymptoms<FP>>().array().template cast<FP>() * y_INS.array())
            .sum();
    auto riskFromInfectedSymptomatic =
        smoother_cosine<FP>(test_and_trace_required, FP(params.template get<TestAndTraceCapacity<FP>>()),
                            params.template get<TestAndTraceCapacity<FP>>() *
                                params.template get<TestAndTraceCapacityMaxRiskSymptoms<FP>>(),
                            p_inf.matrix(), p_inf_max.matrix());

    //set factor for infected
    auto factors = Eigen::VectorX<FP>::Ones(y.rows()).eval();
    slice(factors, {Eigen::Index(InfectionState::InfectedSymptomsNaive), Eigen::Index(size_t(params.get_num_groups())),
                    Eigen::Index(InfectionState::Count)})
        .array() = riskFromInfectedSymptomatic;
    slice(factors, {Eigen::Index(InfectionState::InfectedSymptomsPartialImmunity),
                    Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)})
        .array() = riskFromInfectedSymptomatic;
    slice(factors, {Eigen::Index(InfectionState::InfectedSymptomsImprovedImmunity),
                    Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)})
        .array() = riskFromInfectedSymptomatic;
    return factors;
}

template <typename FP, class Base = mio::Simulation<FP, Model<FP>>>
auto test_commuters(Simulation<FP, Base>& sim, Eigen::Ref<Eigen::VectorX<FP>> mobile_population, FP time)
{
    auto& model     = sim.get_model();
    FP nondetection = 1.0;
    if (time >= model.parameters.get_start_commuter_detection() &&
        time < model.parameters.get_end_commuter_detection()) {
        nondetection = (FP)model.parameters.get_commuter_nondetection();
    }
    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); ++i) {
        auto ISyNi  = model.populations.get_flat_index({i, InfectionState::InfectedSymptomsNaive});
        auto ISyNCi = model.populations.get_flat_index({i, InfectionState::InfectedSymptomsNaiveConfirmed});
        auto INSNi  = model.populations.get_flat_index({i, InfectionState::InfectedNoSymptomsNaive});
        auto INSNCi = model.populations.get_flat_index({i, InfectionState::InfectedNoSymptomsNaiveConfirmed});

        auto ISPIi  = model.populations.get_flat_index({i, InfectionState::InfectedSymptomsPartialImmunity});
        auto ISPICi = model.populations.get_flat_index({i, InfectionState::InfectedSymptomsPartialImmunityConfirmed});
        auto INSPIi = model.populations.get_flat_index({i, InfectionState::InfectedNoSymptomsPartialImmunity});
        auto INSPICi =
            model.populations.get_flat_index({i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed});

        auto ISyIIi  = model.populations.get_flat_index({i, InfectionState::InfectedSymptomsImprovedImmunity});
        auto ISyIICi = model.populations.get_flat_index({i, InfectionState::InfectedSymptomsImprovedImmunityConfirmed});
        auto INSIIi  = model.populations.get_flat_index({i, InfectionState::InfectedNoSymptomsImprovedImmunity});
        auto INSIICi =
            model.populations.get_flat_index({i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed});

        //put detected commuters in their own compartment so they don't contribute to infections in their home node
        sim.get_result().get_last_value()[ISyNi] -= mobile_population[ISyNi] * (1 - nondetection);
        sim.get_result().get_last_value()[ISyNCi] += mobile_population[ISyNi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSNi] -= mobile_population[INSNi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSNCi] += mobile_population[INSNi] * (1 - nondetection);

        sim.get_result().get_last_value()[ISPIi] -= mobile_population[ISPIi] * (1 - nondetection);
        sim.get_result().get_last_value()[ISPICi] += mobile_population[ISPIi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSPIi] -= mobile_population[INSPIi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSPICi] += mobile_population[INSPIi] * (1 - nondetection);

        sim.get_result().get_last_value()[ISyIIi] -= mobile_population[ISyIIi] * (1 - nondetection);
        sim.get_result().get_last_value()[ISyIICi] += mobile_population[ISyIIi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSIIi] -= mobile_population[INSIIi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSIICi] += mobile_population[INSIIi] * (1 - nondetection);

        //reduce the number of commuters
        mobile_population[ISyNi] *= nondetection;
        mobile_population[INSNi] *= nondetection;

        mobile_population[ISPIi] *= nondetection;
        mobile_population[INSPIi] *= nondetection;

        mobile_population[ISyIIi] *= nondetection;
        mobile_population[INSIIi] *= nondetection;
    }
}

} // namespace osecirvvs
} // namespace mio

#endif //MIO_ODE_SECIRVVS_MODEL_H
