/*
* * Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker, Wadim Koslow, Daniel Abele, Martin J. Kühn
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
#ifndef MIO_ODE_SECIRTS_MODEL_H
#define MIO_ODE_SECIRTS_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/epidemiology/populations.h"
#include "ode_secirts/infection_state.h"
#include "ode_secirts/parameters.h"
#include "memilio/math/smoother.h"
#include "memilio/math/eigen_util.h"

#include <numbers>

namespace mio
{
namespace osecirts
{
// clang-format off
using Flows = TypeList<
    //naive
    Flow<InfectionState::SusceptibleNaive,                            InfectionState::ExposedNaive>,
    Flow<InfectionState::SusceptibleNaive,                            InfectionState::TemporaryImmunePartialImmunity>,
    Flow<InfectionState::ExposedNaive,                                InfectionState::InfectedNoSymptomsNaive>,
    Flow<InfectionState::InfectedNoSymptomsNaive,                     InfectionState::InfectedSymptomsNaive>,
    Flow<InfectionState::InfectedNoSymptomsNaive,                     InfectionState::TemporaryImmunePartialImmunity>,
    Flow<InfectionState::InfectedNoSymptomsNaiveConfirmed,            InfectionState::InfectedSymptomsNaiveConfirmed>,
    Flow<InfectionState::InfectedNoSymptomsNaiveConfirmed,            InfectionState::TemporaryImmunePartialImmunity>,
    Flow<InfectionState::InfectedSymptomsNaive,                       InfectionState::InfectedSevereNaive>,
    Flow<InfectionState::InfectedSymptomsNaive,                       InfectionState::TemporaryImmunePartialImmunity>,
    Flow<InfectionState::InfectedSymptomsNaiveConfirmed,              InfectionState::InfectedSevereNaive>,
    Flow<InfectionState::InfectedSymptomsNaiveConfirmed,              InfectionState::TemporaryImmunePartialImmunity>,
    Flow<InfectionState::InfectedSevereNaive,                         InfectionState::InfectedCriticalNaive>,
    Flow<InfectionState::InfectedSevereNaive,                         InfectionState::TemporaryImmunePartialImmunity>,
    Flow<InfectionState::InfectedSevereNaive,                         InfectionState::DeadNaive>,
    Flow<InfectionState::InfectedCriticalNaive,                       InfectionState::DeadNaive>,
    Flow<InfectionState::InfectedCriticalNaive,                       InfectionState::TemporaryImmunePartialImmunity>,
    Flow<InfectionState::TemporaryImmunePartialImmunity,              InfectionState::SusceptiblePartialImmunity>,
    //partial immunity
    Flow<InfectionState::SusceptiblePartialImmunity,                  InfectionState::ExposedPartialImmunity>,
    Flow<InfectionState::SusceptiblePartialImmunity,                  InfectionState::TemporaryImmuneImprovedImmunity>,
    Flow<InfectionState::ExposedPartialImmunity,                      InfectionState::InfectedNoSymptomsPartialImmunity>,
    Flow<InfectionState::InfectedNoSymptomsPartialImmunity,           InfectionState::InfectedSymptomsPartialImmunity>,
    Flow<InfectionState::InfectedNoSymptomsPartialImmunity,           InfectionState::TemporaryImmuneImprovedImmunity>,
    Flow<InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,  InfectionState::InfectedSymptomsPartialImmunityConfirmed>,
    Flow<InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,  InfectionState::TemporaryImmuneImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsPartialImmunity,             InfectionState::InfectedSeverePartialImmunity>,
    Flow<InfectionState::InfectedSymptomsPartialImmunity,             InfectionState::TemporaryImmuneImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsPartialImmunityConfirmed,    InfectionState::InfectedSeverePartialImmunity>,
    Flow<InfectionState::InfectedSymptomsPartialImmunityConfirmed,    InfectionState::TemporaryImmuneImprovedImmunity>,
    Flow<InfectionState::InfectedSeverePartialImmunity,               InfectionState::InfectedCriticalPartialImmunity>,
    Flow<InfectionState::InfectedSeverePartialImmunity,               InfectionState::TemporaryImmuneImprovedImmunity>,
    Flow<InfectionState::InfectedSeverePartialImmunity,               InfectionState::DeadPartialImmunity>,
    Flow<InfectionState::InfectedCriticalPartialImmunity,             InfectionState::DeadPartialImmunity>,
    Flow<InfectionState::InfectedCriticalPartialImmunity,             InfectionState::TemporaryImmuneImprovedImmunity>,
    //improved immunity
    Flow<InfectionState::SusceptibleImprovedImmunity,                 InfectionState::ExposedImprovedImmunity>,
    Flow<InfectionState::SusceptibleImprovedImmunity,                 InfectionState::TemporaryImmuneImprovedImmunity>,
    Flow<InfectionState::ExposedImprovedImmunity,                     InfectionState::InfectedNoSymptomsImprovedImmunity>,
    Flow<InfectionState::InfectedNoSymptomsImprovedImmunity,          InfectionState::InfectedSymptomsImprovedImmunity>,
    Flow<InfectionState::InfectedNoSymptomsImprovedImmunity,          InfectionState::TemporaryImmuneImprovedImmunity>,
    Flow<InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed, InfectionState::InfectedSymptomsImprovedImmunityConfirmed>,
    Flow<InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed, InfectionState::TemporaryImmuneImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsImprovedImmunity,            InfectionState::InfectedSevereImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsImprovedImmunity,            InfectionState::TemporaryImmuneImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsImprovedImmunityConfirmed,   InfectionState::InfectedSevereImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsImprovedImmunityConfirmed,   InfectionState::TemporaryImmuneImprovedImmunity>,
    Flow<InfectionState::InfectedSevereImprovedImmunity,              InfectionState::InfectedCriticalImprovedImmunity>,
    Flow<InfectionState::InfectedSevereImprovedImmunity,              InfectionState::TemporaryImmuneImprovedImmunity>,
    Flow<InfectionState::InfectedSevereImprovedImmunity,              InfectionState::DeadImprovedImmunity>,
    Flow<InfectionState::InfectedCriticalImprovedImmunity,            InfectionState::DeadImprovedImmunity>,
    Flow<InfectionState::InfectedCriticalImprovedImmunity,            InfectionState::TemporaryImmuneImprovedImmunity>,

    // waning
    Flow<InfectionState::TemporaryImmunePartialImmunity,              InfectionState::SusceptiblePartialImmunity>,
    Flow<InfectionState::TemporaryImmuneImprovedImmunity,             InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::SusceptibleImprovedImmunity,                 InfectionState::SusceptiblePartialImmunity>,
    Flow<InfectionState::SusceptiblePartialImmunity,                  InfectionState::SusceptibleNaive>>;
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
        using std::floor;
        using std::min;

        auto const& params   = this->parameters;
        AgeGroup n_agegroups = params.get_num_groups();

        ContactMatrixGroup<FP> const& contact_matrix = params.template get<ContactPatterns<FP>>();

        FP icu_occupancy           = 0.0;
        FP test_and_trace_required = 0.0;
        for (auto i = AgeGroup(0); i < n_agegroups; ++i) {
            // naive flow to symptomatic in unit time
            test_and_trace_required +=
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                params.template get<TimeInfectedNoSymptoms<FP>>()[i] *
                (this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptomsNaive}) +
                 this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptomsNaiveConfirmed}));
            // partial immunity flow to symptomatic in unit time
            test_and_trace_required +=
                (params.template get<ReducInfectedSymptomsPartialImmunity<FP>>()[i] /
                 params.template get<ReducExposedPartialImmunity<FP>>()[i]) *
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedNoSymptoms<FP>>()[i] *
                 params.template get<ReducTimeInfectedMild<FP>>()[i]) *
                (this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptomsPartialImmunity}) +
                 this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}));
            // improved immunity flow to symptomatic in unit time
            test_and_trace_required +=
                (params.template get<ReducInfectedSymptomsImprovedImmunity<FP>>()[i] /
                 params.template get<ReducExposedImprovedImmunity<FP>>()[i]) *
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedNoSymptoms<FP>>()[i] *
                 params.template get<ReducTimeInfectedMild<FP>>()[i]) *
                (this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptomsImprovedImmunity}) +
                 this->populations.get_from(pop, {i, InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}));
            icu_occupancy += this->populations.get_from(pop, {i, InfectionState::InfectedCriticalNaive}) +
                             this->populations.get_from(pop, {i, InfectionState::InfectedCriticalPartialImmunity}) +
                             this->populations.get_from(pop, {i, InfectionState::InfectedCriticalImprovedImmunity});
        }

        // get vaccinations
        auto const partial_vaccination = vaccinations_at(t, params.template get<DailyPartialVaccinations<FP>>());
        auto const full_vaccination    = vaccinations_at(t, params.template get<DailyFullVaccinations<FP>>());
        auto const booster_vaccination = vaccinations_at(t, params.template get<DailyBoosterVaccinations<FP>>());

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
            size_t TImm1 = this->populations.get_flat_index({i, InfectionState::TemporaryImmunePartialImmunity});
            size_t TImm2 = this->populations.get_flat_index({i, InfectionState::TemporaryImmuneImprovedImmunity});

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
                // std::fmod('time', 365.0) is non differentiable. Use std::floor instead to normalize 'time'.
                FP normalized_time = (params.template get<StartDay<FP>>() + t) -
                                     365.0 * floor((params.template get<StartDay<FP>>() + t) / 365.0);
                FP season_val = (1 + params.template get<Seasonality<FP>>() *
                                         sin(std::numbers::pi_v<ScalarType> * (normalized_time / 182.5 + 0.5)));

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

            // vaccinations
            flows[this->template get_flat_flow_index<InfectionState::SusceptibleNaive,
                                                     InfectionState::TemporaryImmunePartialImmunity>({i})] =
                min<FP>(y[SNi] - flows[this->template get_flat_flow_index<InfectionState::SusceptibleNaive,
                                                                          InfectionState::ExposedNaive>({i})],
                        partial_vaccination[static_cast<size_t>(i)]);

            flows[this->template get_flat_flow_index<InfectionState::SusceptiblePartialImmunity,
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
                min<FP>(y[SPIi] -
                            flows[this->template get_flat_flow_index<InfectionState::SusceptiblePartialImmunity,
                                                                     InfectionState::ExposedPartialImmunity>({i})],
                        full_vaccination[static_cast<size_t>(i)]);

            flows[this->template get_flat_flow_index<InfectionState::SusceptibleImprovedImmunity,
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
                min<FP>(y[SIIi] -
                            flows[this->template get_flat_flow_index<InfectionState::SusceptibleImprovedImmunity,
                                                                     InfectionState::ExposedImprovedImmunity>({i})],
                        booster_vaccination[static_cast<size_t>(i)]);

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
                                                     InfectionState::TemporaryImmunePartialImmunity>({i})] =
                params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i] *
                (1 / params.template get<TimeInfectedNoSymptoms<FP>>()[i]) * y[INSNi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsNaive,
                                                     InfectionState::InfectedSymptomsNaive>({i})] =
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                params.template get<TimeInfectedNoSymptoms<FP>>()[i] * y[INSNi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsNaiveConfirmed,
                                                     InfectionState::InfectedSymptomsNaiveConfirmed>({i})] =
                (1 - params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i]) /
                params.template get<TimeInfectedNoSymptoms<FP>>()[i] * y[INSNCi];
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsNaiveConfirmed,
                                                     InfectionState::TemporaryImmunePartialImmunity>({i})] =
                params.template get<RecoveredPerInfectedNoSymptoms<FP>>()[i] *
                (1 / params.template get<TimeInfectedNoSymptoms<FP>>()[i]) * y[INSNCi];

            // // InfectedSymptoms
            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsNaive,
                                                     InfectionState::InfectedSevereNaive>({i})] =
                params.template get<SeverePerInfectedSymptoms<FP>>()[i] /
                params.template get<TimeInfectedSymptoms<FP>>()[i] * y[ISyNi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsNaive,
                                                     InfectionState::TemporaryImmunePartialImmunity>({i})] =
                (1 - params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                params.template get<TimeInfectedSymptoms<FP>>()[i] * y[ISyNi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsNaiveConfirmed,
                                                     InfectionState::InfectedSevereNaive>({i})] =
                params.template get<SeverePerInfectedSymptoms<FP>>()[i] /
                params.template get<TimeInfectedSymptoms<FP>>()[i] * y[ISyNCi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsNaiveConfirmed,
                                                     InfectionState::TemporaryImmunePartialImmunity>({i})] =
                (1 - params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                params.template get<TimeInfectedSymptoms<FP>>()[i] * y[ISyNCi];

            // InfectedSevere
            flows[this->template get_flat_flow_index<InfectionState::InfectedSevereNaive,
                                                     InfectionState::InfectedCriticalNaive>({i})] =
                criticalPerSevereAdjusted / params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevNi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSevereNaive,
                                                     InfectionState::TemporaryImmunePartialImmunity>({i})] =
                (1 - params.template get<CriticalPerSevere<FP>>()[i]) /
                params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevNi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSevereNaive, InfectionState::DeadNaive>(
                {i})] = deathsPerSevereAdjusted / params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevNi];
            // InfectedCritical
            flows[this->template get_flat_flow_index<InfectionState::InfectedCriticalNaive, InfectionState::DeadNaive>(
                {i})] = params.template get<DeathsPerCritical<FP>>()[i] /
                        params.template get<TimeInfectedCritical<FP>>()[i] * y[ICrNi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedCriticalNaive,
                                                     InfectionState::TemporaryImmunePartialImmunity>({i})] =
                (1 - params.template get<DeathsPerCritical<FP>>()[i]) /
                params.template get<TimeInfectedCritical<FP>>()[i] * y[ICrNi];

            // Waning immunity
            flows[this->template get_flat_flow_index<InfectionState::SusceptiblePartialImmunity,
                                                     InfectionState::SusceptibleNaive>({i})] =
                1 / params.template get<TimeWaningPartialImmunity<FP>>()[i] * y[SPIi];
            flows[this->template get_flat_flow_index<InfectionState::TemporaryImmunePartialImmunity,
                                                     InfectionState::SusceptiblePartialImmunity>({i})] =
                1 / params.template get<TimeTemporaryImmunityPI<FP>>()[i] * y[TImm1];

            // /**** path of partially immune ***/

            // Exposed
            flows[this->template get_flat_flow_index<InfectionState::ExposedPartialImmunity,
                                                     InfectionState::InfectedNoSymptomsPartialImmunity>({i})] +=
                y[EPIi] / params.template get<TimeExposed<FP>>()[i];

            // InfectedNoSymptoms
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsPartialImmunity,
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
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
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
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
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity) *
                         params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyPIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsPartialImmunityConfirmed,
                                                     InfectionState::InfectedSeverePartialImmunity>({i})] =
                reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity *
                params.template get<SeverePerInfectedSymptoms<FP>>()[i] /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyPICi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsPartialImmunityConfirmed,
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity) *
                         params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyPICi];

            // InfectedSevere
            flows[this->template get_flat_flow_index<InfectionState::InfectedSeverePartialImmunity,
                                                     InfectionState::InfectedCriticalPartialImmunity>({i})] =
                reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity *
                criticalPerSevereAdjusted / params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevPIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSeverePartialImmunity,
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
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
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                         params.template get<DeathsPerCritical<FP>>()[i]) /
                params.template get<TimeInfectedCritical<FP>>()[i] * y[ICrPIi];

            // /**** path of improved immunity ***/
            // Exposed
            flows[this->template get_flat_flow_index<InfectionState::ExposedImprovedImmunity,
                                                     InfectionState::InfectedNoSymptomsImprovedImmunity>({i})] +=
                y[EIIi] / params.template get<TimeExposed<FP>>()[i];

            // InfectedNoSymptoms
            flows[this->template get_flat_flow_index<InfectionState::InfectedNoSymptomsImprovedImmunity,
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
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
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
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
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
                (1 - (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity) *
                         params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyIIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsImprovedImmunityConfirmed,
                                                     InfectionState::InfectedSevereImprovedImmunity>({i})] =
                reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity *
                params.template get<SeverePerInfectedSymptoms<FP>>()[i] /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyIICi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSymptomsImprovedImmunityConfirmed,
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
                (1 - (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity) *
                         params.template get<SeverePerInfectedSymptoms<FP>>()[i]) /
                (params.template get<TimeInfectedSymptoms<FP>>()[i] * reducTimeInfectedMild) * y[ISyIICi];

            // InfectedSevere
            flows[this->template get_flat_flow_index<InfectionState::InfectedSevereImprovedImmunity,
                                                     InfectionState::InfectedCriticalImprovedImmunity>({i})] =
                reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity *
                criticalPerSevereAdjusted / params.template get<TimeInfectedSevere<FP>>()[i] * y[ISevIIi];

            flows[this->template get_flat_flow_index<InfectionState::InfectedSevereImprovedImmunity,
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
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
                                                     InfectionState::TemporaryImmuneImprovedImmunity>({i})] =
                (1 -
                 (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                     params.template get<DeathsPerCritical<FP>>()[i]) /
                params.template get<TimeInfectedCritical<FP>>()[i] * y[ICrIIi];

            // Waning immunity
            flows[this->template get_flat_flow_index<InfectionState::TemporaryImmuneImprovedImmunity,
                                                     InfectionState::SusceptibleImprovedImmunity>({i})] =
                1 / params.template get<TimeTemporaryImmunityII<FP>>()[i] * y[TImm2];

            flows[this->template get_flat_flow_index<InfectionState::SusceptibleImprovedImmunity,
                                                     InfectionState::SusceptiblePartialImmunity>({i})] =
                1 / params.template get<TimeWaningImprovedImmunity<FP>>()[i] * y[SIIi];
        }
    }

    /**
    * @brief Calculates smoothed vaccinations for a given time point.
    *
    * This function calculates the number of vaccinations for each age group at a given time t,
    * based on daily vaccinations data. The smoothing is done using a cosine function.
    *
    * @param t The time in the simulation.
    * @param daily_vaccinations The daily vaccinations data, indexed by age group and simulation day.
    * @param eps [Default: 0.15] The smoothing factor used in the cosine smoothing function.
    * @return A vector containing the number of vaccinations for each age group at time t.
    */
    Eigen::VectorX<FP> vaccinations_at(const FP t,
                                       const CustomIndexArray<FP, AgeGroup, SimulationDay>& daily_vaccinations,
                                       const FP eps = 0.15) const
    {
        using std::floor;

        auto const& params = this->parameters;
        const FP ub        = floor(t) + 1.0;
        const FP lb        = ub - eps;

        const auto max_time = static_cast<size_t>(daily_vaccinations.template size<SimulationDay>()) - 1;

        Eigen::VectorX<FP> smoothed_vaccinations((size_t)params.get_num_groups());
        smoothed_vaccinations.setZero();

        // if daily_vaccinations is not available for the current time point, we return zero vaccinations.
        if (max_time <= floor(t)) {
            mio::log_warning("Vaccination data not available for time point ", t, ". Returning zero vaccinations.");
            return smoothed_vaccinations;
        }
        if (t >= lb) {
            for (AgeGroup age = AgeGroup(0); age < params.get_num_groups(); age++) {
                // if ub + 1 is out of bounds, we use the value at ub
                FP ubp1 = floor(ub + 1.0);
                if (max_time < ubp1) {
                    ubp1 = floor(ub);
                }
                const auto num_vaccinations_ub = daily_vaccinations[{age, SimulationDay(size_t(floor(ubp1)))}] -
                                                 daily_vaccinations[{age, SimulationDay(size_t(floor(ub)))}];
                const auto num_vaccinations_lb = daily_vaccinations[{age, SimulationDay(size_t(floor(lb + 1.0)))}] -
                                                 daily_vaccinations[{age, SimulationDay(size_t(floor(lb)))}];
                smoothed_vaccinations[(size_t)age] =
                    smoother_cosine<FP>(t, lb, ub, num_vaccinations_lb, num_vaccinations_ub);
            }
        }
        else {
            for (auto age = AgeGroup(0); age < params.get_num_groups(); age++) {
                smoothed_vaccinations[(size_t)age] = daily_vaccinations[{age, SimulationDay(size_t(floor(t + 1)))}] -
                                                     daily_vaccinations[{age, SimulationDay(size_t(floor(t)))}];
            }
        }
        return smoothed_vaccinations;
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
}; // namespace osecirts

//forward declaration, see below.
template <typename FP, class Base = mio::Simulation<FP, Model<FP>>>
class Simulation;

/**
* get percentage of infections per total population.
* @param model the compartment model with initial values.
* @param t current simulation time.
* @param y current value of compartments.
* @tparam Base simulation type that uses the SECIRS-type compartment model. see Simulation.
*/
template <typename FP, class Base = mio::Simulation<FP, Model<FP>>>
FP get_infections_relative(const Simulation<FP, Base>& model, FP t, const Eigen::Ref<const Eigen::VectorX<FP>>& y);

/**
 * specialization of compartment model simulation for the SECIRTS model.
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
    Simulation(mio::osecirts::Model<FP> const& model, FP t0 = 0., FP dt = 0.1)
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
        auto t        = BaseT::get_result().get_last_time();
        const auto dt = dyn_npis.get_interval().get();
        while (t < tmax) {

            auto dt_eff = min<FP>({dt, tmax - t, m_t_last_npi_check + dt - t});
            if (dt_eff >= 1.0) {
                dt_eff = 1.0;
            }

            BaseT::advance(t + dt_eff);
            if (t + 0.5 + dt_eff - floor(t + 0.5) >= 1) {
                this->apply_variant(t, base_infectiousness);
            }

            if (t > 0) {
                delay_npi_implementation =
                    this->get_model().parameters.template get<DynamicNPIsImplementationDelay<FP>>();
            }
            else {
                // DynamicNPIs for t=0 are 'misused' to be from-start NPIs. I.e., do not enforce delay.
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
                            this->get_model().parameters.get_start_commuter_detection() = (FP)t_start;
                            this->get_model().parameters.get_end_commuter_detection()   = (FP)t_end;
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
 * @brief Specialization of simulate for SECIRS-type models using Simulation.
 *
 * @param[in] t0 start time.
 * @param[in] tmax end time.
 * @param[in] dt time step.
 * @param[in] model SECIRS-type model to simulate.
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
 * @brief Specialization of simulate for SECIRS-type models using the FlowSimulation.
 *
 * @param[in] t0 start time.
 * @param[in] tmax end time.
 * @param[in] dt time step.
 * @param[in] model SECIRS-type model to simulate.
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
    FP inf_rel = sum_inf / sim.get_model().populations.get_total();

    return inf_rel;
}

/**
 * Get migration factors.
 * Used by migration graph simulation.
 * Like infection risk, migration of infected individuals is reduced if they are well isolated.
 * @param model the compartment model with initial values.
 * @param t current simulation time.
 * @param y current value of compartments.
 * @return vector expression, same size as y, with migration factors per compartment.
 * @tparam Base simulation type that uses a SECIRS-type compartment model. see Simulation.
 */
template <typename FP, class Base = mio::Simulation<Model<FP>, FP>>
auto get_migration_factors(const Simulation<Base>& sim, FP /*t*/, const Eigen::Ref<const Eigen::VectorX<FP>>& y)
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
                            params.template get<TestAndTraceCapacity<FP>>() * 5, p_inf.matrix(), p_inf_max.matrix());

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

/**
 * @brief Adjusts the state of commuters in a model, accounting for detection and mobility effects.
 *
 * @tparam FP Floating-point type, e.g., double.
 * @tparam Base Simulation type that uses the SECIRTS-type model. see Simulation.
 * @param[in,out] sim Simulation object containing the model and result data.
 * @param[in,out] migrated Vector representing the number of commuters in each state.
 * @param[in] time Current simulation time, used to determine the commuter detection period.
 */
template <typename FP, class Base = mio::Simulation<Model<FP>, FP>>
auto test_commuters(Simulation<FP, Base>& sim, Eigen::Ref<Eigen::VectorX<FP>> migrated, FP time)
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
        sim.get_result().get_last_value()[ISyNi] -= migrated[ISyNi] * (1 - nondetection);
        sim.get_result().get_last_value()[ISyNCi] += migrated[ISyNi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSNi] -= migrated[INSNi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSNCi] += migrated[INSNi] * (1 - nondetection);

        sim.get_result().get_last_value()[ISPIi] -= migrated[ISPIi] * (1 - nondetection);
        sim.get_result().get_last_value()[ISPICi] += migrated[ISPIi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSPIi] -= migrated[INSPIi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSPICi] += migrated[INSPIi] * (1 - nondetection);

        sim.get_result().get_last_value()[ISyIIi] -= migrated[ISyIIi] * (1 - nondetection);
        sim.get_result().get_last_value()[ISyIICi] += migrated[ISyIIi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSIIi] -= migrated[INSIIi] * (1 - nondetection);
        sim.get_result().get_last_value()[INSIICi] += migrated[INSIIi] * (1 - nondetection);

        //reduce the number of commuters
        migrated[ISyNi] *= nondetection;
        migrated[INSNi] *= nondetection;

        migrated[ISPIi] *= nondetection;
        migrated[INSPIi] *= nondetection;

        migrated[ISyIIi] *= nondetection;
        migrated[INSIIi] *= nondetection;
    }
}

} // namespace osecirts
} // namespace mio

#endif //MIO_ODE_SECIRTS_MODEL_H
