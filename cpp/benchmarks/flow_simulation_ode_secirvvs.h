/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Rene Schmieding, Daniel Abele, Martin J. Kuehn
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
#include "memilio/compartments/simulation.h"
#include "models/ode_secirvvs/model.h"

namespace mio
{
namespace benchmark
{

using FlowModel = osecirvvs::Model<>;

// For comparison benchmarks, an old model version that does not provide computation of flows has been reimplemented here.
// For more details see the original implementation in:
// https://github.com/SciCompMod/memilio/blob/13555a6b23177d2d4633c393903461a27ce5762b/cpp/models/ode_secirvvs/model.h
// Updates from Issue/PR 888:
// - Apply_variant function has been adjusted to be more generic
// - Fixed a bug where the transmission probability was always set to zero.

class FlowlessModel : public CompartmentalModel<ScalarType, osecirvvs::InfectionState,
                                                Populations<ScalarType, AgeGroup, osecirvvs::InfectionState>,
                                                osecirvvs::Parameters<ScalarType>>
{
    using InfectionState = osecirvvs::InfectionState;
    using Base           = CompartmentalModel<ScalarType, osecirvvs::InfectionState,
                                    mio::Populations<ScalarType, AgeGroup, osecirvvs::InfectionState>,
                                    osecirvvs::Parameters<ScalarType>>;

public:
    FlowlessModel(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }

    FlowlessModel(int num_agegroups)
        : FlowlessModel(Populations({AgeGroup(num_agegroups), InfectionState::Count}),
                        ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        auto const& params   = this->parameters;
        AgeGroup n_agegroups = params.get_num_groups();

        ContactMatrixGroup const& contact_matrix = params.get<osecirvvs::ContactPatterns<ScalarType>>();

        auto icu_occupancy           = 0.0;
        auto test_and_trace_required = 0.0;
        for (auto i = AgeGroup(0); i < n_agegroups; ++i) {
            test_and_trace_required +=
                (1 - params.get<osecirvvs::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i]) /
                params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] *
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
            size_t DNi  = this->populations.get_flat_index({i, InfectionState::DeadNaive});
            size_t DPIi = this->populations.get_flat_index({i, InfectionState::DeadPartialImmunity});
            size_t DIIi = this->populations.get_flat_index({i, InfectionState::DeadImprovedImmunity});

            dydt[SNi] = 0;
            dydt[ENi] = 0;

            dydt[SPIi] = 0;
            dydt[EPIi] = 0;

            dydt[SIIi] = 0;
            dydt[EIIi] = 0;

            double reducExposedPartialImmunity  = params.get<osecirvvs::ReducExposedPartialImmunity<double>>()[i];
            double reducExposedImprovedImmunity = params.get<osecirvvs::ReducExposedImprovedImmunity<double>>()[i];
            double reducInfectedSymptomsPartialImmunity =
                params.get<osecirvvs::ReducInfectedSymptomsPartialImmunity<ScalarType>>()[i];
            double reducInfectedSymptomsImprovedImmunity =
                params.get<osecirvvs::ReducInfectedSymptomsImprovedImmunity<ScalarType>>()[i];
            double reducInfectedSevereCriticalDeadPartialImmunity =
                params.get<osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[i];
            double reducInfectedSevereCriticalDeadImprovedImmunity =
                params.get<osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<ScalarType>>()[i];
            double reducTimeInfectedMild = params.get<osecirvvs::ReducTimeInfectedMild<ScalarType>>()[i];

            //symptomatic are less well quarantined when testing and tracing is overwhelmed so they infect more people
            auto riskFromInfectedSymptomatic =
                smoother_cosine(test_and_trace_required, params.get<osecirvvs::TestAndTraceCapacity<ScalarType>>(),
                                params.get<osecirvvs::TestAndTraceCapacity<ScalarType>>() * 15,
                                params.get<osecirvvs::RiskOfInfectionFromSymptomatic<ScalarType>>()[i],
                                params.get<osecirvvs::MaxRiskOfInfectionFromSymptomatic<ScalarType>>()[i]);

            auto riskFromInfectedNoSymptoms =
                smoother_cosine(test_and_trace_required, params.get<osecirvvs::TestAndTraceCapacity<ScalarType>>(),
                                params.get<osecirvvs::TestAndTraceCapacity<ScalarType>>() * 2,
                                params.get<osecirvvs::RelativeTransmissionNoSymptoms<ScalarType>>()[i], 1.0);

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
                double season_val =
                    (1 + params.get<osecirvvs::Seasonality<ScalarType>>() *
                             sin(3.141592653589793 * ((params.get<osecirvvs::StartDay>() + t) / 182.5 + 0.5)));
                double cont_freq_eff =
                    season_val * contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                 static_cast<Eigen::Index>((size_t)j));
                // without died people
                double Nj = pop[SNj] + pop[ENj] + pop[INSNj] + pop[ISyNj] + pop[ISevNj] + pop[ICrNj] + pop[INSNCj] +
                            pop[ISyNCj] + pop[SPIj] + pop[EPIj] + pop[INSPIj] + pop[ISyPIj] + pop[ISevPIj] +
                            pop[ICrPIj] + pop[INSPICj] + pop[ISyPICj] + pop[SIIj] + pop[EIIj] + pop[INSIIj] +
                            pop[ISyIIj] + pop[ISevIIj] + pop[ICrIIj] + pop[INSIICj] + pop[ISyIICj];

                double divNj = 1.0 / Nj; // precompute 1.0/Nj

                double ext_inf_force_dummy =
                    cont_freq_eff * divNj *
                    params.template get<osecirvvs::TransmissionProbabilityOnContact<ScalarType>>()[(AgeGroup)i] *
                    (riskFromInfectedNoSymptoms * (pop[INSNj] + pop[INSPIj] + pop[INSIIj]) +
                     riskFromInfectedSymptomatic * (pop[ISyNj] + pop[ISyPIj] + pop[ISyIIj]));

                double dummy_SN = y[SNi] * ext_inf_force_dummy;

                double dummy_SPI = y[SPIi] * reducExposedPartialImmunity * ext_inf_force_dummy;

                double dummy_SII = y[SIIi] * reducExposedImprovedImmunity * ext_inf_force_dummy;

                dydt[SNi] -= dummy_SN;
                dydt[ENi] += dummy_SN;

                dydt[SPIi] -= dummy_SPI;
                dydt[EPIi] += dummy_SPI;

                dydt[SIIi] -= dummy_SII;
                dydt[EIIi] += dummy_SII;
            }

            // ICU capacity shortage is close
            // TODO: if this is used with vaccination model, it has to be adapted if CriticalPerSevere
            // is different for different vaccination status. This is not the case here and in addition, ICUCapacity
            // is set to infinity and this functionality is deactivated, so this is OK for the moment.
            double criticalPerSevereAdjusted =
                smoother_cosine(icu_occupancy, 0.90 * params.get<osecirvvs::ICUCapacity<ScalarType>>(),
                                params.get<osecirvvs::ICUCapacity<ScalarType>>(),
                                params.get<osecirvvs::CriticalPerSevere<ScalarType>>()[i], 0);

            double deathsPerSevereAdjusted =
                params.get<osecirvvs::CriticalPerSevere<ScalarType>>()[i] - criticalPerSevereAdjusted;

            /**** path of immune-naive ***/

            dydt[ENi] -=
                y[ENi] /
                params.get<
                    osecirvvs::TimeExposed<ScalarType>>()[i]; // only exchange of E and InfectedNoSymptoms done here
            dydt[INSNi] = y[ENi] / params.get<osecirvvs::TimeExposed<ScalarType>>()[i] -
                          y[INSNi] / params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i];
            dydt[INSNCi] = -y[INSNCi] / params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i];

            dydt[ISyNi] = (1 - params.get<osecirvvs::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i]) /
                              params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * y[INSNi] -
                          (y[ISyNi] / params.get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[i]);
            dydt[ISyNCi] = (1 - params.get<osecirvvs::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i]) /
                               params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * y[INSNCi] -
                           (y[ISyNCi] / params.get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[i]);

            dydt[ISevNi] = params.get<osecirvvs::SeverePerInfectedSymptoms<ScalarType>>()[i] /
                               params.get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[i] * (y[ISyNi] + y[ISyNCi]) -
                           (1 / params.get<osecirvvs::TimeInfectedSevere<ScalarType>>()[i]) * y[ISevNi];
            dydt[ICrNi] = -y[ICrNi] / params.get<osecirvvs::TimeInfectedCritical<ScalarType>>()[i];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[ICrNi] +=
                criticalPerSevereAdjusted / params.get<osecirvvs::TimeInfectedSevere<ScalarType>>()[i] * y[ISevNi];

            /**** path of partially immune (e.g., one dose of vaccination) ***/

            dydt[EPIi] -=
                y[EPIi] /
                params.get<
                    osecirvvs::TimeExposed<ScalarType>>()[i]; // only exchange of E and InfectedNoSymptoms done here
            dydt[INSPIi] =
                y[EPIi] / params.get<osecirvvs::TimeExposed<ScalarType>>()[i] -
                y[INSPIi] / (params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * reducTimeInfectedMild);
            dydt[INSPICi] =
                -y[INSPICi] / (params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * reducTimeInfectedMild);
            dydt[ISyPIi] =
                (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                    (1 - params.get<osecirvvs::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i]) * y[INSPIi] /
                    (params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * reducTimeInfectedMild) -
                (y[ISyPIi] / (params.get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[i] * reducTimeInfectedMild));
            dydt[ISyPICi] =
                (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                    (1 - params.get<osecirvvs::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i]) /
                    (params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * reducTimeInfectedMild) *
                    y[INSPICi] -
                (y[ISyPICi] / (params.get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[i] * reducTimeInfectedMild));
            dydt[ISevPIi] = reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity *
                                params.get<osecirvvs::SeverePerInfectedSymptoms<ScalarType>>()[i] /
                                (params.get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[i] * reducTimeInfectedMild) *
                                (y[ISyPIi] + y[ISyPICi]) -
                            (1 / params.get<osecirvvs::TimeInfectedSevere<ScalarType>>()[i]) * y[ISevPIi];
            dydt[ICrPIi] = -(1 / params.get<osecirvvs::TimeInfectedCritical<ScalarType>>()[i]) * y[ICrPIi];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[ICrPIi] += reducInfectedSevereCriticalDeadPartialImmunity /
                            reducInfectedSevereCriticalDeadPartialImmunity * criticalPerSevereAdjusted /
                            params.get<osecirvvs::TimeInfectedSevere<ScalarType>>()[i] * y[ISevPIi];

            /**** path of twice vaccinated, here called immune although reinfection is possible now ***/

            dydt[EIIi] -=
                y[EIIi] / params.get<osecirvvs::TimeExposed<ScalarType>>()[i]; // only exchange of E and C done here

            dydt[INSIIi] =
                y[EIIi] / params.get<osecirvvs::TimeExposed<ScalarType>>()[i] -
                y[INSIIi] / (params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * reducTimeInfectedMild);
            dydt[INSIICi] =
                -y[INSIICi] / (params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * reducTimeInfectedMild);

            dydt[ISyIIi] =
                (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                    (1 - params.get<osecirvvs::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i]) * y[INSIIi] /
                    (params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * reducTimeInfectedMild) -
                (1 / (params.get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[i] * reducTimeInfectedMild)) *
                    y[ISyIIi];
            dydt[ISyIICi] =
                (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                    (1 - params.get<osecirvvs::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i]) * y[INSIICi] /
                    (params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * reducTimeInfectedMild) -
                (1 / (params.get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[i] * reducTimeInfectedMild)) *
                    y[ISyIICi];
            dydt[ISevIIi] = reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity *
                                params.get<osecirvvs::SeverePerInfectedSymptoms<ScalarType>>()[i] /
                                (params.get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[i] * reducTimeInfectedMild) *
                                (y[ISyIIi] + y[ISyIICi]) -
                            (1 / params.get<osecirvvs::TimeInfectedSevere<ScalarType>>()[i]) * y[ISevIIi];
            dydt[ICrIIi] = -(1 / params.get<osecirvvs::TimeInfectedCritical<ScalarType>>()[i]) * y[ICrIIi];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[ICrIIi] += reducInfectedSevereCriticalDeadImprovedImmunity /
                            reducInfectedSevereCriticalDeadImprovedImmunity * criticalPerSevereAdjusted /
                            params.get<osecirvvs::TimeInfectedSevere<ScalarType>>()[i] * y[ISevIIi];

            // recovered and deaths from all paths
            dydt[SIIi] += params.get<osecirvvs::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i] /
                              params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * (y[INSNi] + y[INSNCi]) +
                          (1 - params.get<osecirvvs::SeverePerInfectedSymptoms<ScalarType>>()[i]) /
                              params.get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[i] * (y[ISyNi] + y[ISyNCi]) +
                          (1 - params.get<osecirvvs::CriticalPerSevere<ScalarType>>()[i]) /
                              params.get<osecirvvs::TimeInfectedSevere<ScalarType>>()[i] * y[ISevNi] +
                          (1 - params.get<osecirvvs::DeathsPerCritical<ScalarType>>()[i]) /
                              params.get<osecirvvs::TimeInfectedCritical<ScalarType>>()[i] * y[ICrNi];

            dydt[SIIi] +=
                (1 - (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                         (1 - params.get<osecirvvs::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i])) /
                    (params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * reducTimeInfectedMild) *
                    (y[INSPIi] + y[INSPICi]) +
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity) *
                         params.get<osecirvvs::SeverePerInfectedSymptoms<ScalarType>>()[i]) /
                    (params.get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[i] * reducTimeInfectedMild) *
                    (y[ISyPIi] + y[ISyPICi]) +
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                         params.get<osecirvvs::CriticalPerSevere<ScalarType>>()[i]) /
                    params.get<osecirvvs::TimeInfectedSevere<ScalarType>>()[i] * y[ISevPIi] +
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                         params.get<osecirvvs::DeathsPerCritical<ScalarType>>()[i]) /
                    params.get<osecirvvs::TimeInfectedCritical<ScalarType>>()[i] * y[ICrPIi];

            dydt[SIIi] +=
                (1 - (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                         (1 - params.get<osecirvvs::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i])) *
                    (params.get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[i] * reducTimeInfectedMild) *
                    (y[INSIIi] + y[INSIICi]) +
                (1 - (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity) *
                         params.get<osecirvvs::SeverePerInfectedSymptoms<ScalarType>>()[i]) /
                    (params.get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[i] * reducTimeInfectedMild) *
                    (y[ISyIIi] + y[ISyIICi]) +
                (1 -
                 (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                     params.get<osecirvvs::CriticalPerSevere<ScalarType>>()[i]) /
                    params.get<osecirvvs::TimeInfectedSevere<ScalarType>>()[i] * y[ISevIIi] +
                (1 -
                 (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                     params.get<osecirvvs::DeathsPerCritical<ScalarType>>()[i]) /
                    params.get<osecirvvs::TimeInfectedCritical<ScalarType>>()[i] * y[ICrIIi];

            dydt[DNi] = params.get<osecirvvs::DeathsPerCritical<ScalarType>>()[i] /
                        params.get<osecirvvs::TimeInfectedCritical<ScalarType>>()[i] * y[ICrNi];
            dydt[DPIi] = reducInfectedSevereCriticalDeadPartialImmunity /
                         reducInfectedSevereCriticalDeadPartialImmunity *
                         params.get<osecirvvs::DeathsPerCritical<ScalarType>>()[i] /
                         params.get<osecirvvs::TimeInfectedCritical<ScalarType>>()[i] * y[ICrPIi];
            dydt[DIIi] = reducInfectedSevereCriticalDeadImprovedImmunity /
                         reducInfectedSevereCriticalDeadImprovedImmunity *
                         params.get<osecirvvs::DeathsPerCritical<ScalarType>>()[i] /
                         params.get<osecirvvs::TimeInfectedCritical<ScalarType>>()[i] * y[ICrIIi];
            // add potential, additional deaths due to ICU overflow
            dydt[DNi] +=
                deathsPerSevereAdjusted / params.get<osecirvvs::TimeInfectedSevere<ScalarType>>()[i] * y[ISevNi];
            dydt[DPIi] +=
                (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                deathsPerSevereAdjusted / params.get<osecirvvs::TimeInfectedSevere<ScalarType>>()[i] * y[ISevPIi];
            dydt[DIIi] +=
                (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                deathsPerSevereAdjusted / params.get<osecirvvs::TimeInfectedSevere<ScalarType>>()[i] * y[ISevIIi];
        }
    }

    /**
    * serialize this. 
    * @see mio::serialize
    */
    template <class IOContext>
    void serialize(IOContext& io) const
    {
        auto obj = io.create_object("FlowlessModel");
        obj.add_element("Parameters", parameters);
        obj.add_element("Populations", populations);
    }

    /**
    * deserialize an object of this class.
    * @see mio::deserialize
    */
    template <class IOContext>
    static IOResult<FlowlessModel> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("FlowlessModel");
        auto par = obj.expect_element("Parameters", Tag<ParameterSet>{});
        auto pop = obj.expect_element("Populations", Tag<Populations>{});
        return apply(
            io,
            [](auto&& par_, auto&& pop_) {
                return FlowlessModel{pop_, par_};
            },
            par, pop);
    }
};

template <class Base = mio::Simulation<ScalarType, FlowlessModel>>
double get_infections_relative(const Simulation<ScalarType, Base>& model, double t,
                               const Eigen::Ref<const Eigen::VectorXd>& y);

template <class Base>
class Simulation : public Base
{
public:
    /**
    * construct a simulation.
    * @param model the model to simulate.
    * @param t0 start time
    * @param dt time steps
    */
    Simulation(FlowlessModel const& model, double t0 = 0., double dt = 0.1)
        : Base(model, t0, dt)
        , m_t_last_npi_check(t0)
    {
    }

    void apply_variant(const double t, const CustomIndexArray<UncertainValue<ScalarType>, AgeGroup> base_infectiousness)
    {
        auto start_day             = this->get_model().parameters.template get<osecirvvs::StartDay>();
        auto start_day_new_variant = this->get_model().parameters.template get<osecirvvs::StartDayNewVariant>();

        if (start_day + t >= start_day_new_variant - 1e-10) {
            const double days_variant      = t - (start_day_new_variant - start_day);
            const double share_new_variant = std::min(1.0, 0.01 * pow(2, (1. / 7) * days_variant));
            const auto num_groups          = this->get_model().parameters.get_num_groups();
            for (auto i = AgeGroup(0); i < num_groups; ++i) {
                double new_transmission =
                    (1 - share_new_variant) * base_infectiousness[i] +
                    share_new_variant * base_infectiousness[i] *
                        this->get_model().parameters.template get<osecirvvs::InfectiousnessNewVariant<ScalarType>>()[i];
                this->get_model()
                    .parameters.template get<osecirvvs::TransmissionProbabilityOnContact<ScalarType>>()[i] =
                    new_transmission;
            }
        }
    }

    void apply_vaccination(double t)
    {
        auto t_idx        = SimulationDay((size_t)t);
        auto& params      = this->get_model().parameters;
        size_t num_groups = (size_t)params.get_num_groups();
        auto last_value   = this->get_result().get_last_value();

        auto count = (size_t)osecirvvs::InfectionState::Count;
        auto S     = (size_t)osecirvvs::InfectionState::SusceptibleNaive;
        auto SV    = (size_t)osecirvvs::InfectionState::SusceptiblePartialImmunity;
        auto R     = (size_t)osecirvvs::InfectionState::SusceptibleImprovedImmunity;

        for (size_t i = 0; i < num_groups; ++i) {

            double first_vacc;
            double full_vacc;
            if (t_idx == SimulationDay(0)) {
                first_vacc = params.template get<osecirvvs::DailyFirstVaccination<ScalarType>>()[{(AgeGroup)i, t_idx}];
                full_vacc  = params.template get<osecirvvs::DailyFullVaccination<ScalarType>>()[{(AgeGroup)i, t_idx}];
            }
            else {
                first_vacc = params.template get<osecirvvs::DailyFirstVaccination<ScalarType>>()[{(AgeGroup)i, t_idx}] -
                             params.template get<osecirvvs::DailyFirstVaccination<ScalarType>>()[{
                                 (AgeGroup)i, t_idx - SimulationDay(1)}];
                full_vacc = params.template get<osecirvvs::DailyFullVaccination<ScalarType>>()[{(AgeGroup)i, t_idx}] -
                            params.template get<osecirvvs::DailyFullVaccination<ScalarType>>()[{
                                (AgeGroup)i, t_idx - SimulationDay(1)}];
            }

            if (last_value(count * i + S) - first_vacc < 0) {
                auto corrected = 0.99 * last_value(count * i + S);
                log_warning("too many first vaccinated at time {}: setting first_vacc from {} to {}", t, first_vacc,
                            corrected);
                first_vacc = corrected;
            }

            last_value(count * i + S) -= first_vacc;
            last_value(count * i + SV) += first_vacc;

            if (last_value(count * i + SV) - full_vacc < 0) {
                auto corrected = 0.99 * last_value(count * i + SV);
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
    Eigen::Ref<Eigen::VectorXd> advance(double tmax)
    {
        auto& t_end_dyn_npis = this->get_model().parameters.get_end_dynamic_npis();
        auto& dyn_npis =
            this->get_model().parameters.template get<osecirvvs::DynamicNPIsInfectedSymptoms<ScalarType>>();
        auto& contact_patterns  = this->get_model().parameters.template get<osecirvvs::ContactPatterns<ScalarType>>();
        const size_t num_groups = (size_t)this->get_model().parameters.get_num_groups();

        auto base_infectiousness =
            this->get_model().parameters.template get<osecirvvs::TransmissionProbabilityOnContact<ScalarType>>();

        double delay_npi_implementation;
        auto t        = Base::get_result().get_last_time();
        const auto dt = dyn_npis.get_interval().get();
        while (t < tmax) {

            auto dt_eff = std::min({dt, tmax - t, m_t_last_npi_check + dt - t});
            if (dt_eff >= 1.0) {
                dt_eff = 1.0;
            }

            if (t == 0) {
                //this->apply_vaccination(t); // done in init now?
                this->apply_variant(t, base_infectiousness);
            }
            Base::advance(t + dt_eff);
            if (t + 0.5 + dt_eff - std::floor(t + 0.5) >= 1) {
                this->apply_vaccination(t + 0.5 + dt_eff);
                this->apply_variant(t, base_infectiousness);
            }

            if (t > 0) {
                delay_npi_implementation = 7;
            }
            else {
                delay_npi_implementation = 0;
            }
            t = t + dt_eff;

            if (dyn_npis.get_thresholds().size() > 0) {
                if (floating_point_greater_equal(t, m_t_last_npi_check + dt)) {
                    if (t < t_end_dyn_npis) {
                        auto inf_rel = get_infections_relative(*this, t, this->get_result().get_last_value()) *
                                       dyn_npis.get_base_value();
                        auto exceeded_threshold = dyn_npis.get_max_exceeded_threshold(inf_rel);
                        if (exceeded_threshold != dyn_npis.get_thresholds().end() &&
                            (exceeded_threshold->first > m_dynamic_npi.first ||
                             t > double(m_dynamic_npi.second))) { //old npi was weaker or is expired

                            auto t_start = SimulationTime(t + delay_npi_implementation);
                            auto t_end   = t_start + SimulationTime(dyn_npis.get_duration());
                            this->get_model().parameters.get_start_commuter_detection() = (double)t_start;
                            this->get_model().parameters.get_end_commuter_detection()   = (double)t_end;
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

        this->get_model().parameters.template get<osecirvvs::TransmissionProbabilityOnContact<ScalarType>>() =
            base_infectiousness;

        return this->get_result().get_last_value();
    }

private:
    double m_t_last_npi_check;
    std::pair<double, SimulationTime> m_dynamic_npi = {-std::numeric_limits<double>::max(), SimulationTime(0)};
};

template <class Base>
double get_infections_relative(const Simulation<Base>& sim, double /*t*/, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    double sum_inf = 0;
    for (auto i = AgeGroup(0); i < sim.get_model().parameters.get_num_groups(); ++i) {
        sum_inf += sim.get_model().populations.get_from(y, {i, osecirvvs::InfectionState::InfectedSymptomsNaive});
        sum_inf +=
            sim.get_model().populations.get_from(y, {i, osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed});
        sum_inf +=
            sim.get_model().populations.get_from(y, {i, osecirvvs::InfectionState::InfectedSymptomsPartialImmunity});
        sum_inf +=
            sim.get_model().populations.get_from(y, {i, osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity});
        sum_inf += sim.get_model().populations.get_from(
            y, {i, osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed});
        sum_inf += sim.get_model().populations.get_from(
            y, {i, osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed});
    }
    auto inf_rel = sum_inf / sim.get_model().populations.get_total();

    return inf_rel;
}

// setup taken from examples/secirvvs.cpp:main
template <class Model>
void setup_model(Model& model)
{
    for (AgeGroup i = 0; i < model.parameters.get_num_groups(); i++) {
        model.populations[{i, osecirvvs::InfectionState::ExposedNaive}]                                = 10;
        model.populations[{i, osecirvvs::InfectionState::ExposedImprovedImmunity}]                     = 11;
        model.populations[{i, osecirvvs::InfectionState::ExposedPartialImmunity}]                      = 12;
        model.populations[{i, osecirvvs::InfectionState::InfectedNoSymptomsNaive}]                     = 13;
        model.populations[{i, osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}]            = 13;
        model.populations[{i, osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}]           = 14;
        model.populations[{i, osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]  = 14;
        model.populations[{i, osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}]          = 15;
        model.populations[{i, osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 15;
        model.populations[{i, osecirvvs::InfectionState::InfectedSymptomsNaive}]                       = 5;
        model.populations[{i, osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]              = 5;
        model.populations[{i, osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}]             = 6;
        model.populations[{i, osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]    = 6;
        model.populations[{i, osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]            = 7;
        model.populations[{i, osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]   = 7;
        model.populations[{i, osecirvvs::InfectionState::InfectedSevereNaive}]                         = 8;
        model.populations[{i, osecirvvs::InfectionState::InfectedSevereImprovedImmunity}]              = 1;
        model.populations[{i, osecirvvs::InfectionState::InfectedSeverePartialImmunity}]               = 2;
        model.populations[{i, osecirvvs::InfectionState::InfectedCriticalNaive}]                       = 3;
        model.populations[{i, osecirvvs::InfectionState::InfectedCriticalPartialImmunity}]             = 4;
        model.populations[{i, osecirvvs::InfectionState::InfectedCriticalImprovedImmunity}]            = 5;
        model.populations[{i, osecirvvs::InfectionState::SusceptibleImprovedImmunity}]                 = 6;
        model.populations[{i, osecirvvs::InfectionState::SusceptiblePartialImmunity}]                  = 7;
        model.populations[{(AgeGroup)0, osecirvvs::InfectionState::DeadNaive}]                         = 0;
        model.populations[{(AgeGroup)0, osecirvvs::InfectionState::DeadPartialImmunity}]               = 0;
        model.populations[{(AgeGroup)0, osecirvvs::InfectionState::DeadImprovedImmunity}]              = 0;
        model.populations.template set_difference_from_group_total<AgeGroup>(
            {i, osecirvvs::InfectionState::SusceptibleNaive}, 1000);
    }

    model.parameters.template get<osecirvvs::ICUCapacity<ScalarType>>()          = 100;
    model.parameters.template get<osecirvvs::TestAndTraceCapacity<ScalarType>>() = 0.0143;
    model.parameters.template get<osecirvvs::DailyFirstVaccination<ScalarType>>().resize(SimulationDay(size_t(1000)));
    model.parameters.template get<osecirvvs::DailyFirstVaccination<ScalarType>>().array().setConstant(5);
    model.parameters.template get<osecirvvs::DailyFullVaccination<ScalarType>>().resize(SimulationDay(size_t(1000)));
    model.parameters.template get<osecirvvs::DailyFullVaccination<ScalarType>>().array().setConstant(3);

    auto& contacts       = model.parameters.template get<osecirvvs::ContactPatterns<ScalarType>>();
    auto& contact_matrix = contacts.get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(0.5);
    contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
    contact_matrix[0].add_damping(0.3, SimulationTime(5.0));

    //times
    model.parameters.template get<osecirvvs::TimeExposed<ScalarType>>()[AgeGroup(0)]            = 3.33;
    model.parameters.template get<osecirvvs::TimeInfectedNoSymptoms<ScalarType>>()[AgeGroup(0)] = 1.87;
    model.parameters.template get<osecirvvs::TimeInfectedSymptoms<ScalarType>>()[AgeGroup(0)]   = 7;
    model.parameters.template get<osecirvvs::TimeInfectedSevere<ScalarType>>()[AgeGroup(0)]     = 6;
    model.parameters.template get<osecirvvs::TimeInfectedCritical<ScalarType>>()[AgeGroup(0)]   = 7;

    //probabilities
    model.parameters.template get<osecirvvs::TransmissionProbabilityOnContact<ScalarType>>()[AgeGroup(0)] = 0.15;
    model.parameters.template get<osecirvvs::RelativeTransmissionNoSymptoms<ScalarType>>()[AgeGroup(0)]   = 0.5;
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    model.parameters.template get<osecirvvs::RiskOfInfectionFromSymptomatic<ScalarType>>()[AgeGroup(0)]    = 0.0;
    model.parameters.template get<osecirvvs::MaxRiskOfInfectionFromSymptomatic<ScalarType>>()[AgeGroup(0)] = 0.4;
    model.parameters.template get<osecirvvs::RecoveredPerInfectedNoSymptoms<ScalarType>>()[AgeGroup(0)]    = 0.2;
    model.parameters.template get<osecirvvs::SeverePerInfectedSymptoms<ScalarType>>()[AgeGroup(0)]         = 0.1;
    model.parameters.template get<osecirvvs::CriticalPerSevere<ScalarType>>()[AgeGroup(0)]                 = 0.1;
    model.parameters.template get<osecirvvs::DeathsPerCritical<ScalarType>>()[AgeGroup(0)]                 = 0.1;

    model.parameters.template get<osecirvvs::ReducExposedPartialImmunity<ScalarType>>()[AgeGroup(0)]           = 0.8;
    model.parameters.template get<osecirvvs::ReducExposedImprovedImmunity<ScalarType>>()[AgeGroup(0)]          = 0.331;
    model.parameters.template get<osecirvvs::ReducInfectedSymptomsPartialImmunity<ScalarType>>()[AgeGroup(0)]  = 0.65;
    model.parameters.template get<osecirvvs::ReducInfectedSymptomsImprovedImmunity<ScalarType>>()[AgeGroup(0)] = 0.243;
    model.parameters
        .template get<osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<ScalarType>>()[AgeGroup(0)] = 0.1;
    model.parameters
        .template get<osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<ScalarType>>()[AgeGroup(0)] = 0.091;
    model.parameters.template get<osecirvvs::ReducTimeInfectedMild<ScalarType>>()[AgeGroup(0)]               = 0.9;

    model.parameters.template get<osecirvvs::Seasonality<ScalarType>>() = 0.2;

    model.apply_constraints();
}

} // namespace benchmark
} // namespace mio
