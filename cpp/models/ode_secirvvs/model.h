/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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
#ifndef ODESECIRVVS_MODEL_H
#define ODESECIRVVS_MODEL_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/compartments/simulation.h"
#include "memilio/epidemiology/populations.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/parameters.h"
#include "memilio/math/smoother.h"
#include "memilio/math/eigen_util.h"

namespace mio
{
namespace osecirvvs
{
class Model : public CompartmentalModel<InfectionState, Populations<AgeGroup, InfectionState>, Parameters>
{
    using Base = CompartmentalModel<InfectionState, mio::Populations<AgeGroup, InfectionState>, Parameters>;

public:
    Model(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }

    Model(int num_agegroups)
        : Model(Populations({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        auto const& params   = this->parameters;
        AgeGroup n_agegroups = params.get_num_groups();

        ContactMatrixGroup const& contact_matrix = params.get<ContactPatterns>();

        auto icu_occupancy           = 0.0;
        auto test_and_trace_required = 0.0;
        for (auto i = AgeGroup(0); i < n_agegroups; ++i) {
            auto rateINS = 0.5 / (params.get<IncubationTime>()[i] - params.get<SerialInterval>()[i]);
            test_and_trace_required +=
                (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * rateINS *
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

            size_t ITi = this->populations.get_flat_index({i, InfectionState::TotalInfections});

            dydt[SNi] = 0;
            dydt[ENi] = 0;

            dydt[SPIi] = 0;
            dydt[EPIi] = 0;

            dydt[SIIi] = 0;
            dydt[EIIi] = 0;

            double rateE   = 1.0 / (2 * params.get<SerialInterval>()[i] - params.get<IncubationTime>()[i]);
            double rateINS = 0.5 / (params.get<IncubationTime>()[i] - params.get<SerialInterval>()[i]);

            double reducExposedPartialImmunity           = params.get<ReducExposedPartialImmunity>()[i];
            double reducExposedImprovedImmunity          = params.get<ReducExposedImprovedImmunity>()[i];
            double reducInfectedSymptomsPartialImmunity  = params.get<ReducInfectedSymptomsPartialImmunity>()[i];
            double reducInfectedSymptomsImprovedImmunity = params.get<ReducInfectedSymptomsImprovedImmunity>()[i];
            double reducInfectedSevereCriticalDeadPartialImmunity =
                params.get<ReducInfectedSevereCriticalDeadPartialImmunity>()[i];
            double reducInfectedSevereCriticalDeadImprovedImmunity =
                params.get<ReducInfectedSevereCriticalDeadImprovedImmunity>()[i];
            double reducTimeInfectedMild = params.get<ReducTimeInfectedMild>()[i];

            //symptomatic are less well quarantined when testing and tracing is overwhelmed so they infect more people
            auto riskFromInfectedSymptomatic = smoother_cosine(
                test_and_trace_required, params.get<TestAndTraceCapacity>(), params.get<TestAndTraceCapacity>() * 15,
                params.get<RiskOfInfectionFromSymptomatic>()[i], params.get<MaxRiskOfInfectionFromSymptomatic>()[i]);

            auto riskFromInfectedNoSymptoms = smoother_cosine(
                test_and_trace_required, params.get<TestAndTraceCapacity>(), params.get<TestAndTraceCapacity>() * 2,
                params.get<RelativeTransmissionNoSymptoms>()[i], 1.0);

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
                    (1 + params.get<Seasonality>() *
                             sin(3.141592653589793 * (std::fmod((params.get<StartDay>() + t), 365.0) / 182.5 + 0.5)));

                double cont_freq_eff =
                    season_val * contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                 static_cast<Eigen::Index>((size_t)j));
                // without died people
                double Nj = pop[SNj] + pop[ENj] + pop[INSNj] + pop[ISyNj] + pop[ISevNj] + pop[ICrNj] + pop[INSNCj] +
                            pop[ISyNCj] + pop[SPIj] + pop[EPIj] + pop[INSPIj] + pop[ISyPIj] + pop[ISevPIj] +
                            pop[ICrPIj] + pop[INSPICj] + pop[ISyPICj] + pop[SIIj] + pop[EIIj] + pop[INSIIj] +
                            pop[ISyIIj] + pop[ISevIIj] + pop[ICrIIj] + pop[INSIICj] + pop[ISyIICj];
                double divNj = 0;
                // only needed in case of mobility nodes since some age groups are fordbidden to change regions.
                if (Nj > 0) {
                    divNj = 1.0 / Nj; // precompute 1.0/Nj
                }

                double ext_inf_force_dummy = cont_freq_eff * divNj *
                                             params.template get<TransmissionProbabilityOnContact>()[(AgeGroup)i] *
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
                smoother_cosine(icu_occupancy, 0.90 * params.get<ICUCapacity>(), params.get<ICUCapacity>(),
                                params.get<CriticalPerSevere>()[i], 0);

            double deathsPerSevereAdjusted = params.get<CriticalPerSevere>()[i] - criticalPerSevereAdjusted;

            /**** path of immune-naive ***/

            dydt[ENi] -= rateE * y[ENi]; // only exchange of E and InfectedNoSymptoms done here
            dydt[INSNi]  = rateE * y[ENi] - rateINS * y[INSNi];
            dydt[INSNCi] = -rateINS * y[INSNCi];

            dydt[ISyNi] = (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * rateINS * y[INSNi] -
                          (y[ISyNi] / params.get<TimeInfectedSymptoms>()[i]);
            dydt[ISyNCi] = (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * rateINS * y[INSNCi] -
                           (y[ISyNCi] / params.get<TimeInfectedSymptoms>()[i]);

            dydt[ISevNi] = params.get<SeverePerInfectedSymptoms>()[i] / params.get<TimeInfectedSymptoms>()[i] *
                               (y[ISyNi] + y[ISyNCi]) -
                           (1 / params.get<TimeInfectedSevere>()[i]) * y[ISevNi];
            dydt[ICrNi] = -y[ICrNi] / params.get<TimeInfectedCritical>()[i];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[ICrNi] += criticalPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevNi];

            /**** path of partially immune (e.g., one dose of vaccination) ***/

            dydt[EPIi] -= rateE * y[EPIi]; // only exchange of E and InfectedNoSymptoms done here
            dydt[INSPIi]  = rateE * y[EPIi] - (rateINS / reducTimeInfectedMild) * y[INSPIi];
            dydt[INSPICi] = -(rateINS / reducTimeInfectedMild) * y[INSPICi];
            dydt[ISyPIi]  = (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                               (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) *
                               (rateINS / reducTimeInfectedMild) * y[INSPIi] -
                           (y[ISyPIi] / (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild));
            dydt[ISyPICi] = (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                                (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) *
                                (rateINS / reducTimeInfectedMild) * y[INSPICi] -
                            (y[ISyPICi] / (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild));
            dydt[ISevPIi] = reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity *
                                params.get<SeverePerInfectedSymptoms>()[i] /
                                (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild) *
                                (y[ISyPIi] + y[ISyPICi]) -
                            (1 / params.get<TimeInfectedSevere>()[i]) * y[ISevPIi];
            dydt[ICrPIi] = -(1 / params.get<TimeInfectedCritical>()[i]) * y[ICrPIi];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[ICrPIi] += reducInfectedSevereCriticalDeadPartialImmunity /
                            reducInfectedSevereCriticalDeadPartialImmunity * criticalPerSevereAdjusted /
                            params.get<TimeInfectedSevere>()[i] * y[ISevPIi];

            /**** path of twice vaccinated, here called immune although reinfection is possible now ***/

            dydt[EIIi] -= rateE * y[EIIi]; // only exchange of E and C done here

            dydt[INSIIi]  = rateE * y[EIIi] - (rateINS / reducTimeInfectedMild) * y[INSIIi];
            dydt[INSIICi] = -(rateINS / reducTimeInfectedMild) * y[INSIICi];

            dydt[ISyIIi] = (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                               (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) *
                               (rateINS / reducTimeInfectedMild) * y[INSIIi] -
                           (1 / (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild)) * y[ISyIIi];
            dydt[ISyIICi] = (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                                (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) *
                                (rateINS / reducTimeInfectedMild) * y[INSIICi] -
                            (1 / (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild)) * y[ISyIICi];
            dydt[ISevIIi] = reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity *
                                params.get<SeverePerInfectedSymptoms>()[i] /
                                (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild) *
                                (y[ISyIIi] + y[ISyIICi]) -
                            (1 / params.get<TimeInfectedSevere>()[i]) * y[ISevIIi];
            dydt[ICrIIi] = -(1 / params.get<TimeInfectedCritical>()[i]) * y[ICrIIi];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[ICrIIi] += reducInfectedSevereCriticalDeadImprovedImmunity /
                            reducInfectedSevereCriticalDeadImprovedImmunity * criticalPerSevereAdjusted /
                            params.get<TimeInfectedSevere>()[i] * y[ISevIIi];

            // compute auxiliary compartment of all past infections
            dydt[ITi] =
                (1 / params.get<TimeInfectedSymptoms>()[i]) * (y[ISyNi] + y[ISyNCi]) +
                (1 / (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild)) * (y[ISyPIi] + y[ISyPICi]) +
                (1 / (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild)) * (y[ISyIIi] + y[ISyIICi]);

            // recovered and deaths from all paths
            dydt[SIIi] += params.get<RecoveredPerInfectedNoSymptoms>()[i] * rateINS * (y[INSNi] + y[INSNCi]) +
                          (1 - params.get<SeverePerInfectedSymptoms>()[i]) / params.get<TimeInfectedSymptoms>()[i] *
                              (y[ISyNi] + y[ISyNCi]) +
                          (1 - params.get<CriticalPerSevere>()[i]) / params.get<TimeInfectedSevere>()[i] * y[ISevNi] +
                          (1 - params.get<DeathsPerCritical>()[i]) / params.get<TimeInfectedCritical>()[i] * y[ICrNi];

            dydt[SIIi] +=
                (1 - (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                         (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i])) *
                    rateINS / reducTimeInfectedMild * (y[INSPIi] + y[INSPICi]) +
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity) *
                         params.get<SeverePerInfectedSymptoms>()[i]) /
                    (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild) * (y[ISyPIi] + y[ISyPICi]) +
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                         params.get<CriticalPerSevere>()[i]) /
                    params.get<TimeInfectedSevere>()[i] * y[ISevPIi] +
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                         params.get<DeathsPerCritical>()[i]) /
                    params.get<TimeInfectedCritical>()[i] * y[ICrPIi];

            dydt[SIIi] +=
                (1 - (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                         (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i])) *
                    rateINS / reducTimeInfectedMild * (y[INSIIi] + y[INSIICi]) +
                (1 - (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity) *
                         params.get<SeverePerInfectedSymptoms>()[i]) /
                    (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild) * (y[ISyIIi] + y[ISyIICi]) +
                (1 -
                 (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                     params.get<CriticalPerSevere>()[i]) /
                    params.get<TimeInfectedSevere>()[i] * y[ISevIIi] +
                (1 -
                 (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                     params.get<DeathsPerCritical>()[i]) /
                    params.get<TimeInfectedCritical>()[i] * y[ICrIIi];

            dydt[DNi]  = params.get<DeathsPerCritical>()[i] / params.get<TimeInfectedCritical>()[i] * y[ICrNi];
            dydt[DPIi] = reducInfectedSevereCriticalDeadPartialImmunity /
                         reducInfectedSevereCriticalDeadPartialImmunity * params.get<DeathsPerCritical>()[i] /
                         params.get<TimeInfectedCritical>()[i] * y[ICrPIi];
            dydt[DIIi] = reducInfectedSevereCriticalDeadImprovedImmunity /
                         reducInfectedSevereCriticalDeadImprovedImmunity * params.get<DeathsPerCritical>()[i] /
                         params.get<TimeInfectedCritical>()[i] * y[ICrIIi];
            // add potential, additional deaths due to ICU overflow
            dydt[DNi] += deathsPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevNi];
            dydt[DPIi] +=
                (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                deathsPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevPIi];
            dydt[DIIi] +=
                (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                deathsPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevIIi];
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
        obj.add_element("Parameters", parameters);
        obj.add_element("Populations", populations);
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
template <class Base = mio::Simulation<Model>>
class Simulation;

/**
* get percentage of infections per total population.
* @param model the compartment model with initial values.
* @param t current simulation time.
* @param y current value of compartments.
* @tparam Base simulation type that uses a secir compartment model. see Simulation.
*/
template <class Base = mio::Simulation<Model>>
double get_infections_relative(const Simulation<Base>& model, double t, const Eigen::Ref<const Eigen::VectorXd>& y);

/**
    * specialization of compartment model simulation for the SECIRVVS model.
    * @tparam Base simulation type, default mio::Simulation. For testing purposes only!
    */
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
    Simulation(Model const& model, double t0 = 0., double dt = 0.1)
        : Base(model, t0, dt)
        , m_t_last_npi_check(t0)
    {
    }

    void apply_b161(double t)
    {

        auto start_day   = this->get_model().parameters.template get<StartDay>();
        auto b161_growth = (start_day - get_day_in_year(Date(2021, 6, 6))) * 0.1666667;
        // 2 equal to the share of the delta variant on June 6
        double share_new_variant = std::min(1.0, pow(2, t * 0.1666667 + b161_growth) * 0.01);
        size_t num_groups        = (size_t)this->get_model().parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; ++i) {
            double new_transmission =
                (1 - share_new_variant) *
                    this->get_model().parameters.template get<BaseInfectiousnessB117>()[(AgeGroup)i] +
                share_new_variant * this->get_model().parameters.template get<BaseInfectiousnessB161>()[(AgeGroup)i];
            this->get_model().parameters.template get<TransmissionProbabilityOnContact>()[(AgeGroup)i] =
                new_transmission;
        }
    }

    void apply_vaccination(double t)
    {
        auto t_idx        = SimulationDay((size_t)t);
        auto& params      = this->get_model().parameters;
        size_t num_groups = (size_t)params.get_num_groups();
        auto last_value   = this->get_result().get_last_value();

        auto count = (size_t)InfectionState::Count;
        auto S     = (size_t)InfectionState::SusceptibleNaive;
        auto SV    = (size_t)InfectionState::SusceptiblePartialImmunity;
        auto R     = (size_t)InfectionState::SusceptibleImprovedImmunity;

        for (size_t i = 0; i < num_groups; ++i) {

            double first_vacc;
            double full_vacc;
            if (t_idx == SimulationDay(0)) {
                first_vacc = params.template get<DailyFirstVaccination>()[{(AgeGroup)i, t_idx}];
                full_vacc  = params.template get<DailyFullVaccination>()[{(AgeGroup)i, t_idx}];
            }
            else {
                first_vacc = params.template get<DailyFirstVaccination>()[{(AgeGroup)i, t_idx}] -
                             params.template get<DailyFirstVaccination>()[{(AgeGroup)i, t_idx - SimulationDay(1)}];
                full_vacc = params.template get<DailyFullVaccination>()[{(AgeGroup)i, t_idx}] -
                            params.template get<DailyFullVaccination>()[{(AgeGroup)i, t_idx - SimulationDay(1)}];
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
        auto& t_end_dyn_npis   = this->get_model().parameters.get_end_dynamic_npis();
        auto& dyn_npis         = this->get_model().parameters.template get<DynamicNPIsInfectedSymptoms>();
        auto& contact_patterns = this->get_model().parameters.template get<ContactPatterns>();

        double delay_lockdown;
        auto t        = Base::get_result().get_last_time();
        const auto dt = dyn_npis.get_interval().get();
        while (t < tmax) {

            auto dt_eff = std::min({dt, tmax - t, m_t_last_npi_check + dt - t});
            if (dt_eff >= 1.0) {
                dt_eff = 1.0;
            }

            if (t == 0) {
                //this->apply_vaccination(t); // done in init now?
                this->apply_b161(t);
            }
            Base::advance(t + dt_eff);
            if (t + 0.5 + dt_eff - std::floor(t + 0.5) >= 1) {
                // this->apply_vaccination(t + 0.5 + dt_eff);
                this->apply_b161(t);
            }

            if (t > 0) {
                delay_lockdown = 7;
            }
            else {
                delay_lockdown = 0;
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

                            auto t_start = SimulationTime(t + delay_lockdown);
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
        return this->get_result().get_last_value();
    }

private:
    double m_t_last_npi_check;
    std::pair<double, SimulationTime> m_dynamic_npi = {-std::numeric_limits<double>::max(), SimulationTime(0)};
};

/**
* Run simulation using a SECIRVVS model.
* @param t0 start time.
* @param tmax end time.
* @param dt time step.
* @param model secir model to simulate.
* @param integrator optional integrator, uses rk45 if nullptr.
*/
inline auto simulate(double t0, double tmax, double dt, const Model& model,
                     std::shared_ptr<IntegratorCore> integrator = nullptr)
{
    return mio::simulate<Model, Simulation<>>(t0, tmax, dt, model, integrator);
}

//see declaration above.
template <class Base>
double get_infections_relative(const Simulation<Base>& sim, double /*t*/, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    double sum_inf = 0;
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
* Get migration factors.
* Used by migration graph simulation.
* Like infection risk, migration of infected individuals is reduced if they are well isolated.
* @param model the compartment model with initial values.
* @param t current simulation time.
* @param y current value of compartments.
* @return vector expression, same size as y, with migration factors per compartment.
* @tparam Base simulation type that uses a secir compartment model. see Simulation.
*/
template <class Base = mio::Simulation<Model>>
auto get_migration_factors(const Simulation<Base>& sim, double /*t*/, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    auto& params = sim.get_model().parameters;
    //parameters as arrays
    auto& t_inc     = params.template get<IncubationTime>().array().template cast<double>();
    auto& t_ser     = params.template get<SerialInterval>().array().template cast<double>();
    auto& p_asymp   = params.template get<RecoveredPerInfectedNoSymptoms>().array().template cast<double>();
    auto& p_inf     = params.template get<RiskOfInfectionFromSymptomatic>().array().template cast<double>();
    auto& p_inf_max = params.template get<MaxRiskOfInfectionFromSymptomatic>().array().template cast<double>();
    //slice of InfectedNoSymptoms
    auto y_car = slice(y, {Eigen::Index(InfectionState::InfectedNoSymptomsNaive),
                           Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)}) +
                 slice(y, {Eigen::Index(InfectionState::InfectedNoSymptomsPartialImmunity),
                           Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)}) +
                 slice(y, {Eigen::Index(InfectionState::InfectedNoSymptomsImprovedImmunity),
                           Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)});

    //compute isolation, same as infection risk from main model
    auto R3                      = 0.5 / (t_inc - t_ser);
    auto test_and_trace_required = ((1 - p_asymp) * R3 * y_car.array()).sum();
    auto riskFromInfectedSymptomatic =
        smoother_cosine(test_and_trace_required, double(params.template get<TestAndTraceCapacity>()),
                        params.template get<TestAndTraceCapacity>() * 5, p_inf.matrix(), p_inf_max.matrix());

    //set factor for infected
    auto factors = Eigen::VectorXd::Ones(y.rows()).eval();
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

template <class Base = mio::Simulation<Model>>
auto test_commuters(Simulation<Base>& sim, Eigen::Ref<Eigen::VectorXd> migrated, double time)
{
    auto& model       = sim.get_model();
    auto nondetection = 1.0;
    if (time >= model.parameters.get_start_commuter_detection() &&
        time < model.parameters.get_end_commuter_detection()) {
        nondetection = (double)model.parameters.get_commuter_nondetection();
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

} // namespace osecirvvs
} // namespace mio

#endif //ODESECIRVVS_MODEL_H
