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

#include "memilio/compartments/flow_model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/populations.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/parameters.h"
#include "memilio/math/smoother.h"
#include "memilio/math/eigen_util.h"
#include <cstddef>
#include <vector>

namespace mio
{
namespace osecirvvs
{
// clang-format off
using Flows = TypeList<
    //naive
    Flow<InfectionState::SusceptibleNaive,                            InfectionState::ExposedNaive>, 
    Flow<InfectionState::SusceptibleNaive,                            InfectionState::TemporaryImmunity1>, 
    Flow<InfectionState::ExposedNaive,                                InfectionState::InfectedNoSymptomsNaive>,
    Flow<InfectionState::InfectedNoSymptomsNaive,                     InfectionState::InfectedSymptomsNaive>,
    Flow<InfectionState::InfectedNoSymptomsNaive,                     InfectionState::TemporaryImmunity1>,
    Flow<InfectionState::InfectedNoSymptomsNaiveConfirmed,            InfectionState::InfectedSymptomsNaiveConfirmed>,
    Flow<InfectionState::InfectedNoSymptomsNaiveConfirmed,            InfectionState::TemporaryImmunity1>,
    Flow<InfectionState::InfectedSymptomsNaive,                       InfectionState::InfectedSevereNaive>,
    Flow<InfectionState::InfectedSymptomsNaive,                       InfectionState::TemporaryImmunity1>,
    Flow<InfectionState::InfectedSymptomsNaiveConfirmed,              InfectionState::InfectedSevereNaive>,
    Flow<InfectionState::InfectedSymptomsNaiveConfirmed,              InfectionState::TemporaryImmunity1>,
    Flow<InfectionState::InfectedSevereNaive,                         InfectionState::InfectedCriticalNaive>,
    Flow<InfectionState::InfectedSevereNaive,                         InfectionState::TemporaryImmunity1>, 
    Flow<InfectionState::InfectedSevereNaive,                         InfectionState::DeadNaive>,
    Flow<InfectionState::InfectedCriticalNaive,                       InfectionState::DeadNaive>,
    Flow<InfectionState::InfectedCriticalNaive,                       InfectionState::TemporaryImmunity1>,
    Flow<InfectionState::TemporaryImmunity1,                          InfectionState::SusceptiblePartialImmunity>,
    //partial immunity
    Flow<InfectionState::SusceptiblePartialImmunity,                  InfectionState::ExposedPartialImmunity>,
    Flow<InfectionState::SusceptiblePartialImmunity,                  InfectionState::TemporaryImmunity2>,
    Flow<InfectionState::ExposedPartialImmunity,                      InfectionState::InfectedNoSymptomsPartialImmunity>,
    Flow<InfectionState::InfectedNoSymptomsPartialImmunity,           InfectionState::InfectedSymptomsPartialImmunity>,
    Flow<InfectionState::InfectedNoSymptomsPartialImmunity,           InfectionState::TemporaryImmunity2>,
    Flow<InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,  InfectionState::InfectedSymptomsPartialImmunityConfirmed>,
    Flow<InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,  InfectionState::TemporaryImmunity2>,
    Flow<InfectionState::InfectedSymptomsPartialImmunity,             InfectionState::InfectedSeverePartialImmunity>,
    Flow<InfectionState::InfectedSymptomsPartialImmunity,             InfectionState::TemporaryImmunity2>,
    Flow<InfectionState::InfectedSymptomsPartialImmunityConfirmed,    InfectionState::InfectedSeverePartialImmunity>,
    Flow<InfectionState::InfectedSymptomsPartialImmunityConfirmed,    InfectionState::TemporaryImmunity2>,
    Flow<InfectionState::InfectedSeverePartialImmunity,               InfectionState::InfectedCriticalPartialImmunity>,
    Flow<InfectionState::InfectedSeverePartialImmunity,               InfectionState::TemporaryImmunity2>,
    Flow<InfectionState::InfectedSeverePartialImmunity,               InfectionState::DeadPartialImmunity>,
    Flow<InfectionState::InfectedCriticalPartialImmunity,             InfectionState::DeadPartialImmunity>,
    Flow<InfectionState::InfectedCriticalPartialImmunity,             InfectionState::TemporaryImmunity2>,
    //improved immunity
    Flow<InfectionState::SusceptibleImprovedImmunity,                 InfectionState::ExposedImprovedImmunity>,
    Flow<InfectionState::SusceptibleImprovedImmunity,                 InfectionState::TemporaryImmunity2>,
    Flow<InfectionState::ExposedImprovedImmunity,                     InfectionState::InfectedNoSymptomsImprovedImmunity>,
    Flow<InfectionState::InfectedNoSymptomsImprovedImmunity,          InfectionState::InfectedSymptomsImprovedImmunity>,
    Flow<InfectionState::InfectedNoSymptomsImprovedImmunity,          InfectionState::TemporaryImmunity2>,
    Flow<InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed, InfectionState::InfectedSymptomsImprovedImmunityConfirmed>,
    Flow<InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed, InfectionState::TemporaryImmunity2>,
    Flow<InfectionState::InfectedSymptomsImprovedImmunity,            InfectionState::InfectedSevereImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsImprovedImmunity,            InfectionState::TemporaryImmunity2>,
    Flow<InfectionState::InfectedSymptomsImprovedImmunityConfirmed,   InfectionState::InfectedSevereImprovedImmunity>,
    Flow<InfectionState::InfectedSymptomsImprovedImmunityConfirmed,   InfectionState::TemporaryImmunity2>,
    Flow<InfectionState::InfectedSevereImprovedImmunity,              InfectionState::InfectedCriticalImprovedImmunity>,
    Flow<InfectionState::InfectedSevereImprovedImmunity,              InfectionState::TemporaryImmunity2>,
    Flow<InfectionState::InfectedSevereImprovedImmunity,              InfectionState::DeadImprovedImmunity>,
    Flow<InfectionState::InfectedCriticalImprovedImmunity,            InfectionState::DeadImprovedImmunity>,
    Flow<InfectionState::InfectedCriticalImprovedImmunity,            InfectionState::TemporaryImmunity2>,
    
    // waning
    Flow<InfectionState::TemporaryImmunity1,                          InfectionState::SusceptiblePartialImmunity>,
    Flow<InfectionState::TemporaryImmunity2,                          InfectionState::SusceptibleImprovedImmunity>,
    Flow<InfectionState::SusceptibleImprovedImmunity,                 InfectionState::SusceptiblePartialImmunity>,
    Flow<InfectionState::SusceptiblePartialImmunity,                  InfectionState::SusceptibleNaive>>;

// clang-format on

class Model : public FlowModel<InfectionState, Populations<AgeGroup, InfectionState>, Parameters, Flows>
{
    using Base = FlowModel<InfectionState, mio::Populations<AgeGroup, InfectionState>, Parameters, Flows>;

public:
    Model(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }

    Model(int num_agegroups)
        : Model(Populations({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                   Eigen::Ref<Eigen::VectorXd> flows) const override
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
            size_t TImm1 = this->populations.get_flat_index({i, InfectionState::TemporaryImmunity1});
            size_t TImm2 = this->populations.get_flat_index({i, InfectionState::TemporaryImmunity2});

            size_t SIIi    = this->populations.get_flat_index({i, InfectionState::SusceptibleImprovedImmunity});
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

                flows[get_flat_flow_index<InfectionState::SusceptibleNaive, InfectionState::ExposedNaive>({i})] +=
                    dummy_SN;
                flows[get_flat_flow_index<InfectionState::SusceptiblePartialImmunity,
                                          InfectionState::ExposedPartialImmunity>({i})] += dummy_SPI;
                flows[get_flat_flow_index<InfectionState::SusceptibleImprovedImmunity,
                                          InfectionState::ExposedImprovedImmunity>({i})] += dummy_SII;
            }

            // vaccinations
            auto t_idx = SimulationDay((size_t)t);
            if (t_idx > SimulationDay(0)) {
                auto t_m1 = SimulationDay((size_t)t - 1);
                double first_vacc;
                double full_vacc;
                double booster_vacc;

                first_vacc = params.template get<DailyPartialVaccination>()[{(AgeGroup)i, t_idx}] -
                             params.template get<DailyPartialVaccination>()[{(AgeGroup)i, t_m1}];
                full_vacc = params.template get<DailyFullVaccination>()[{(AgeGroup)i, t_idx}] -
                            params.template get<DailyFullVaccination>()[{(AgeGroup)i, t_m1}];
                booster_vacc = params.template get<DailyBoosterVaccination>()[{(AgeGroup)i, t_idx}] -
                               params.template get<DailyBoosterVaccination>()[{(AgeGroup)i, t_m1}];

                double first_vaccinations =
                    (y[SNi] - flows[get_flat_flow_index<InfectionState::SusceptibleNaive, InfectionState::ExposedNaive>(
                                  {i})] <
                     first_vacc)
                        ? y[SNi] -
                              flows[get_flat_flow_index<InfectionState::SusceptibleNaive, InfectionState::ExposedNaive>(
                                  {i})]
                        : first_vacc;

                double second_vaccinations =
                    (y[SPIi] - flows[get_flat_flow_index<InfectionState::SusceptiblePartialImmunity,
                                                         InfectionState::ExposedPartialImmunity>({i})] <
                     full_vacc)
                        ? y[SPIi] - flows[get_flat_flow_index<InfectionState::SusceptiblePartialImmunity,
                                                              InfectionState::ExposedPartialImmunity>({i})]
                        : full_vacc;
                flows[get_flat_flow_index<InfectionState::SusceptibleNaive, InfectionState::TemporaryImmunity1>({i})] =
                    first_vaccinations;

                flows[get_flat_flow_index<InfectionState::SusceptiblePartialImmunity,
                                          InfectionState::TemporaryImmunity2>({i})] = second_vaccinations;

                double third_vaccinations =
                    (y[SIIi] - flows[get_flat_flow_index<InfectionState::SusceptibleImprovedImmunity,
                                                         InfectionState::ExposedImprovedImmunity>({i})] <
                     booster_vacc)
                        ? y[SIIi] - flows[get_flat_flow_index<InfectionState::SusceptibleImprovedImmunity,
                                                              InfectionState::ExposedImprovedImmunity>({i})]
                        : booster_vacc;

                flows[get_flat_flow_index<InfectionState::SusceptibleImprovedImmunity,
                                          InfectionState::TemporaryImmunity2>({i})] = third_vaccinations;
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
            // Exposed
            flows[get_flat_flow_index<InfectionState::ExposedNaive, InfectionState::InfectedNoSymptomsNaive>({i})] +=
                rateE * y[ENi];

            // // InfectedNoSymptoms
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsNaive, InfectionState::TemporaryImmunity1>(
                {i})] = params.get<RecoveredPerInfectedNoSymptoms>()[i] * rateINS * y[INSNi];
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsNaive, InfectionState::InfectedSymptomsNaive>(
                {i})] = (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * rateINS * y[INSNi];
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsNaiveConfirmed,
                                      InfectionState::InfectedSymptomsNaiveConfirmed>({i})] =
                (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * rateINS * y[INSNCi];
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsNaiveConfirmed,
                                      InfectionState::TemporaryImmunity1>({i})] =
                params.get<RecoveredPerInfectedNoSymptoms>()[i] * rateINS * y[INSNCi];

            // // InfectedSymptoms
            flows[get_flat_flow_index<InfectionState::InfectedSymptomsNaive, InfectionState::InfectedSevereNaive>(
                {i})] = params.get<SeverePerInfectedSymptoms>()[i] / params.get<TimeInfectedSymptoms>()[i] * y[ISyNi];

            flows[get_flat_flow_index<InfectionState::InfectedSymptomsNaive, InfectionState::TemporaryImmunity1>({i})] =
                (1 - params.get<SeverePerInfectedSymptoms>()[i]) / params.get<TimeInfectedSymptoms>()[i] * y[ISyNi];

            flows[get_flat_flow_index<InfectionState::InfectedSymptomsNaiveConfirmed,
                                      InfectionState::InfectedSevereNaive>({i})] =
                params.get<SeverePerInfectedSymptoms>()[i] / params.get<TimeInfectedSymptoms>()[i] * y[ISyNCi];

            flows[get_flat_flow_index<InfectionState::InfectedSymptomsNaiveConfirmed,
                                      InfectionState::TemporaryImmunity1>({i})] =
                (1 - params.get<SeverePerInfectedSymptoms>()[i]) / params.get<TimeInfectedSymptoms>()[i] * y[ISyNCi];

            // InfectedSevere
            flows[get_flat_flow_index<InfectionState::InfectedSevereNaive, InfectionState::InfectedCriticalNaive>(
                {i})] = criticalPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevNi];

            flows[get_flat_flow_index<InfectionState::InfectedSevereNaive, InfectionState::TemporaryImmunity1>({i})] =
                (1 - params.get<CriticalPerSevere>()[i]) / params.get<TimeInfectedSevere>()[i] * y[ISevNi];

            flows[get_flat_flow_index<InfectionState::InfectedSevereNaive, InfectionState::DeadNaive>({i})] =
                deathsPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevNi];
            // InfectedCritical
            flows[get_flat_flow_index<InfectionState::InfectedCriticalNaive, InfectionState::DeadNaive>({i})] =
                params.get<DeathsPerCritical>()[i] / params.get<TimeInfectedCritical>()[i] * y[ICrNi];

            flows[get_flat_flow_index<InfectionState::InfectedCriticalNaive, InfectionState::TemporaryImmunity1>({i})] =
                (1 - params.get<DeathsPerCritical>()[i]) / params.get<TimeInfectedCritical>()[i] * y[ICrNi];

            // Timm2
            flows[get_flat_flow_index<InfectionState::TemporaryImmunity1, InfectionState::SusceptiblePartialImmunity>(
                {i})] = 1 / params.get<ImmunityInterval1>()[i] * y[TImm1];

            // /**** path of partially immune ***/

            // Exposed
            flows[get_flat_flow_index<InfectionState::ExposedPartialImmunity,
                                      InfectionState::InfectedNoSymptomsPartialImmunity>({i})] += rateE * y[EPIi];

            // InfectedNoSymptoms
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsPartialImmunity,
                                      InfectionState::TemporaryImmunity2>({i})] =
                (1 - (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                         (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i])) *
                rateINS / reducTimeInfectedMild * y[INSPIi];
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsPartialImmunity,
                                      InfectionState::InfectedSymptomsPartialImmunity>({i})] =
                (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * (rateINS / reducTimeInfectedMild) * y[INSPIi];
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,
                                      InfectionState::InfectedSymptomsPartialImmunityConfirmed>({i})] =
                (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * (rateINS / reducTimeInfectedMild) * y[INSPICi];
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsPartialImmunityConfirmed,
                                      InfectionState::TemporaryImmunity2>({i})] =
                (1 - (reducInfectedSymptomsPartialImmunity / reducExposedPartialImmunity) *
                         (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i])) *
                rateINS / reducTimeInfectedMild * y[INSPICi];

            // InfectedSymptoms
            flows[get_flat_flow_index<InfectionState::InfectedSymptomsPartialImmunity,
                                      InfectionState::InfectedSeverePartialImmunity>({i})] =
                reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity *
                params.get<SeverePerInfectedSymptoms>()[i] /
                (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild) * y[ISyPIi];

            flows[get_flat_flow_index<InfectionState::InfectedSymptomsPartialImmunity,
                                      InfectionState::TemporaryImmunity2>({i})] =
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity) *
                         params.get<SeverePerInfectedSymptoms>()[i]) /
                (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild) * y[ISyPIi];

            flows[get_flat_flow_index<InfectionState::InfectedSymptomsPartialImmunityConfirmed,
                                      InfectionState::InfectedSeverePartialImmunity>({i})] =
                reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity *
                params.get<SeverePerInfectedSymptoms>()[i] /
                (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild) * y[ISyPICi];

            flows[get_flat_flow_index<InfectionState::InfectedSymptomsPartialImmunityConfirmed,
                                      InfectionState::TemporaryImmunity2>({i})] =
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSymptomsPartialImmunity) *
                         params.get<SeverePerInfectedSymptoms>()[i]) /
                (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild) * y[ISyPICi];

            // InfectedSevere
            flows[get_flat_flow_index<InfectionState::InfectedSeverePartialImmunity,
                                      InfectionState::InfectedCriticalPartialImmunity>({i})] =
                reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity *
                criticalPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevPIi];

            flows[get_flat_flow_index<InfectionState::InfectedSeverePartialImmunity,
                                      InfectionState::TemporaryImmunity2>({i})] =
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                         params.get<CriticalPerSevere>()[i]) /
                params.get<TimeInfectedSevere>()[i] * y[ISevPIi];

            flows[get_flat_flow_index<InfectionState::InfectedSeverePartialImmunity,
                                      InfectionState::DeadPartialImmunity>({i})] =
                (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                deathsPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevPIi];
            // InfectedCritical
            flows[get_flat_flow_index<InfectionState::InfectedCriticalPartialImmunity,
                                      InfectionState::DeadPartialImmunity>({i})] =
                (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                params.get<DeathsPerCritical>()[i] / params.get<TimeInfectedCritical>()[i] * y[ICrPIi];

            flows[get_flat_flow_index<InfectionState::InfectedCriticalPartialImmunity,
                                      InfectionState::TemporaryImmunity2>({i})] =
                (1 - (reducInfectedSevereCriticalDeadPartialImmunity / reducInfectedSevereCriticalDeadPartialImmunity) *
                         params.get<DeathsPerCritical>()[i]) /
                params.get<TimeInfectedCritical>()[i] * y[ICrPIi];

            // /**** path of improved immunity ***/
            // Exposed
            flows[get_flat_flow_index<InfectionState::ExposedImprovedImmunity,
                                      InfectionState::InfectedNoSymptomsImprovedImmunity>({i})] += rateE * y[EIIi];

            // InfectedNoSymptoms
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsImprovedImmunity,
                                      InfectionState::TemporaryImmunity2>({i})] =
                (1 - (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                         (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i])) *
                rateINS / reducTimeInfectedMild * y[INSIIi];
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsImprovedImmunity,
                                      InfectionState::InfectedSymptomsImprovedImmunity>({i})] =
                (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * (rateINS / reducTimeInfectedMild) * y[INSIIi];
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed,
                                      InfectionState::InfectedSymptomsImprovedImmunityConfirmed>({i})] =
                (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i]) * (rateINS / reducTimeInfectedMild) * y[INSIICi];
            flows[get_flat_flow_index<InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed,
                                      InfectionState::TemporaryImmunity2>({i})] =
                (1 - (reducInfectedSymptomsImprovedImmunity / reducExposedImprovedImmunity) *
                         (1 - params.get<RecoveredPerInfectedNoSymptoms>()[i])) *
                rateINS / reducTimeInfectedMild * y[INSIICi];

            // InfectedSymptoms
            flows[get_flat_flow_index<InfectionState::InfectedSymptomsImprovedImmunity,
                                      InfectionState::InfectedSevereImprovedImmunity>({i})] =
                reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity *
                params.get<SeverePerInfectedSymptoms>()[i] /
                (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild) * y[ISyIIi];

            flows[get_flat_flow_index<InfectionState::InfectedSymptomsImprovedImmunity,
                                      InfectionState::TemporaryImmunity2>({i})] =
                (1 - (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity) *
                         params.get<SeverePerInfectedSymptoms>()[i]) /
                (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild) * y[ISyIIi];

            flows[get_flat_flow_index<InfectionState::InfectedSymptomsImprovedImmunityConfirmed,
                                      InfectionState::InfectedSevereImprovedImmunity>({i})] =
                reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity *
                params.get<SeverePerInfectedSymptoms>()[i] /
                (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild) * y[ISyIICi];

            flows[get_flat_flow_index<InfectionState::InfectedSymptomsImprovedImmunityConfirmed,
                                      InfectionState::TemporaryImmunity2>({i})] =
                (1 - (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSymptomsImprovedImmunity) *
                         params.get<SeverePerInfectedSymptoms>()[i]) /
                (params.get<TimeInfectedSymptoms>()[i] * reducTimeInfectedMild) * y[ISyIICi];

            // InfectedSevere
            flows[get_flat_flow_index<InfectionState::InfectedSevereImprovedImmunity,
                                      InfectionState::InfectedCriticalImprovedImmunity>({i})] =
                reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity *
                criticalPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevIIi];

            flows[get_flat_flow_index<InfectionState::InfectedSevereImprovedImmunity,
                                      InfectionState::TemporaryImmunity2>({i})] =
                (1 -
                 (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                     params.get<CriticalPerSevere>()[i]) /
                params.get<TimeInfectedSevere>()[i] * y[ISevIIi];

            flows[get_flat_flow_index<InfectionState::InfectedSevereImprovedImmunity,
                                      InfectionState::DeadImprovedImmunity>({i})] =
                (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                deathsPerSevereAdjusted / params.get<TimeInfectedSevere>()[i] * y[ISevIIi];

            // InfectedCritical
            flows[get_flat_flow_index<InfectionState::InfectedCriticalImprovedImmunity,
                                      InfectionState::DeadImprovedImmunity>({i})] =
                (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                params.get<DeathsPerCritical>()[i] / params.get<TimeInfectedCritical>()[i] * y[ICrIIi];

            flows[get_flat_flow_index<InfectionState::InfectedCriticalImprovedImmunity,
                                      InfectionState::TemporaryImmunity2>({i})] =
                (1 -
                 (reducInfectedSevereCriticalDeadImprovedImmunity / reducInfectedSevereCriticalDeadImprovedImmunity) *
                     params.get<DeathsPerCritical>()[i]) /
                params.get<TimeInfectedCritical>()[i] * y[ICrIIi];

            // Timm2
            flows[get_flat_flow_index<InfectionState::TemporaryImmunity2, InfectionState::SusceptibleImprovedImmunity>(
                {i})] = 1 / params.get<ImmunityInterval2>()[i] * y[TImm2];

            // waning
            flows[get_flat_flow_index<InfectionState::SusceptiblePartialImmunity, InfectionState::SusceptibleNaive>(
                {i})] = 1 / params.get<WaningPartialImmunity>()[i] * y[SPIi];
            flows[get_flat_flow_index<InfectionState::SusceptibleImprovedImmunity,
                                      InfectionState::SusceptiblePartialImmunity>({i})] =
                1 / params.get<WaningImprovedImmunity>()[i] * y[SIIi];
        }
        auto flow_index_C =
            get_flat_flow_index<InfectionState::InfectedNoSymptomsNaive, InfectionState::TemporaryImmunity1>(
                {mio::AgeGroup(4)});
        auto flow_index_Cp =
            get_flat_flow_index<InfectionState::InfectedNoSymptomsNaiveConfirmed, InfectionState::TemporaryImmunity1>(
                {mio::AgeGroup(4)});
        auto flow_index_I =
            get_flat_flow_index<InfectionState::InfectedSymptomsNaive, InfectionState::TemporaryImmunity1>(
                {mio::AgeGroup(4)});
        auto flow_index_Ip =
            get_flat_flow_index<InfectionState::InfectedNoSymptomsNaiveConfirmed, InfectionState::TemporaryImmunity1>(
                {mio::AgeGroup(4)});
        auto flow_index_H =
            get_flat_flow_index<InfectionState::InfectedSevereNaive, InfectionState::TemporaryImmunity1>(
                {mio::AgeGroup(4)});
        auto flow_index_ICU =
            get_flat_flow_index<InfectionState::InfectedCriticalNaive, InfectionState::TemporaryImmunity1>(
                {mio::AgeGroup(4)});
        auto flow_index_Vacc =
            get_flat_flow_index<InfectionState::SusceptibleNaive, InfectionState::TemporaryImmunity1>(
                {mio::AgeGroup(4)});
        auto flow_index_Waning_ti =
            get_flat_flow_index<InfectionState::TemporaryImmunity1, InfectionState::SusceptiblePartialImmunity>(
                {mio::AgeGroup(4)});

        std::vector<double> indx_timms1;
        indx_timms1.push_back(flows[flow_index_C]);
        indx_timms1.push_back(flows[flow_index_Cp]);
        indx_timms1.push_back(flows[flow_index_I]);
        indx_timms1.push_back(flows[flow_index_Ip]);
        indx_timms1.push_back(flows[flow_index_H]);
        indx_timms1.push_back(flows[flow_index_ICU]);
        indx_timms1.push_back(flows[flow_index_Vacc]);
        indx_timms1.push_back(flows[flow_index_Waning_ti]);
        std::cout << "";
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

    // /**
    //  * @brief advance simulation to tmax.
    //  * Overwrites Simulation::advance and includes a check for dynamic NPIs in regular intervals.
    //  * @see Simulation::advance
    //  * @param tmax next stopping point of simulation
    //  * @return value at tmax
    //  */
    // Eigen::Ref<Eigen::VectorXd> advance(double tmax)
    // {
    //     auto& t_end_dyn_npis   = this->get_model().parameters.get_end_dynamic_npis();
    //     auto& dyn_npis         = this->get_model().parameters.template get<DynamicNPIsInfectedSymptoms>();
    //     auto& contact_patterns = this->get_model().parameters.template get<ContactPatterns>();

    //     double delay_lockdown;
    //     auto t        = Base::get_result().get_last_time();
    //     const auto dt = dyn_npis.get_interval().get();
    //     while (t < tmax) {

    //         auto dt_eff = std::min({dt, tmax - t, m_t_last_npi_check + dt - t});
    //         if (dt_eff >= 1.0) {
    //             dt_eff = 1.0;
    //         }

    //         if (t == 0) {
    //             //this->apply_vaccination(t); // done in init now?
    //             this->apply_b161(t);
    //         }
    //         Base::advance(t + dt_eff);
    //         if (t + 0.5 + dt_eff - std::floor(t + 0.5) >= 1) {
    //             this->apply_vaccination(t + 0.5 + dt_eff);
    //             this->apply_b161(t);
    //         }

    //         if (t > 0) {
    //             delay_lockdown = 7;
    //         }
    //         else {
    //             delay_lockdown = 0;
    //         }
    //         t = t + dt_eff;

    //         if (dyn_npis.get_thresholds().size() > 0) {
    //             if (floating_point_greater_equal(t, m_t_last_npi_check + dt)) {
    //                 if (t < t_end_dyn_npis) {
    //                     auto inf_rel = get_infections_relative(*this, t, this->get_result().get_last_value()) *
    //                                    dyn_npis.get_base_value();
    //                     auto exceeded_threshold = dyn_npis.get_max_exceeded_threshold(inf_rel);
    //                     if (exceeded_threshold != dyn_npis.get_thresholds().end() &&
    //                         (exceeded_threshold->first > m_dynamic_npi.first ||
    //                          t > double(m_dynamic_npi.second))) { //old npi was weaker or is expired

    //                         auto t_start = SimulationTime(t + delay_lockdown);
    //                         auto t_end   = t_start + SimulationTime(dyn_npis.get_duration());
    //                         this->get_model().parameters.get_start_commuter_detection() = (double)t_start;
    //                         this->get_model().parameters.get_end_commuter_detection()   = (double)t_end;
    //                         m_dynamic_npi = std::make_pair(exceeded_threshold->first, t_end);
    //                         implement_dynamic_npis(contact_patterns.get_cont_freq_mat(), exceeded_threshold->second,
    //                                                t_start, t_end, [](auto& g) {
    //                                                    return make_contact_damping_matrix(g);
    //                                                });
    //                     }
    //                 }
    //                 m_t_last_npi_check = t;
    //             }
    //         }
    //         else {
    //             m_t_last_npi_check = t;
    //         }
    //     }
    //     return this->get_result().get_last_value();
    // }

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

template <class M>
auto set_contact_mobility(M& model, Eigen::Ref<Eigen::MatrixXd> cmatrix, double /**/)
{
    auto& contact_pattern            = model.parameters.template get<ContactPatterns>();
    auto& contact_matrix             = contact_pattern.get_cont_freq_mat();
    contact_matrix[0].get_baseline() = cmatrix;
    Eigen::MatrixXd cmatrix2         = Eigen::MatrixXd::Zero(6, 6);
    contact_matrix[1].get_baseline() = cmatrix2;
    contact_matrix[2].get_baseline() = cmatrix2;
    contact_matrix[3].get_baseline() = cmatrix2;
}

template <class M>
void add_infections_from_flow(M model, double& infections_any, double& infecions_symptomatic_mobility)
{
    auto flows_step = model.get_flow_values();
    std::vector<size_t> indx_infections;
    indx_infections.reserve(18);
    std::vector<size_t> indx_symp;
    indx_symp.reserve(18);
    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); ++i) {
        indx_infections.push_back(
            model.template get_flat_flow_index<InfectionState::SusceptibleNaive, InfectionState::ExposedNaive>({i}));
        indx_infections.push_back(model.template get_flat_flow_index<InfectionState::SusceptiblePartialImmunity,
                                                                     InfectionState::ExposedPartialImmunity>({i}));
        indx_infections.push_back(model.template get_flat_flow_index<InfectionState::SusceptibleImprovedImmunity,
                                                                     InfectionState::ExposedImprovedImmunity>({i}));
        indx_symp.push_back(model.template get_flat_flow_index<InfectionState::InfectedNoSymptomsNaive,
                                                               InfectionState::InfectedSymptomsNaive>({i}));
        indx_symp.push_back(model.template get_flat_flow_index<InfectionState::InfectedNoSymptomsPartialImmunity,
                                                               InfectionState::InfectedSymptomsPartialImmunity>({i}));
        indx_symp.push_back(model.template get_flat_flow_index<InfectionState::InfectedNoSymptomsImprovedImmunity,
                                                               InfectionState::InfectedSymptomsImprovedImmunity>({i}));
    }
    auto sum_infections =
        std::accumulate(indx_infections.begin(), indx_infections.end(), 0.0, [&flows_step](double sum, int i) {
            return sum + flows_step(i);
        });
    infections_any += sum_infections;

    auto sum_symp = std::accumulate(indx_symp.begin(), indx_symp.end(), 0.0, [&flows_step](double sum, int i) {
        return sum + flows_step(i);
    });
    infecions_symptomatic_mobility += sum_symp;
}
} // namespace osecirvvs
} // namespace mio

#endif //ODESECIRVVS_MODEL_H
