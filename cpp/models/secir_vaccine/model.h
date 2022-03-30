/* 
* Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
*
* Authors: Wadim Koslow, Daniel Abele
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
#ifndef SECIR_H
#define SECIR_H

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/compartments/simulation.h"
#include "memilio/epidemiology/populations.h"
#include "secir_vaccine/infection_state.h"
#include "secir_vaccine/parameters.h"
#include "memilio/math/smoother.h"
#include "memilio/math/eigen_util.h"

namespace mio
{
namespace secirv
{
    class SecirModel : public CompartmentalModel<InfectionState, Populations<AgeGroup, InfectionState>, SecirParams>
    {
        using Base = CompartmentalModel<InfectionState, mio::Populations<AgeGroup, InfectionState>, SecirParams>;
        using Pa   = Base::ParameterSet;
        using Po   = Base::Populations;

    public:
        SecirModel(const Populations& pop, const ParameterSet& params)
            : Base(pop, params)
        {
        }

        SecirModel(int num_agegroups)
            : SecirModel(Po({AgeGroup(num_agegroups), InfectionState::Count}), Pa(AgeGroup(num_agegroups)))
        {
        }

        void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                             Eigen::Ref<Eigen::VectorXd> dydt) const override
        {
            // alpha  // percentage of asymptomatic cases
            // beta // risk of infection from the infected symptomatic patients
            // rho   // hospitalized per infectious
            // theta // icu per hospitalized
            // delta  // deaths per ICUs
            // 0: S,      1: E,     2: C,     3: I,     4: H,     5: U,     6: R,     7: D
            auto const& params   = this->parameters;
            AgeGroup n_agegroups = params.get_num_groups();

            ContactMatrixGroup const& contact_matrix = params.get<ContactPatterns>();

            auto icu_occupancy           = 0.0;
            auto test_and_trace_required = 0.0;
            for (auto i = AgeGroup(0); i < n_agegroups; ++i) {
                auto dummy_R3 = 0.5 / (params.get<IncubationTime>()[i] - params.get<SerialInterval>()[i]);
                test_and_trace_required +=
                    (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 *
                    (this->populations.get_from(pop, {i, InfectionState::Carrier}) +
                     this->populations.get_from(pop, {i, InfectionState::CarrierPartiallyImmune}));
                icu_occupancy += this->populations.get_from(pop, {i, InfectionState::ICU}) +
                                 this->populations.get_from(pop, {i, InfectionState::ICUPartiallyImmune});
            }

            for (auto i = AgeGroup(0); i < n_agegroups; i++) {

                size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
                size_t Ei = this->populations.get_flat_index({i, InfectionState::Exposed});
                size_t Ci = this->populations.get_flat_index({i, InfectionState::Carrier});
                size_t Ii = this->populations.get_flat_index({i, InfectionState::Infected});
                size_t Hi = this->populations.get_flat_index({i, InfectionState::Hospitalized});
                size_t Ui = this->populations.get_flat_index({i, InfectionState::ICU});

                size_t CDi = this->populations.get_flat_index({i, InfectionState::CarrierTested});
                size_t IDi = this->populations.get_flat_index({i, InfectionState::InfectedTested});

                size_t SVi = this->populations.get_flat_index({i, InfectionState::SusceptiblePartiallyImmune});
                size_t EVi = this->populations.get_flat_index({i, InfectionState::ExposedPartiallyImmune});
                size_t CVi = this->populations.get_flat_index({i, InfectionState::CarrierPartiallyImmune});
                size_t IVi = this->populations.get_flat_index({i, InfectionState::InfectedPartiallyImmune});
                size_t HVi = this->populations.get_flat_index({i, InfectionState::HospitalizedPartiallyImmune});
                size_t UVi = this->populations.get_flat_index({i, InfectionState::ICUPartiallyImmune});

                size_t CDVi = this->populations.get_flat_index({i, InfectionState::CarrierPartiallyImmuneTested});
                size_t IDVi = this->populations.get_flat_index({i, InfectionState::InfectedPartiallyImmuneTested});

                size_t EV2i = this->populations.get_flat_index({i, InfectionState::ExposedFullyImmune});
                size_t CV2i = this->populations.get_flat_index({i, InfectionState::CarrierFullyImmune});
                size_t IV2i = this->populations.get_flat_index({i, InfectionState::InfectedFullyImmune});
                size_t HV2i = this->populations.get_flat_index({i, InfectionState::HospitalizedFullyImmune});
                size_t UV2i = this->populations.get_flat_index({i, InfectionState::ICUFullyImmune});

                size_t CDV2i = this->populations.get_flat_index({i, InfectionState::CarrierFullyImmuneTested});
                size_t IDV2i = this->populations.get_flat_index({i, InfectionState::InfectedFullyImmuneTested});

                size_t Ri = this->populations.get_flat_index({i, InfectionState::Recovered});
                size_t Di = this->populations.get_flat_index({i, InfectionState::Dead});

                size_t ITi = this->populations.get_flat_index({i, InfectionState::InfTotal});

                dydt[Si] = 0;
                dydt[Ei] = 0;

                dydt[SVi] = 0;
                dydt[EVi] = 0;

                dydt[Ri]   = 0;
                dydt[EV2i] = 0;

                double dummy_R2 =
                    1.0 / (2 * params.get<SerialInterval>()[i] - params.get<IncubationTime>()[i]); // R2 = 1/(2SI-TINC)
                double dummy_R3 =
                    0.5 / (params.get<IncubationTime>()[i] - params.get<SerialInterval>()[i]); // R3 = 1/(2(TINC-SI))

                double risk_from_vacc    = 1.0;
                double risk_from_immune  = 1.0;
                double reduc_vacc_exp    = params.get<ExposedFactorPartiallyImmune>()[i];
                double reduc_immune_exp  = params.get<ExposedFactorFullyImmune>()[i];
                double reduc_vacc_inf    = params.get<InfectedFactorPartiallyImmune>()[i];
                double reduc_immune_inf  = params.get<InfectedFactorFullyImmune>()[i];
                double reduc_vacc_hosp   = params.get<HospitalizedFactorPartiallyImmune>()[i];
                double reduc_vacc_icu    = params.get<HospitalizedFactorPartiallyImmune>()[i];
                double reduc_vacc_dead   = params.get<HospitalizedFactorPartiallyImmune>()[i];
                double reduc_immune_hosp = params.get<HospitalizedFactorFullyImmune>()[i];
                double reduc_immune_icu  = params.get<HospitalizedFactorFullyImmune>()[i];
                double reduc_immune_dead = params.get<HospitalizedFactorFullyImmune>()[i];

                double reduc_t = params.get<InfectiousTimeFactorImmune>()[i];

                //symptomatic are less well quarantined when testing and tracing is overwhelmed so they infect more people
                auto risk_from_symptomatic = smoother_cosine(
                    test_and_trace_required, params.get<TestAndTraceCapacity>(),
                    params.get<TestAndTraceCapacity>() * 15, params.get<RiskOfInfectionFromSympomatic>()[i],
                    params.get<MaxRiskOfInfectionFromSympomatic>()[i]);

                auto risk_from_carrier = smoother_cosine(test_and_trace_required, params.get<TestAndTraceCapacity>(),
                                                         params.get<TestAndTraceCapacity>() * 2,
                                                         params.get<RelativeCarrierInfectability>()[i], 1.0);

                for (auto j = AgeGroup(0); j < n_agegroups; j++) {
                    size_t Sj = this->populations.get_flat_index({j, InfectionState::Susceptible});
                    size_t Ej = this->populations.get_flat_index({j, InfectionState::Exposed});
                    size_t Cj = this->populations.get_flat_index({j, InfectionState::Carrier});
                    size_t Ij = this->populations.get_flat_index({j, InfectionState::Infected});
                    size_t Hj = this->populations.get_flat_index({j, InfectionState::Hospitalized});
                    size_t Uj = this->populations.get_flat_index({j, InfectionState::ICU});
                    size_t Rj = this->populations.get_flat_index({j, InfectionState::Recovered});

                    size_t SVj = this->populations.get_flat_index({j, InfectionState::SusceptiblePartiallyImmune});
                    size_t EVj = this->populations.get_flat_index({j, InfectionState::ExposedPartiallyImmune});
                    size_t CVj = this->populations.get_flat_index({j, InfectionState::CarrierPartiallyImmune});
                    size_t IVj = this->populations.get_flat_index({j, InfectionState::InfectedPartiallyImmune});
                    size_t HVj = this->populations.get_flat_index({j, InfectionState::HospitalizedPartiallyImmune});
                    size_t UVj = this->populations.get_flat_index({j, InfectionState::ICUPartiallyImmune});

                    size_t CDj = this->populations.get_flat_index({j, InfectionState::CarrierTested});
                    size_t IDj = this->populations.get_flat_index({j, InfectionState::InfectedTested});

                    size_t CDVj = this->populations.get_flat_index({j, InfectionState::CarrierPartiallyImmuneTested});
                    size_t IDVj = this->populations.get_flat_index({j, InfectionState::InfectedPartiallyImmuneTested});

                    size_t EV2j = this->populations.get_flat_index({j, InfectionState::ExposedFullyImmune});
                    size_t CV2j = this->populations.get_flat_index({j, InfectionState::CarrierFullyImmune});
                    size_t IV2j = this->populations.get_flat_index({j, InfectionState::InfectedFullyImmune});
                    size_t HV2j = this->populations.get_flat_index({j, InfectionState::HospitalizedFullyImmune});
                    size_t UV2j = this->populations.get_flat_index({j, InfectionState::ICUFullyImmune});

                    size_t CDV2j = this->populations.get_flat_index({j, InfectionState::CarrierFullyImmuneTested});
                    size_t IDV2j = this->populations.get_flat_index({j, InfectionState::InfectedFullyImmuneTested});

                    // effective contact rate by contact rate between groups i and j and damping j
                    double season_val = (1 + params.get<Seasonality>() *
                                                 sin(3.141592653589793 *
                                                     (std::fmod((params.get<StartDay>() + t), 365.0) / 182.5 + 0.5)));
                    double cont_freq_eff =
                        season_val * contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                     static_cast<Eigen::Index>((size_t)j));

                    double Nj = pop[Sj] + pop[Ej] + pop[Cj] + pop[Ij] + pop[Hj] + pop[Uj] + pop[Rj] + pop[CDj] +
                                pop[IDj] + pop[SVj] + pop[EVj] + pop[CVj] + pop[IVj] + pop[HVj] + pop[UVj] + pop[CDVj] +
                                pop[IDVj] + pop[EV2j] + pop[CV2j] + pop[IV2j] + pop[HV2j] + pop[UV2j] + pop[CDV2j] +
                                pop[IDV2j]; // without died people

                    double divNj = 1.0 / Nj; // precompute 1.0/Nj

                    /*if (params.get<DynamicInfectionFromContact>()[i].size() > 0) {
                    infection_rate = params.get<DynamicInfectionFromContact>()[i][(size_t)t];
                }
                else {
                    infection_rate = params.get<InfectionProbabilityFromContact>()[i];
                }*/

                    double ext_inf_force_dummy =
                        cont_freq_eff * divNj * params.template get<InfectionProbabilityFromContact>()[(AgeGroup)i] *
                        (risk_from_carrier * (pop[Cj] + risk_from_vacc * pop[CVj] + risk_from_immune * pop[CV2j]) +
                         risk_from_symptomatic * (pop[Ij] + risk_from_vacc * pop[IVj] + risk_from_immune * pop[IV2j]));

                    double dummy_S = y[Si] * ext_inf_force_dummy;

                    double dummy_SV = y[SVi] * reduc_vacc_exp * ext_inf_force_dummy;

                    double dummy_R = y[Ri] * reduc_immune_exp * ext_inf_force_dummy;

                    dydt[Si] -= dummy_S; // -R1*(C+beta*I)*S/N0
                    dydt[Ei] += dummy_S; // R1*(C+beta*I)*S/N0-R2*E

                    dydt[SVi] -= dummy_SV; // -R1*(C+beta*I)*S/N0
                    dydt[EVi] += dummy_SV; // R1*(C+beta*I)*S/N0-R2*E

                    dydt[Ri] -= dummy_R; // -R1*(C+beta*I)*S/N0
                    dydt[EV2i] += dummy_R; // R1*(C+beta*I)*S/N0-R2*E
                }

                // ICU capacity shortage is close
                // TODO: if this is used with vaccination model, it has to be adapted if ICUCasesPerHospitalized
                // is different for different vaccination status. This is not the case here and in addition, ICUCapacity
                // is set to infinity and this functionality is deactivated, so this is OK for the moment.
                double prob_hosp2icu =
                    smoother_cosine(icu_occupancy, 0.90 * params.get<ICUCapacity>(), params.get<ICUCapacity>(),
                                    params.get<ICUCasesPerHospitalized>()[i], 0);

                double prob_hosp2dead = params.get<ICUCasesPerHospitalized>()[i] - prob_hosp2icu;

                dydt[Ei] -= dummy_R2 * y[Ei]; // only exchange of E and C done here
                dydt[Ci] = dummy_R2 * y[Ei] - ((1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 +
                                               params.get<AsymptoticCasesPerInfectious>()[i] /
                                                   params.get<InfectiousTimeAsymptomatic>()[i]) *
                                                  y[Ci];
                dydt[CDi] =
                    -((1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 +
                      params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i]) *
                    y[CDi];

                dydt[Ii] =
                    (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 * y[Ci] -
                    ((1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] +
                     params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i]) *
                        y[Ii];
                dydt[IDi] =
                    (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 * y[CDi] -
                    ((1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] +
                     params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i]) *
                        y[IDi];

                dydt[Hi] =
                    params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i] * y[Ii] +
                    params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i] * y[IDi] -
                    ((1 - params.get<ICUCasesPerHospitalized>()[i]) / params.get<HospitalizedToHomeTime>()[i] +
                     params.get<ICUCasesPerHospitalized>()[i] / params.get<HospitalizedToICUTime>()[i]) *
                        y[Hi];
                dydt[Ui] = -((1 - params.get<DeathsPerICU>()[i]) / params.get<ICUToHomeTime>()[i] +
                             params.get<DeathsPerICU>()[i] / params.get<ICUToDeathTime>()[i]) *
                           y[Ui];
                // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
                dydt[Ui] += prob_hosp2icu / params.get<HospitalizedToICUTime>()[i] * y[Hi];

                /**** path of (partially, i.e., one dose) vaccinated ***/

                dydt[EVi] -= dummy_R2 * y[EVi]; // only exchange of E and C done here
                dydt[CVi] =
                    dummy_R2 * y[EVi] -
                    ((reduc_vacc_inf / reduc_vacc_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i]) *
                         dummy_R3 +
                     (1 - (reduc_vacc_inf / reduc_vacc_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                         (params.get<InfectiousTimeAsymptomatic>()[i] * reduc_t)) *
                        y[CVi];
                dydt[CDVi] =
                    -((reduc_vacc_inf / reduc_vacc_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i]) *
                          dummy_R3 +
                      (1 - (reduc_vacc_inf / reduc_vacc_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                          (params.get<InfectiousTimeAsymptomatic>()[i] * reduc_t)) *
                    y[CDVi];
                dydt[IVi] = (reduc_vacc_inf / reduc_vacc_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i]) *
                                dummy_R3 * y[CVi] -
                            ((1 - reduc_vacc_hosp / reduc_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i]) /
                                 (params.get<InfectiousTimeMild>()[i] * reduc_t) +
                             reduc_vacc_hosp / reduc_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                                 params.get<HomeToHospitalizedTime>()[i]) *
                                y[IVi];
                dydt[IDVi] = (reduc_vacc_inf / reduc_vacc_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i]) *
                                 dummy_R3 * y[CDVi] -
                             ((1 - reduc_vacc_hosp / reduc_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i]) /
                                  (params.get<InfectiousTimeMild>()[i] * reduc_t) +
                              reduc_vacc_hosp / reduc_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                                  params.get<HomeToHospitalizedTime>()[i]) *
                                 y[IDVi];
                dydt[HVi] = reduc_vacc_hosp / reduc_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                                params.get<HomeToHospitalizedTime>()[i] * y[IVi] +
                            reduc_vacc_hosp / reduc_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                                params.get<HomeToHospitalizedTime>()[i] * y[IDVi] -
                            ((1 - reduc_vacc_icu / reduc_vacc_hosp * params.get<ICUCasesPerHospitalized>()[i]) /
                                 params.get<HospitalizedToHomeTime>()[i] +
                             reduc_vacc_icu / reduc_vacc_hosp * params.get<ICUCasesPerHospitalized>()[i] /
                                 params.get<HospitalizedToICUTime>()[i]) *
                                y[HVi];
                dydt[UVi] = -((1 - reduc_vacc_dead / reduc_vacc_icu * params.get<DeathsPerICU>()[i]) /
                                  params.get<ICUToHomeTime>()[i] +
                              reduc_vacc_dead / reduc_vacc_icu * params.get<DeathsPerICU>()[i] /
                                  params.get<ICUToDeathTime>()[i]) *
                            y[UVi];
                // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
                dydt[UVi] +=
                    reduc_vacc_icu / reduc_vacc_hosp * prob_hosp2icu / params.get<HospitalizedToICUTime>()[i] * y[HVi];

                /**** path of twice vaccinated, here called immune although reinfection is possible now ***/

                dydt[EV2i] -= dummy_R2 * y[EV2i]; // only exchange of E and C done here

                dydt[CV2i] =
                    dummy_R2 * y[EV2i] -
                    ((reduc_immune_inf / reduc_immune_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i]) *
                         dummy_R3 +
                     (1 - (reduc_immune_inf / reduc_immune_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                         (params.get<InfectiousTimeAsymptomatic>()[i] * reduc_t)) *
                        y[CV2i];
                dydt[CDV2i] = -((reduc_immune_inf / reduc_immune_exp) *
                                    (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 +
                                (1 - (reduc_immune_inf / reduc_immune_exp) *
                                         (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                                    (params.get<InfectiousTimeAsymptomatic>()[i] * reduc_t)) *
                              y[CDV2i];

                dydt[IV2i] =
                    (reduc_immune_inf / reduc_immune_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i]) *
                        dummy_R3 * y[CV2i] -
                    ((1 - reduc_immune_hosp / reduc_immune_inf * params.get<HospitalizedCasesPerInfectious>()[i]) /
                         (params.get<InfectiousTimeMild>()[i] * reduc_t) +
                     reduc_immune_hosp / reduc_immune_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                         params.get<HomeToHospitalizedTime>()[i]) *
                        y[IV2i];
                dydt[IDV2i] =
                    (reduc_immune_inf / reduc_immune_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i]) *
                        dummy_R3 * y[CDV2i] -
                    ((1 - reduc_immune_hosp / reduc_immune_inf * params.get<HospitalizedCasesPerInfectious>()[i]) /
                         (params.get<InfectiousTimeMild>()[i] * reduc_t) +
                     reduc_immune_hosp / reduc_immune_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                         params.get<HomeToHospitalizedTime>()[i]) *
                        y[IDV2i];
                dydt[HV2i] = reduc_immune_hosp / reduc_immune_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                                 params.get<HomeToHospitalizedTime>()[i] * y[IV2i] +
                             reduc_immune_hosp / reduc_immune_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                                 params.get<HomeToHospitalizedTime>()[i] * y[IDV2i] -
                             ((1 - reduc_immune_icu / reduc_immune_hosp * params.get<ICUCasesPerHospitalized>()[i]) /
                                  params.get<HospitalizedToHomeTime>()[i] +
                              reduc_immune_icu / reduc_immune_hosp * params.get<ICUCasesPerHospitalized>()[i] /
                                  params.get<HospitalizedToICUTime>()[i]) *
                                 y[HV2i];
                dydt[UV2i] = -((1 - reduc_immune_dead / reduc_immune_icu * params.get<DeathsPerICU>()[i]) /
                                   params.get<ICUToHomeTime>()[i] +
                               reduc_immune_dead / reduc_immune_icu * params.get<DeathsPerICU>()[i] /
                                   params.get<ICUToDeathTime>()[i]) *
                             y[UV2i];
                // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
                dydt[UV2i] += reduc_immune_icu / reduc_immune_hosp * prob_hosp2icu /
                              params.get<HospitalizedToICUTime>()[i] * y[HV2i];

                // compute auxiliary compartment of all past infections
                dydt[ITi] =
                    ((1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] +
                     params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i]) *
                        y[Ii] +
                    ((1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] +
                     params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i]) *
                        y[IDi] +
                    ((1 - reduc_vacc_hosp / reduc_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i]) /
                         (params.get<InfectiousTimeMild>()[i] * reduc_t) +
                     reduc_vacc_hosp / reduc_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                         params.get<HomeToHospitalizedTime>()[i]) *
                        y[IDVi] +
                    ((1 - reduc_vacc_hosp / reduc_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i]) /
                         (params.get<InfectiousTimeMild>()[i] * reduc_t) +
                     reduc_vacc_hosp / reduc_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                         params.get<HomeToHospitalizedTime>()[i]) *
                        y[IVi] +
                    ((1 - reduc_immune_hosp / reduc_immune_inf * params.get<HospitalizedCasesPerInfectious>()[i]) /
                         (params.get<InfectiousTimeMild>()[i] * reduc_t) +
                     reduc_immune_hosp / reduc_immune_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                         params.get<HomeToHospitalizedTime>()[i]) *
                        y[IDV2i] +
                    ((1 - reduc_immune_hosp / reduc_immune_inf * params.get<HospitalizedCasesPerInfectious>()[i]) /
                         (params.get<InfectiousTimeMild>()[i] * reduc_t) +
                     reduc_immune_hosp / reduc_immune_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                         params.get<HomeToHospitalizedTime>()[i]) *
                        y[IV2i];

                // recovered and deaths from all paths
                dydt[Ri] +=
                    params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i] *
                        y[Ci] +
                    params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i] *
                        y[CDi] +
                    (1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] *
                        y[Ii] +
                    (1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] *
                        y[IDi] +
                    (1 - params.get<ICUCasesPerHospitalized>()[i]) / params.get<HospitalizedToHomeTime>()[i] * y[Hi] +
                    (1 - params.get<DeathsPerICU>()[i]) / params.get<ICUToHomeTime>()[i] * y[Ui];

                dydt[Ri] +=
                    (1 - (reduc_vacc_inf / reduc_vacc_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                        (params.get<InfectiousTimeAsymptomatic>()[i] * reduc_t) * y[CVi] +
                    (1 - (reduc_vacc_inf / reduc_vacc_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                        (params.get<InfectiousTimeAsymptomatic>()[i] * reduc_t) * y[CDVi] +
                    (1 - (reduc_vacc_hosp / reduc_vacc_inf) * params.get<HospitalizedCasesPerInfectious>()[i]) /
                        (params.get<InfectiousTimeMild>()[i] * reduc_t) * y[IVi] +
                    (1 - (reduc_vacc_hosp / reduc_vacc_inf) * params.get<HospitalizedCasesPerInfectious>()[i]) /
                        (params.get<InfectiousTimeMild>()[i] * reduc_t) * y[IDVi] +
                    (1 - (reduc_vacc_icu / reduc_vacc_hosp) * params.get<ICUCasesPerHospitalized>()[i]) /
                        params.get<HospitalizedToHomeTime>()[i] * y[HVi] +
                    (1 - (reduc_vacc_dead / reduc_vacc_icu) * params.get<DeathsPerICU>()[i]) /
                        params.get<ICUToHomeTime>()[i] * y[UVi];

                dydt[Ri] +=
                    (1 - (reduc_immune_inf / reduc_immune_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                        (params.get<InfectiousTimeAsymptomatic>()[i] * reduc_t) * y[CV2i] +
                    (1 - (reduc_immune_inf / reduc_immune_exp) * (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                        (params.get<InfectiousTimeAsymptomatic>()[i] * reduc_t) * y[CDV2i] +
                    (1 - (reduc_immune_hosp / reduc_immune_inf) * params.get<HospitalizedCasesPerInfectious>()[i]) /
                        (params.get<InfectiousTimeMild>()[i] * reduc_t) * y[IV2i] +
                    (1 - (reduc_immune_hosp / reduc_immune_inf) * params.get<HospitalizedCasesPerInfectious>()[i]) /
                        (params.get<InfectiousTimeMild>()[i] * reduc_t) * y[IDV2i] +
                    (1 - (reduc_immune_icu / reduc_immune_hosp) * params.get<ICUCasesPerHospitalized>()[i]) /
                        params.get<HospitalizedToHomeTime>()[i] * y[HV2i] +
                    (1 - (reduc_immune_dead / reduc_immune_icu) * params.get<DeathsPerICU>()[i]) /
                        params.get<ICUToHomeTime>()[i] * y[UV2i];

                dydt[Di] = params.get<DeathsPerICU>()[i] / params.get<ICUToDeathTime>()[i] * y[Ui] +
                           reduc_vacc_dead / reduc_vacc_icu * params.get<DeathsPerICU>()[i] /
                               params.get<ICUToDeathTime>()[i] * y[UVi] +
                           reduc_immune_dead / reduc_immune_icu * params.get<DeathsPerICU>()[i] /
                               params.get<ICUToDeathTime>()[i] * y[UV2i];
                ;
                // add potential, additional deaths due to ICU overflow
                dydt[Di] += prob_hosp2dead / params.get<HospitalizedToICUTime>()[i] * y[Hi];

                dydt[Di] += (reduc_vacc_dead / reduc_vacc_hosp) * prob_hosp2dead /
                            params.get<HospitalizedToICUTime>()[i] * y[HVi];

                dydt[Di] += (reduc_immune_dead / reduc_immune_hosp) * prob_hosp2dead /
                            params.get<HospitalizedToICUTime>()[i] * y[HV2i];
            }
        }

        /**
        * serialize this. 
        * @see mio::serialize
        */
        template <class IOContext>
        void serialize(IOContext& io) const
        {
            auto obj = io.create_object("Secir");
            obj.add_element("Parameters", parameters);
            obj.add_element("Populations", populations);
        }

        /**
        * deserialize an object of this class.
        * @see mio::deserialize
        */
        template <class IOContext>
        static IOResult<SecirModel> deserialize(IOContext& io)
        {
            auto obj = io.expect_object("Secir");
            auto par = obj.expect_element("Parameters", Tag<ParameterSet>{});
            auto pop = obj.expect_element("Populations", Tag<Populations>{});
            return apply(
                io,
                [](auto&& par_, auto&& pop_) {
                    return SecirModel{pop_, par_};
                },
                par, pop);
        }
    };

    //forward declaration, see below.
    template <class Base = Simulation<SecirModel>>
    class SecirSimulation;

    /**
    * get percentage of infections per total population.
    * @param model the compartment model with initial values.
    * @param t current simulation time.
    * @param y current value of compartments.
    * @tparam Base simulation type that uses a secir compartment model. see SecirSimulation.
    */
    template <class Base = Simulation<SecirModel>>
    double get_infections_relative(const SecirSimulation<Base>& model, double t,
                                   const Eigen::Ref<const Eigen::VectorXd>& y);

    /**
    * specialization of compartment model simulation for secir models.
    * @tparam Base simulation type that uses a secir compartment model. default mio::SecirSimulation. For testing purposes only!
    */
    template <class Base>
    class SecirSimulation : public Base
    {
    public:
        /**
        * construct a simulation.
        * @param model the model to simulate.
        * @param t0 start time
        * @param dt time steps
        */
        SecirSimulation(SecirModel const& model, double t0 = 0., double dt = 0.1)
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
                    share_new_variant *
                        this->get_model().parameters.template get<BaseInfectiousnessB161>()[(AgeGroup)i];
                this->get_model().parameters.template get<InfectionProbabilityFromContact>()[(AgeGroup)i] =
                    new_transmission;
            }
        }

        void apply_vaccination(double t)
        {
            auto t_idx        = SimulationDay(t);
            auto& params      = this->get_model().parameters;
            size_t num_groups = (size_t)params.get_num_groups();
            auto last_value   = this->get_result().get_last_value();

            auto count = (size_t)InfectionState::Count;
            auto S     = (size_t)InfectionState::Susceptible;
            auto SV    = (size_t)InfectionState::SusceptiblePartiallyImmune;
            auto R     = (size_t)InfectionState::Recovered;

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

            auto& start_day        = this->get_model().parameters.template get<StartDay>();
            auto& start_summer     = this->get_model().parameters.template get<StartSummer>();
            auto& dyn_npis         = this->get_model().parameters.template get<DynamicNPIsInfected>();
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
                    this->apply_vaccination(t + 0.5 + dt_eff);
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
                        if (Base::get_result().get_last_time() < start_summer - start_day) {
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
    * specialization of simulate for secir models using SecirSimulation.
    * @param t0 start time.
    * @param tmax end time.
    * @param dt time step.
    * @param model secir model to simulate.
    * @param integrator optional integrator, uses rk45 if nullptr.
    */
    inline auto simulate(double t0, double tmax, double dt, const SecirModel& model,
                         std::shared_ptr<IntegratorCore> integrator = nullptr)
    {
        return mio::simulate<SecirModel, SecirSimulation<>>(t0, tmax, dt, model, integrator);
    }

    //see declaration above.
    template <class Base>
    double get_infections_relative(const SecirSimulation<Base>& sim, double /*t*/,
                                   const Eigen::Ref<const Eigen::VectorXd>& y)
    {
        double sum_inf = 0;
        for (auto i = AgeGroup(0); i < sim.get_model().parameters.get_num_groups(); ++i) {
            sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::Infected});
            sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedTested});
            sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedPartiallyImmune});
            sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedFullyImmune});
            sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedPartiallyImmuneTested});
            sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedFullyImmuneTested});
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
    * @tparam Base simulation type that uses a secir compartment model. see SecirSimulation.
    */
    template <class Base = Simulation<SecirModel>>
    auto get_migration_factors(const SecirSimulation<Base>& sim, double /*t*/,
                               const Eigen::Ref<const Eigen::VectorXd>& y)
    {
        auto& params = sim.get_model().parameters;
        //parameters as arrays
        auto& t_inc     = params.template get<IncubationTime>().array().template cast<double>();
        auto& t_ser     = params.template get<SerialInterval>().array().template cast<double>();
        auto& p_asymp   = params.template get<AsymptoticCasesPerInfectious>().array().template cast<double>();
        auto& p_inf     = params.template get<RiskOfInfectionFromSympomatic>().array().template cast<double>();
        auto& p_inf_max = params.template get<MaxRiskOfInfectionFromSympomatic>().array().template cast<double>();
        //slice of carriers
        auto y_car = slice(y, {Eigen::Index(InfectionState::Carrier), Eigen::Index(size_t(params.get_num_groups())),
                               Eigen::Index(InfectionState::Count)}) +
                     slice(y, {Eigen::Index(InfectionState::CarrierPartiallyImmune),
                               Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)}) +
                     slice(y, {Eigen::Index(InfectionState::CarrierFullyImmune),
                               Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)});

        //compute isolation, same as infection risk from main model
        auto R3                      = 0.5 / (t_inc - t_ser);
        auto test_and_trace_required = ((1 - p_asymp) * R3 * y_car.array()).sum();
        auto risk_from_symptomatic =
            smoother_cosine(test_and_trace_required, double(params.template get<TestAndTraceCapacity>()),
                            params.template get<TestAndTraceCapacity>() * 5, p_inf.matrix(), p_inf_max.matrix());

        //set factor for infected
        auto factors = Eigen::VectorXd::Ones(y.rows()).eval();
        slice(factors, {Eigen::Index(InfectionState::Infected), Eigen::Index(size_t(params.get_num_groups())),
                        Eigen::Index(InfectionState::Count)})
            .array() = risk_from_symptomatic;
        slice(factors, {Eigen::Index(InfectionState::InfectedPartiallyImmune),
                        Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)})
            .array() = risk_from_symptomatic;
        slice(factors, {Eigen::Index(InfectionState::InfectedFullyImmune),
                        Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)})
            .array() = risk_from_symptomatic;
        return factors;
    }

    template <class Base = Simulation<SecirModel>>
    auto test_commuters(SecirSimulation<Base>& sim, Eigen::Ref<Eigen::VectorXd> migrated, double time)
    {
        auto& model       = sim.get_model();
        auto nondetection = 1.0;
        if (time >= model.parameters.get_start_commuter_detection() &&
            time < model.parameters.get_end_commuter_detection()) {
            nondetection = (double)model.parameters.get_commuter_nondetection();
        }
        for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); ++i) {
            auto Ii  = model.populations.get_flat_index({i, InfectionState::Infected});
            auto IDi = model.populations.get_flat_index({i, InfectionState::InfectedTested});
            auto Ci  = model.populations.get_flat_index({i, InfectionState::Carrier});
            auto CDi = model.populations.get_flat_index({i, InfectionState::CarrierTested});

            auto IV1i  = model.populations.get_flat_index({i, InfectionState::InfectedPartiallyImmune});
            auto IDV1i = model.populations.get_flat_index({i, InfectionState::InfectedPartiallyImmuneTested});
            auto CV1i  = model.populations.get_flat_index({i, InfectionState::CarrierPartiallyImmune});
            auto CDV1i = model.populations.get_flat_index({i, InfectionState::CarrierPartiallyImmuneTested});

            auto IV2i  = model.populations.get_flat_index({i, InfectionState::InfectedFullyImmune});
            auto IDV2i = model.populations.get_flat_index({i, InfectionState::InfectedFullyImmuneTested});
            auto CV2i  = model.populations.get_flat_index({i, InfectionState::CarrierFullyImmune});
            auto CDV2i = model.populations.get_flat_index({i, InfectionState::CarrierFullyImmuneTested});

            //put detected commuters in their own compartment so they don't contribute to infections in their home node
            sim.get_result().get_last_value()[Ii] -= migrated[Ii] * (1 - nondetection);
            sim.get_result().get_last_value()[IDi] += migrated[Ii] * (1 - nondetection);
            sim.get_result().get_last_value()[Ci] -= migrated[Ci] * (1 - nondetection);
            sim.get_result().get_last_value()[CDi] += migrated[Ci] * (1 - nondetection);

            sim.get_result().get_last_value()[IV1i] -= migrated[IV1i] * (1 - nondetection);
            sim.get_result().get_last_value()[IDV1i] += migrated[IV1i] * (1 - nondetection);
            sim.get_result().get_last_value()[CV1i] -= migrated[CV1i] * (1 - nondetection);
            sim.get_result().get_last_value()[CDV1i] += migrated[CV1i] * (1 - nondetection);

            sim.get_result().get_last_value()[IV2i] -= migrated[IV2i] * (1 - nondetection);
            sim.get_result().get_last_value()[IDV2i] += migrated[IV2i] * (1 - nondetection);
            sim.get_result().get_last_value()[CV2i] -= migrated[CV2i] * (1 - nondetection);
            sim.get_result().get_last_value()[CDV2i] += migrated[CV2i] * (1 - nondetection);

            //reduce the number of commuters
            migrated[Ii] *= nondetection;
            migrated[Ci] *= nondetection;

            migrated[IV1i] *= nondetection;
            migrated[CV1i] *= nondetection;

            migrated[IV2i] *= nondetection;
            migrated[CV2i] *= nondetection;
        }
    }

} // namespace secirv
} // namespace mio

#endif
