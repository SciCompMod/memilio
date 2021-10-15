/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
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
#ifndef SECIR_VACCINATED_H
#define SECIR_VACCINATED_H

#include "epidemiology/model/compartmentalmodel.h"
#include "epidemiology/model/simulation.h"
#include "epidemiology/model/populations.h"
#include "epidemiology/secir/infection_state.h"
#include "epidemiology/secir/secir_params.h"
#include "epidemiology/math/smoother.h"
#include "epidemiology/utils/eigen_util.h"
#include "epidemiology/utils/date.h"

namespace epi
{

// Create template specializations for the age resolved
// SECIHURD model

class SecirModelV : public CompartmentalModel<Populations<AgeGroup, InfectionStateV>, SecirParams>
{
    using Base = CompartmentalModel<epi::Populations<AgeGroup, InfectionStateV>, SecirParams>;
    using Pa   = Base::ParameterSet;
    using Po   = Base::Populations;

public:
    SecirModelV(const Populations& pop, const ParameterSet& params)
        : Base(pop, params)
    {

#if !USE_DERIV_FUNC
        size_t n_agegroups = (size_t)AgeGroup::Count;
        for (size_t i = 0; i < n_agegroups; i++) {
            for (size_t j = 0; j < n_agegroups; j++) {

                // Si to Ei individually for each age group j
                this->add_flow(
                    std::make_tuple((AgeGroup)i, InfectionState::S), std::make_tuple((AgeGroup)i, InfectionState::E),
                    [i, j](Pa const& p, Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y,
                           double t) {
                        //symptomatic are less well quarantined when testing and tracing is overwhelmed so they infect more people
                        auto test_and_trace_required = (1 - params.probabilities[i].get_asymp_per_infectious()) *
                                                       dummy_R3 * Po::get_from(pop, (AgeGroup)i, InfectionState::C);
                        auto risk_from_symptomatic =
                            smoother_cosine(test_and_trace_required, params.get_test_and_trace_capacity(),
                                            params.get_test_and_trace_capacity() * 5,
                                            params.probabilities[i].get_risk_from_symptomatic(),
                                            params.probabilities[i].get_test_and_trace_max_risk_from_symptomatic());

                        // effective contact rate by contact rate between groups i and j and damping j
                        ScalarType season_val =
                            (1 + p.get_seasonality() * sin(3.141592653589793 *
                                                           (std::fmod((p.get_start_day() + t), 365.0) / 182.5 + 0.5)));
                        ScalarType cont_freq_eff =
                            season_val *
                            contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j));
                        ScalarType Nj = Po::get_from(pop, (AgeGroup)j, InfectionState::S) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::E) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::C) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::I) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::H) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::U) +
                                        Po::get_from(pop, (AgeGroup)j, InfectionState::R); // without died people
                        ScalarType divNj = 1.0 / Nj; // precompute 1.0/Nj
                        ScalarType Si    = Po::get_from(y, (AgeGroup)i, InfectionState::S);
                        ScalarType Cj    = Po::get_from(pop, (AgeGroup)j, InfectionState::C);
                        ScalarType Ij    = Po::get_from(pop, (AgeGroup)j, InfectionState::I);
                        return Si * cont_freq_eff * divNj * p.probabilities[i].get_infection_from_contact() *
                               (p.probabilities[j].get_carrier_infectability() * Cj + risk_from_symptomatic * Ij);
                    });
            }

            // Ei to Ci
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::E),
                           std::make_tuple((AgeGroup)i, InfectionState::C),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionState::E) /
                                      (2 * p.times[i].get_serialinterval() - p.times[i].get_incubation());
                           });

            // Ci to Ii
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::C),
                           std::make_tuple((AgeGroup)i, InfectionState::I),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               double dummy_R3 = 0.5 / (p.times[i].get_incubation() - p.times[i].get_serialinterval());
                               double alpha    = p.probabilities[i].get_asymp_per_infectious();
                               return ((1 - alpha) * dummy_R3) * Po::get_from(y, (AgeGroup)i, InfectionState::C);
                           });

            // Ci to Ri
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::C),
                           std::make_tuple((AgeGroup)i, InfectionState::R),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               double alpha = p.probabilities[i].get_asymp_per_infectious();
                               return (alpha / p.times[i].get_infectious_asymp()) *
                                      Po::get_from(y, (AgeGroup)i, InfectionState::C);
                           });

            // Ii to Ri
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::I),
                           std::make_tuple((AgeGroup)i, InfectionState::R),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionState::I) *
                                      (1 - p.probabilities[i].get_hospitalized_per_infectious()) /
                                      p.times[i].get_infectious_mild();
                           });

            // Ii to Hi
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionState::I), std::make_tuple((AgeGroup)i, InfectionState::H),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y,
                    double /*t*/) {
                    return Po::get_from(y, (AgeGroup)i, InfectionState::I) *
                           p.probabilities[i].get_hospitalized_per_infectious() / p.times[i].get_home_to_hospitalized();
                });

            // Hi to Ui
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::H),
                           std::make_tuple((AgeGroup)i, InfectionState::U),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               ScalarType icu_occupancy = 0;
                               for (size_t j = 0; j < (size_t)AgeGroup::Count; ++j) {
                                   icu_occupancy += Po::get_from(y, (AgeGroup)j, InfectionState::U);
                               }

                               ScalarType prob_hosp2icu =
                                   smoother_cosine(icu_occupancy, 0.90 * p.get_icu_capacity(), p.get_icu_capacity(),
                                                   p.probabilities[i].get_icu_per_hospitalized(), 0);

                               return Po::get_from(y, (AgeGroup)i, InfectionState::H) * prob_hosp2icu /
                                      p.times[i].get_hospitalized_to_icu();
                           });

            // Hi to Di
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionState::H), std::make_tuple((AgeGroup)i, InfectionState::D),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y,
                    double /*t*/) {
                    ScalarType icu_occupancy = 0;
                    for (size_t j = 0; j < (size_t)AgeGroup::Count; ++j) {
                        icu_occupancy += Po::get_from(y, (AgeGroup)j, InfectionState::U);
                    }

                    ScalarType prob_hosp2icu =
                        smoother_cosine(icu_occupancy, 0.90 * p.get_icu_capacity(), p.get_icu_capacity(),
                                        p.probabilities[i].get_icu_per_hospitalized(), 0);
                    ScalarType prob_hosp2dead = p.probabilities[i].get_icu_per_hospitalized() - prob_hosp2icu;

                    return Po::get_from(y, (AgeGroup)i, InfectionState::H) * prob_hosp2dead /
                           p.times[i].get_hospitalized_to_icu();
                });

            // Hi to Ri
            this->add_flow(
                std::make_tuple((AgeGroup)i, InfectionState::H), std::make_tuple((AgeGroup)i, InfectionState::R),
                [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/, Eigen::Ref<const Eigen::VectorXd> y,
                    double /*t*/) {
                    return Po::get_from(y, (AgeGroup)i, InfectionState::H) *
                           (1 - p.probabilities[i].get_icu_per_hospitalized()) / p.times[i].get_hospitalized_to_home();
                });

            // Ui to Ri
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::U),
                           std::make_tuple((AgeGroup)i, InfectionState::R),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionState::U) *
                                      (1 - p.probabilities[i].get_dead_per_icu()) / p.times[i].get_icu_to_home();
                           });

            // Ui to Di
            this->add_flow(std::make_tuple((AgeGroup)i, InfectionState::U),
                           std::make_tuple((AgeGroup)i, InfectionState::D),
                           [i](Pa const& p, Eigen::Ref<const Eigen::VectorXd> /*pop*/,
                               Eigen::Ref<const Eigen::VectorXd> y, double /*t*/) {
                               return Po::get_from(y, (AgeGroup)i, InfectionState::U) *
                                      p.probabilities[i].get_dead_per_icu() / p.times[i].get_icu_to_dead();
                           });
        }
#endif
    }

    SecirModelV(int num_agegroups)
        : SecirModelV(Po({AgeGroup(num_agegroups), InfectionStateV::Count}), Pa(AgeGroup(num_agegroups)))
    {
    }

#if USE_DERIV_FUNC

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

        ContactMatrixGroup const& contact_matrix = params.get<epi::ContactPatterns>();

        auto icu_occupancy           = 0.0;
        auto test_and_trace_required = 0.0;
        for (auto i = AgeGroup(0); i < n_agegroups; ++i) {
            auto dummy_R3 = 0.5 / (params.get<IncubationTime>()[i] - params.get<SerialInterval>()[i]);
            test_and_trace_required += (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 *
                                       (this->populations.get_from(pop, {i, InfectionStateV::Carrier}) +
                                        this->populations.get_from(pop, {i, InfectionStateV::CarrierV1}));
            icu_occupancy += this->populations.get_from(pop, {i, InfectionStateV::ICU}) +
                             this->populations.get_from(pop, {i, InfectionStateV::ICUV1});
        }

        for (auto i = AgeGroup(0); i < n_agegroups; i++) {

            size_t Si = this->populations.get_flat_index({i, InfectionStateV::Susceptible});
            size_t Ei = this->populations.get_flat_index({i, InfectionStateV::Exposed});
            size_t Ci = this->populations.get_flat_index({i, InfectionStateV::Carrier});
            size_t Ii = this->populations.get_flat_index({i, InfectionStateV::Infected});
            size_t Hi = this->populations.get_flat_index({i, InfectionStateV::Hospitalized});
            size_t Ui = this->populations.get_flat_index({i, InfectionStateV::ICU});
            size_t Ri = this->populations.get_flat_index({i, InfectionStateV::Recovered});
            size_t Di = this->populations.get_flat_index({i, InfectionStateV::Dead});

            size_t SVi = this->populations.get_flat_index({i, InfectionStateV::SusceptibleV1});
            size_t EVi = this->populations.get_flat_index({i, InfectionStateV::ExposedV1});
            size_t CVi = this->populations.get_flat_index({i, InfectionStateV::CarrierV1});
            size_t IVi = this->populations.get_flat_index({i, InfectionStateV::InfectedV1});
            size_t HVi = this->populations.get_flat_index({i, InfectionStateV::HospitalizedV1});
            size_t UVi = this->populations.get_flat_index({i, InfectionStateV::ICUV1});

            size_t CDi = this->populations.get_flat_index({i, InfectionStateV::CarrierT});
            size_t IDi = this->populations.get_flat_index({i, InfectionStateV::InfectedT});

            size_t CDVi = this->populations.get_flat_index({i, InfectionStateV::CarrierTV1});
            size_t IDVi = this->populations.get_flat_index({i, InfectionStateV::InfectedTV1});

            dydt[Si] = 0;
            dydt[Ei] = 0;

            dydt[SVi] = 0;
            dydt[EVi] = 0;

            double dummy_R2 =
                1.0 / (2 * params.get<SerialInterval>()[i] - params.get<IncubationTime>()[i]); // R2 = 1/(2SI-TINC)
            double dummy_R3 =
                0.5 / (params.get<IncubationTime>()[i] - params.get<SerialInterval>()[i]); // R3 = 1/(2(TINC-SI))

            double risk_from_vacc = 1.0;
            double risk_vacc_inf  = 0.5;
            double reduc_inf_hosp = 0.25; // TODO
            double reduc_hosp_icu = 0.4; // TODO
            double reduc_icu_dead = 0.4; // TODO

            //symptomatic are less well quarantined when testing and tracing is overwhelmed so they infect more people
            auto risk_from_symptomatic = smoother_cosine(
                test_and_trace_required, params.get<TestAndTraceCapacity>(), params.get<TestAndTraceCapacity>() * 15,
                params.get<RiskOfInfectionFromSympomatic>()[i], params.get<MaxRiskOfInfectionFromSympomatic>()[i]);

            auto risk_from_carrier = smoother_cosine(test_and_trace_required, params.get<TestAndTraceCapacity>(),
                                                     params.get<TestAndTraceCapacity>() * 2,
                                                     params.get<RelativeCarrierInfectability>()[i], 1.0);

            for (auto j = AgeGroup(0); j < n_agegroups; j++) {
                size_t Sj = this->populations.get_flat_index({j, InfectionStateV::Susceptible});
                size_t Ej = this->populations.get_flat_index({j, InfectionStateV::Exposed});
                size_t Cj = this->populations.get_flat_index({j, InfectionStateV::Carrier});
                size_t Ij = this->populations.get_flat_index({j, InfectionStateV::Infected});
                size_t Hj = this->populations.get_flat_index({j, InfectionStateV::Hospitalized});
                size_t Uj = this->populations.get_flat_index({j, InfectionStateV::ICU});
                size_t Rj = this->populations.get_flat_index({j, InfectionStateV::Recovered});

                size_t SVj = this->populations.get_flat_index({j, InfectionStateV::SusceptibleV1});
                size_t EVj = this->populations.get_flat_index({j, InfectionStateV::ExposedV1});
                size_t CVj = this->populations.get_flat_index({j, InfectionStateV::CarrierV1});
                size_t IVj = this->populations.get_flat_index({j, InfectionStateV::InfectedV1});
                size_t HVj = this->populations.get_flat_index({j, InfectionStateV::HospitalizedV1});
                size_t UVj = this->populations.get_flat_index({j, InfectionStateV::ICUV1});

                size_t CDj = this->populations.get_flat_index({j, InfectionStateV::CarrierT});
                size_t IDj = this->populations.get_flat_index({j, InfectionStateV::InfectedT});

                size_t CDVj = this->populations.get_flat_index({j, InfectionStateV::CarrierTV1});
                size_t IDVj = this->populations.get_flat_index({j, InfectionStateV::InfectedTV1});

                // effective contact rate by contact rate between groups i and j and damping j
                double season_val = (1 + params.get<epi::Seasonality>() *
                                             sin(3.141592653589793 *
                                                 (std::fmod((params.get<epi::StartDay>() + t), 365.0) / 182.5 + 0.5)));
                double cont_freq_eff =
                    season_val * contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                 static_cast<Eigen::Index>((size_t)j));

                double Nj = pop[Sj] + pop[Ej] + pop[Cj] + pop[Ij] + pop[Hj] + pop[Uj] + pop[Rj] + pop[CDj] + pop[IDj] +
                            pop[SVj] + pop[EVj] + pop[CVj] + pop[IVj] + pop[HVj] + pop[UVj] + pop[CDVj] +
                            pop[IDVj]; // without died people

                double divNj = 1.0 / Nj; // precompute 1.0/Nj

                /*if (params.get<DynamicInfectionFromContact>()[i].size() > 0) {
                    infection_rate = params.get<DynamicInfectionFromContact>()[i][(size_t)t];
                }
                else {
                    infection_rate = params.get<InfectionProbabilityFromContact>()[i];
                }*/

                auto start_day   = params.template get<StartDay>();
                auto b161_growth = (start_day - get_day_in_year(Date(2021, 6, 6))) * 0.1666667;
                // 2 equal to the share of the delta variant on June 6
                double share_new_variant = std::min(1.0, pow(2, t * 0.1666667 + b161_growth) * 0.01);
                double transmission_prob = (1 - share_new_variant) * params.template get<BaseInfB117>()[(AgeGroup)i] +
                                           share_new_variant * params.template get<BaseInfB161>()[(AgeGroup)i];
                unused(transmission_prob);

                double dummy_S = y[Si] * cont_freq_eff * divNj *
                                 params.template get<InfectionProbabilityFromContact>()[(AgeGroup)i] *
                                 (risk_from_carrier * (pop[Cj] + risk_from_vacc * pop[CVj]) +
                                  risk_from_symptomatic * (pop[Ij] + risk_from_vacc * pop[IVj]));

                double dummy_SV = y[SVi] * cont_freq_eff * divNj *
                                  params.template get<InfectionProbabilityFromContact>()[(AgeGroup)i] * risk_vacc_inf *
                                  (risk_from_carrier * (pop[Cj] + risk_from_vacc * pop[CVj]) +
                                   risk_from_symptomatic * (pop[Ij] + risk_from_vacc * pop[IVj]));

                dydt[Si] -= dummy_S; // -R1*(C+beta*I)*S/N0
                dydt[Ei] += dummy_S; // R1*(C+beta*I)*S/N0-R2*E

                dydt[SVi] -= dummy_SV; // -R1*(C+beta*I)*S/N0
                dydt[EVi] += dummy_SV; // R1*(C+beta*I)*S/N0-R2*E
            }

            std::cout << "Infection from Contact at time " << t << ": "
                      << params.template get<InfectionProbabilityFromContact>()[AgeGroup(0)] << std::endl;
            // ICU capacity shortage is close
            double prob_hosp2icu =
                smoother_cosine(icu_occupancy, 0.90 * params.get<epi::ICUCapacity>(), params.get<epi::ICUCapacity>(),
                                params.get<ICUCasesPerHospitalized>()[i], 0);

            double prob_hosp2dead = params.get<ICUCasesPerHospitalized>()[i] - prob_hosp2icu;

            dydt[Ei] -= dummy_R2 * y[Ei]; // only exchange of E and C done here
            dydt[Ci] = dummy_R2 * y[Ei] -
                       ((1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 +
                        params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i]) *
                           y[Ci];
            dydt[CDi] = -((1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 +
                          params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i]) *
                        y[CDi];

            dydt[Ii] = (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 * y[Ci] -
                       ((1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] +
                        params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i]) *
                           y[Ii];
            dydt[IDi] = (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 * y[CDi] -
                        ((1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] +
                         params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i]) *
                            y[IDi];

            dydt[Hi] =
                params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i] * y[Ii] +
                params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i] * y[IDi] -
                ((1 - params.get<ICUCasesPerHospitalized>()[i]) / params.get<HospitalizedToHomeTime>()[i] +
                 params.get<ICUCasesPerHospitalized>()[i] / params.get<HospitalizedToICUTime>()[i]) *
                    y[Hi];
            dydt[Ui] = -((1 - params.get<DeathsPerHospitalized>()[i]) / params.get<ICUToHomeTime>()[i] +
                         params.get<DeathsPerHospitalized>()[i] / params.get<ICUToDeathTime>()[i]) *
                       y[Ui];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[Ui] += prob_hosp2icu / params.get<HospitalizedToICUTime>()[i] * y[Hi];

            /**** path of vaccinated ***/

            dydt[EVi] -= dummy_R2 * y[EVi]; // only exchange of E and C done here
            dydt[CVi] = dummy_R2 * y[EVi] -
                        ((1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 +
                         params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i]) *
                            y[CVi];
            dydt[CDVi] =
                -((1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 +
                  params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i]) *
                y[CDVi];
            dydt[IVi] = (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 * y[CVi] -
                        ((1 - reduc_inf_hosp / risk_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i]) /
                             params.get<InfectiousTimeMild>()[i] +
                         reduc_inf_hosp / risk_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                             params.get<HomeToHospitalizedTime>()[i]) *
                            y[IVi];
            dydt[IDVi] = (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 * y[CDVi] -
                         ((1 - reduc_inf_hosp / risk_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i]) /
                              params.get<InfectiousTimeMild>()[i] +
                          reduc_inf_hosp / risk_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                              params.get<HomeToHospitalizedTime>()[i]) *
                             y[IDVi];
            dydt[HVi] = reduc_inf_hosp / risk_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                            params.get<HomeToHospitalizedTime>()[i] * y[IVi] +
                        reduc_inf_hosp / risk_vacc_inf * params.get<HospitalizedCasesPerInfectious>()[i] /
                            params.get<HomeToHospitalizedTime>()[i] * y[IDVi] -
                        ((1 - reduc_hosp_icu / reduc_inf_hosp * params.get<ICUCasesPerHospitalized>()[i]) /
                             params.get<HospitalizedToHomeTime>()[i] +
                         reduc_hosp_icu / reduc_inf_hosp * params.get<ICUCasesPerHospitalized>()[i] /
                             params.get<HospitalizedToICUTime>()[i]) *
                            y[HVi];
            dydt[UVi] = -((1 - reduc_icu_dead / reduc_hosp_icu * params.get<DeathsPerHospitalized>()[i]) /
                              params.get<ICUToHomeTime>()[i] +
                          reduc_icu_dead / reduc_hosp_icu * params.get<DeathsPerHospitalized>()[i] /
                              params.get<ICUToDeathTime>()[i]) *
                        y[UVi];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[UVi] +=
                reduc_hosp_icu / reduc_inf_hosp * prob_hosp2icu / params.get<HospitalizedToICUTime>()[i] * y[HVi];

            dydt[Ri] =
                params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i] * y[Ci] +
                params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i] * y[CDi] +
                (1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] * y[Ii] +
                (1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] * y[IDi] +
                params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i] * y[CVi] +
                params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i] * y[CDVi] +
                (1 - (reduc_inf_hosp / risk_vacc_inf) * params.get<HospitalizedCasesPerInfectious>()[i]) /
                    params.get<InfectiousTimeMild>()[i] * y[IVi] +
                (1 - (reduc_inf_hosp / risk_vacc_inf) * params.get<HospitalizedCasesPerInfectious>()[i]) /
                    params.get<InfectiousTimeMild>()[i] * y[IDVi] +
                (1 - params.get<ICUCasesPerHospitalized>()[i]) / params.get<HospitalizedToHomeTime>()[i] * y[Hi] +
                (1 - params.get<DeathsPerHospitalized>()[i]) / params.get<ICUToHomeTime>()[i] * y[Ui] +
                (1 - (reduc_hosp_icu / reduc_inf_hosp) * params.get<ICUCasesPerHospitalized>()[i]) /
                    params.get<HospitalizedToHomeTime>()[i] * y[HVi] +
                (1 - (reduc_icu_dead / reduc_hosp_icu) * params.get<DeathsPerHospitalized>()[i]) /
                    params.get<ICUToHomeTime>()[i] * y[UVi];

            dydt[Di] = params.get<DeathsPerHospitalized>()[i] / params.get<ICUToDeathTime>()[i] * y[Ui] +
                       reduc_icu_dead / reduc_hosp_icu * params.get<DeathsPerHospitalized>()[i] /
                           params.get<ICUToDeathTime>()[i] * y[UVi];
            // add potential, additional deaths due to ICU overflow
            dydt[Di] += prob_hosp2dead / params.get<HospitalizedToICUTime>()[i] * y[Hi];
            dydt[Di] += (params.get<ICUCasesPerHospitalized>()[i] - reduc_hosp_icu / reduc_inf_hosp * prob_hosp2icu) /
                        params.get<HospitalizedToICUTime>()[i] * y[HVi];
        }
    }

#endif // USE_DERIV_FUNC

    /**
     * serialize this. 
     * @see epi::serialize
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
     * @see epi::deserialize
     */
    template <class IOContext>
    static IOResult<SecirModelV> deserialize(IOContext& io)
    {
        auto obj = io.expect_object("Secir");
        auto par = obj.expect_element("Parameters", Tag<ParameterSet>{});
        auto pop = obj.expect_element("Populations", Tag<Populations>{});
        return apply(
            io,
            [](auto&& par_, auto&& pop_) {
                return SecirModelV{pop_, par_};
            },
            par, pop);
    }
};

//forward declaration, see below.
template <class Base = Simulation<SecirModelV>>
class SecirSimulationV;

/**
 * get percentage of infections per total population.
 * @param model the compartment model with initial values.
 * @param t current simulation time.
 * @param y current value of compartments.
 * @tparam Base simulation type that uses a secir compartment model. see SecirSimulationV.
 */
template <class Base = Simulation<SecirModelV>>
double get_infections_relative(const SecirSimulationV<Base>& model, double t,
                               const Eigen::Ref<const Eigen::VectorXd>& y);

/**
 * specialization of compartment model simulation for secir models.
 * @tparam Base simulation type that uses a secir compartment model. default epi::SecirSimulationV. For testing purposes only!
 */
template <class Base>
class SecirSimulationV : public Base
{
public:
    /**
     * construct a simulation.
     * @param model the model to simulate.
     * @param t0 start time
     * @param dt time steps
     */
    SecirSimulationV(SecirModelV const& model, double t0 = 0., double dt = 0.1)
        : Base(model, t0, dt)
        , m_t_last_npi_check(t0)
    {
    }

    void apply_b161(double t)
    {

        size_t num_groups = this->get_model().parameters.get_num_groups();
        for (size_t i = 0; i < num_groups; ++i) {
            double share_new_variant = std::min(1.0, pow(2, t / 7) / 100.0);
            double new_transmission =
                (1 - share_new_variant) * this->get_model().parameters.template get<BaseInfB117>()[(AgeGroup)i] +
                share_new_variant * this->get_model().parameters.template get<BaseInfB161>()[(AgeGroup)i];
            this->get_model().parameters.template get<InfectionProbabilityFromContact>()[(AgeGroup)i] =
                new_transmission;
        }
    }

    void apply_vaccination()
    {
        double threshhold = 0.8;
        auto& m_params    = this->get_model().parameters;
        auto& population  = this->get_model().populations;
        size_t num_groups = m_params.get_num_groups();
        auto last_value   = this->get_result().get_last_value();

        std::vector<bool> saturated(num_groups, false);

        double leftover = 0;
        size_t count    = 0;
        for (size_t i = 2; i < num_groups; ++i) {
            std::vector<double> daily_first_vaccinations = m_params.template get<DailyFirstVaccination>()[(AgeGroup)i];
            std::vector<double> daily_full_vaccinations  = m_params.template get<DailyFullVaccination>()[(AgeGroup)i];
            double new_vaccinated_full =
                daily_first_vaccinations[daily_first_vaccinations.size() -
                                         m_params.template get<VaccinationGap>()[(AgeGroup)i] -
                                         m_params.template get<DaysUntilEffective>()[(AgeGroup)i]];

            /*std::cout << "size of ag " << i << ": " << daily_first_vaccinations.size() << ", value: "
                      << daily_first_vaccinations[daily_first_vaccinations.size() -
                                                  m_params.template get<VaccinationGap>()[(AgeGroup)i] -
                                                  m_params.template get<DaysUntilEffective>()[(AgeGroup)i]]
                      << std::endl;*/
            double new_vaccinated_first =
                daily_first_vaccinations[daily_first_vaccinations.size() -
                                         m_params.template get<DaysUntilEffective>()[(AgeGroup)i]] +
                daily_full_vaccinations[daily_first_vaccinations.size() -
                                        m_params.template get<DaysUntilEffective>()[(AgeGroup)i]] -
                new_vaccinated_full;
            /*if (new_vaccinated_first < 0) {
                new_vaccinated_first = 0;
            }*/
            if (last_value((size_t)InfectionStateV::Count * i + (size_t)InfectionStateV::Susceptible) -
                    new_vaccinated_first >
                (1 - threshhold) * population.get_group_total((AgeGroup)i)) {
                last_value((size_t)InfectionStateV::Count * i + (size_t)InfectionStateV::Susceptible) -=
                    new_vaccinated_first;
                last_value((size_t)InfectionStateV::Count * i + (size_t)InfectionStateV::SusceptibleV1) +=
                    new_vaccinated_first;
                m_params.template get<DailyFirstVaccination>()[(AgeGroup)i].push_back(new_vaccinated_first);
            }
            else {
                m_params.template get<DailyFirstVaccination>()[(AgeGroup)i].push_back(0);
                saturated[i] = true;
                count++;
                leftover += new_vaccinated_first;
            }

            if (last_value((size_t)InfectionStateV::Count * i + (size_t)InfectionStateV::SusceptibleV1) -
                    new_vaccinated_full >
                0) {
                last_value((size_t)InfectionStateV::Count * i + (size_t)InfectionStateV::SusceptibleV1) -=
                    new_vaccinated_full;
                last_value((size_t)InfectionStateV::Count * i + (size_t)InfectionStateV::Recovered) +=
                    new_vaccinated_full;

                m_params.template get<DailyFullVaccination>()[(AgeGroup)i].push_back(new_vaccinated_full);
            }
            else {
                leftover += new_vaccinated_full;
                m_params.template get<DailyFullVaccination>()[(AgeGroup)i].push_back(0);
            }
        }
        for (size_t i = 2; i < num_groups; ++i) {
            if (!saturated[i] && count > 0) {
                last_value((size_t)InfectionStateV::Count * i + (size_t)InfectionStateV::Susceptible) -=
                    leftover / count;
                last_value((size_t)InfectionStateV::Count * i + (size_t)InfectionStateV::SusceptibleV1) +=
                    leftover / count;
                m_params.template get<DailyFirstVaccination>()[(
                    AgeGroup)i][m_params.template get<DailyFirstVaccination>()[(AgeGroup)i].size() - 1] +=
                    leftover / count;
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
            Base::advance(t + dt_eff);
            if (t + dt_eff - std::floor(t) >= 1) {
                this->apply_vaccination();
                //this->apply_b161(t);
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
                            auto t_end    = epi::SimulationTime(t + double(dyn_npis.get_duration()) + delay_lockdown);
                            m_dynamic_npi = std::make_pair(exceeded_threshold->first, t_end);
                            epi::implement_dynamic_npis(contact_patterns.get_cont_freq_mat(),
                                                        exceeded_threshold->second, SimulationTime(t + delay_lockdown),
                                                        t_end, [this](auto& g) {
                                                            return epi::make_contact_damping_matrix(g);
                                                        });
                        }
                    }
                }
            }
            m_t_last_npi_check = t;
        }
        return this->get_result().get_last_value();
    }

private:
    double m_t_last_npi_check;
    std::pair<double, SimulationTime> m_dynamic_npi = {-std::numeric_limits<double>::max(), epi::SimulationTime(0)};
};

/**
 * specialization of simulate for secir models using SecirSimulationV.
 * @param t0 start time.
 * @param tmax end time.
 * @param dt time step.
 * @param model secir model to simulate.
 * @param integrator optional integrator, uses rk45 if nullptr.
 */
inline auto simulate(double t0, double tmax, double dt, const SecirModelV& model,
                     std::shared_ptr<IntegratorCore> integrator = nullptr)
{
    return simulate<SecirModelV, SecirSimulationV<>>(t0, tmax, dt, model, integrator);
}

//see declaration above.
template <class Base>
double get_infections_relative(const SecirSimulationV<Base>& sim, double /*t*/,
                               const Eigen::Ref<const Eigen::VectorXd>& y)
{
    double sum_inf = 0;
    for (auto i = AgeGroup(0); i < sim.get_model().parameters.get_num_groups(); ++i) {
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionStateV::Infected});
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionStateV::InfectedT});
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionStateV::InfectedV1});
        //sum_inf += sim.get_model().populations.get_from(y, {i, InfectionStateV::InfectedV2});
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionStateV::InfectedTV1});
        //sum_inf += sim.get_model().populations.get_from(y, {i, InfectionStateV::InfectedTV2});
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
 * @tparam Base simulation type that uses a secir compartment model. see SecirSimulationV.
 */
template <class Base = Simulation<SecirModelV>>
auto get_migration_factors(const SecirSimulationV<Base>& sim, double /*t*/, const Eigen::Ref<const Eigen::VectorXd>& y)
{
    auto& params = sim.get_model().parameters;
    //parameters as arrays
    auto& t_inc     = params.template get<IncubationTime>().array().template cast<double>();
    auto& t_ser     = params.template get<SerialInterval>().array().template cast<double>();
    auto& p_asymp   = params.template get<AsymptoticCasesPerInfectious>().array().template cast<double>();
    auto& p_inf     = params.template get<RiskOfInfectionFromSympomatic>().array().template cast<double>();
    auto& p_inf_max = params.template get<MaxRiskOfInfectionFromSympomatic>().array().template cast<double>();
    //slice of carriers
    auto y_car = slice(y, {Eigen::Index(InfectionStateV::Carrier), Eigen::Index(size_t(params.get_num_groups())),
                           Eigen::Index(InfectionStateV::Count)}) +
                 slice(y, {Eigen::Index(InfectionStateV::CarrierV1), Eigen::Index(size_t(params.get_num_groups())),
                           Eigen::Index(InfectionStateV::Count)});

    //compute isolation, same as infection risk from main model
    auto R3                      = 0.5 / (t_inc - t_ser);
    auto test_and_trace_required = ((1 - p_asymp) * R3 * y_car.array()).sum();
    auto risk_from_symptomatic =
        smoother_cosine(test_and_trace_required, double(params.template get<TestAndTraceCapacity>()),
                        params.template get<TestAndTraceCapacity>() * 5, p_inf.matrix(), p_inf_max.matrix());

    //set factor for infected
    auto factors = Eigen::VectorXd::Ones(y.rows()).eval();
    slice(factors, {Eigen::Index(InfectionStateV::Infected), Eigen::Index(size_t(params.get_num_groups())),
                    Eigen::Index(InfectionStateV::Count)})
        .array() = risk_from_symptomatic;
    slice(factors, {Eigen::Index(InfectionStateV::InfectedV1), Eigen::Index(size_t(params.get_num_groups())),
                    Eigen::Index(InfectionStateV::Count)})
        .array() = risk_from_symptomatic;
    return factors;
}

/*void set_incidence_rate(SecirModelV& model, double incidence_rate)
{
    size_t n_agegroups = (size_t)model.parameters.get_num_groups();

    std::vector<size_t> indices;

    auto params = model.parameters;
    for (size_t j = 0; j < n_agegroups; j++) {

        double Nj = model.populations.get_group_total((AgeGroup)j);

        double dummy_I                                              = incidence_rate * Nj;
        model.populations[{(AgeGroup)j, InfectionStateV::Infected}] = dummy_I;

        double T_E_C = 2 * params.get<SerialInterval>()[(AgeGroup)j] -
                       params.get<IncubationTime>()[(AgeGroup)j]; // (T_E^C)^{-1} = R2 = 1/(2SI-TINC)
        double T_C_I = 2 * (params.get<IncubationTime>()[(AgeGroup)j] -
                            params.get<SerialInterval>()[(AgeGroup)j]); // (T_C^I)^{-1} = R3 = 1/(2(TINC-SI))

        double dummy_C = T_C_I *
                         (1 -
                          params.get<AsymptoticCasesPerInfectious>()[(AgeGroup)j] *
                              ((1 - params.get<HospitalizedCasesPerInfectious>()[(AgeGroup)j])) /
                              params.get<InfectiousTimeMild>()[(AgeGroup)j] +
                          params.get<HospitalizedCasesPerInfectious>()[(AgeGroup)j] /
                              params.get<HomeToHospitalizedTime>()[(AgeGroup)j]) *
                         dummy_I;
        model.populations[{(AgeGroup)j, InfectionStateV::Carrier}] = dummy_C;

        double dummy_E = T_E_C *
                         ((1 - params.get<AsymptoticCasesPerInfectious>()[(AgeGroup)j]) / T_C_I +
                          params.get<AsymptoticCasesPerInfectious>()[(AgeGroup)j] /
                              params.get<InfectiousTimeAsymptomatic>()[(AgeGroup)j]) *
                         dummy_C;
        model.populations[{(AgeGroup)j, InfectionStateV::Exposed}] = dummy_E;
    }
}*/

template <class Base = Simulation<SecirModelV>>
auto test_commuters(SecirSimulationV<Base>& sim, Eigen::Ref<Eigen::VectorXd> migrated, double time)
{
    auto& model       = sim.get_model();
    auto nondetection = 1.0;
    if (time >= model.parameters.get_start_commuter_nondetection() &&
        time < model.parameters.get_end_commuter_nondetection()) {
        nondetection = (double)model.parameters.get_commuter_nondetection();
    }
    for (auto i = AgeGroup(0); i < model.parameters.get_num_groups(); ++i) {
        auto Ii  = model.populations.get_flat_index({i, InfectionStateV::Infected});
        auto IDi = model.populations.get_flat_index({i, InfectionStateV::InfectedT});
        auto Ci  = model.populations.get_flat_index({i, InfectionStateV::Carrier});
        auto CDi = model.populations.get_flat_index({i, InfectionStateV::CarrierT});

        auto IV1i  = model.populations.get_flat_index({i, InfectionStateV::InfectedV1});
        auto IDV1i = model.populations.get_flat_index({i, InfectionStateV::InfectedTV1});
        auto CV1i  = model.populations.get_flat_index({i, InfectionStateV::CarrierV1});
        auto CDV1i = model.populations.get_flat_index({i, InfectionStateV::CarrierTV1});

        /*auto IV2i  = model.populations.get_flat_index({i, size_t(InfectionStateV::InfectedV2)});
        auto IDV2i = model.populations.get_flat_index({i, size_t(InfectionStateV::InfectedTV2)});
        auto CV2i  = model.populations.get_flat_index({i, size_t(InfectionStateV::CarrierV2)});
        auto CDV2i = model.populations.get_flat_index({i, size_t(InfectionStateV::CarrierTV2)});*/

        //put detected commuters in their own compartment so they don't contribute to infections in their home node
        sim.get_result().get_last_value()[Ii] -= migrated[Ii] * (1 - nondetection);
        sim.get_result().get_last_value()[IDi] += migrated[Ii] * (1 - nondetection);
        sim.get_result().get_last_value()[Ci] -= migrated[Ci] * (1 - nondetection);
        sim.get_result().get_last_value()[CDi] += migrated[Ci] * (1 - nondetection);

        sim.get_result().get_last_value()[IV1i] -= migrated[IV1i] * (1 - nondetection);
        sim.get_result().get_last_value()[IDV1i] += migrated[IV1i] * (1 - nondetection);
        sim.get_result().get_last_value()[CV1i] -= migrated[CV1i] * (1 - nondetection);
        sim.get_result().get_last_value()[CDV1i] += migrated[CV1i] * (1 - nondetection);

        /*sim.get_result().get_last_value()[IV2i] -= migrated[IV2i] * (1 - nondetection);
        sim.get_result().get_last_value()[IDV2i] += migrated[IV2i] * (1 - nondetection);
        sim.get_result().get_last_value()[CV2i] -= migrated[CV2i] * (1 - nondetection);
        sim.get_result().get_last_value()[CDV2i] += migrated[CV2i] * (1 - nondetection);*/

        //reduce the number of commuters
        migrated[Ii] *= nondetection;
        migrated[Ci] *= nondetection;

        migrated[IV1i] *= nondetection;
        migrated[CV1i] *= nondetection;

        /*migrated[IV2i] *= nondetection;
        migrated[CV2i] *= nondetection;*/
    }
}

} // namespace epi

#endif
