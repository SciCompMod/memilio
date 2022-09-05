/* 
* Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
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
            auto dummy_R3 = 0.5 / (params.get<IncubationTime>()[i] - params.get<SerialInterval>()[i]);
            test_and_trace_required +=
                (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 *
                (this->populations.get_from(pop, {i, InfectionState::CarrierNaive}) +
                 this->populations.get_from(pop, {i, InfectionState::CarrierPartialImmunity}) +
                 this->populations.get_from(pop, {i, InfectionState::CarrierImprovedImmunity}) +
                 this->populations.get_from(pop, {i, InfectionState::CarrierNaiveConfirmed}) +
                 this->populations.get_from(pop, {i, InfectionState::CarrierPartialImmunityConfirmed}) +
                 this->populations.get_from(pop, {i, InfectionState::CarrierImprovedImmunityConfirmed}));
            icu_occupancy += this->populations.get_from(pop, {i, InfectionState::ICUNaive}) +
                             this->populations.get_from(pop, {i, InfectionState::ICUPartialImmunity}) +
                             this->populations.get_from(pop, {i, InfectionState::ICUImprovedImmunity});
        }

        for (auto i = AgeGroup(0); i < n_agegroups; i++) {

            size_t Si = this->populations.get_flat_index({i, InfectionState::SusceptibleNaive});
            size_t Ei = this->populations.get_flat_index({i, InfectionState::ExposedNaive});
            size_t Ci = this->populations.get_flat_index({i, InfectionState::CarrierNaive});
            size_t Ii = this->populations.get_flat_index({i, InfectionState::InfectedNaive});
            size_t Hi = this->populations.get_flat_index({i, InfectionState::HospitalizedNaive});
            size_t Ui = this->populations.get_flat_index({i, InfectionState::ICUNaive});

            size_t CDi = this->populations.get_flat_index({i, InfectionState::CarrierNaiveConfirmed});
            size_t IDi = this->populations.get_flat_index({i, InfectionState::InfectedNaiveConfirmed});

            size_t SVi = this->populations.get_flat_index({i, InfectionState::SusceptiblePartialImmunity});
            size_t EVi = this->populations.get_flat_index({i, InfectionState::ExposedPartialImmunity});
            size_t CVi = this->populations.get_flat_index({i, InfectionState::CarrierPartialImmunity});
            size_t IVi = this->populations.get_flat_index({i, InfectionState::InfectedPartialImmunity});
            size_t HVi = this->populations.get_flat_index({i, InfectionState::HospitalizedPartialImmunity});
            size_t UVi = this->populations.get_flat_index({i, InfectionState::ICUPartialImmunity});

            size_t CDVi = this->populations.get_flat_index({i, InfectionState::CarrierPartialImmunityConfirmed});
            size_t IDVi = this->populations.get_flat_index({i, InfectionState::InfectedPartialImmunityConfirmed});

            size_t EV2i = this->populations.get_flat_index({i, InfectionState::ExposedImprovedImmunity});
            size_t CV2i = this->populations.get_flat_index({i, InfectionState::CarrierImprovedImmunity});
            size_t IV2i = this->populations.get_flat_index({i, InfectionState::InfectedImprovedImmunity});
            size_t HV2i = this->populations.get_flat_index({i, InfectionState::HospitalizedImprovedImmunity});
            size_t UV2i = this->populations.get_flat_index({i, InfectionState::ICUImprovedImmunity});

            size_t CDV2i = this->populations.get_flat_index({i, InfectionState::CarrierImprovedImmunityConfirmed});
            size_t IDV2i = this->populations.get_flat_index({i, InfectionState::InfectedImprovedImmunityConfirmed});

            size_t TImm1 = this->populations.get_flat_index({i, InfectionState::TemporaryImmunity1});
            size_t TImm2 = this->populations.get_flat_index({i, InfectionState::TemporaryImmunity2});

            size_t Ri = this->populations.get_flat_index({i, InfectionState::Recovered});
            size_t Di = this->populations.get_flat_index({i, InfectionState::Dead});

            size_t ITi = this->populations.get_flat_index({i, InfectionState::TotalInfections});

            dydt[Si] = 0;
            dydt[Ei] = 0;

            dydt[SVi] = 0;
            dydt[EVi] = 0;

            dydt[TImm1] = 0;
            dydt[TImm2] = 0;

            dydt[Ri]   = 0;
            dydt[EV2i] = 0;

            double dummy_R2 = 1.0 / (2 * params.get<SerialInterval>()[i] - params.get<IncubationTime>()[i]);
            double dummy_R3 = 0.5 / (params.get<IncubationTime>()[i] - params.get<SerialInterval>()[i]);

            double exp_fac_part_immune    = params.get<ExposedFactorPartialImmunity>()[i];
            double exp_fac_impr_immune    = params.get<ExposedFactorImprovedImmunity>()[i];
            double inf_fac_part_immune    = params.get<InfectedFactorPartialImmunity>()[i];
            double inf_fac_impr_immune    = params.get<InfectedFactorImprovedImmunity>()[i];
            double hosp_fac_part_immune   = params.get<HospitalizedFactorPartialImmunity>()[i];
            double icu_fac_part_immune    = params.get<HospitalizedFactorPartialImmunity>()[i];
            double death_fac_part_immune  = params.get<HospitalizedFactorPartialImmunity>()[i];
            double hosp_fac_impr_immune   = params.get<HospitalizedFactorImprovedImmunity>()[i];
            double icu_fac_impr_immune    = params.get<HospitalizedFactorImprovedImmunity>()[i];
            double death_fac_impr_immune  = params.get<HospitalizedFactorImprovedImmunity>()[i];
            double inf_time_factor_immune = params.get<InfectiousTimeFactorImmune>()[i];

            //symptomatic are less well quarantined when testing and tracing is overwhelmed so they infect more people
            auto risk_from_symptomatic = smoother_cosine(
                test_and_trace_required, params.get<TestAndTraceCapacity>(), params.get<TestAndTraceCapacity>() * 15,
                params.get<RiskOfInfectionFromSympomatic>()[i], params.get<MaxRiskOfInfectionFromSympomatic>()[i]);

            auto risk_from_carrier = smoother_cosine(test_and_trace_required, params.get<TestAndTraceCapacity>(),
                                                     params.get<TestAndTraceCapacity>() * 2,
                                                     params.get<RelativeCarrierInfectability>()[i], 1.0);

            for (auto j = AgeGroup(0); j < n_agegroups; j++) {
                size_t Sj = this->populations.get_flat_index({j, InfectionState::SusceptibleNaive});
                size_t Ej = this->populations.get_flat_index({j, InfectionState::ExposedNaive});
                size_t Cj = this->populations.get_flat_index({j, InfectionState::CarrierNaive});
                size_t Ij = this->populations.get_flat_index({j, InfectionState::InfectedNaive});
                size_t Hj = this->populations.get_flat_index({j, InfectionState::HospitalizedNaive});
                size_t Uj = this->populations.get_flat_index({j, InfectionState::ICUNaive});
                size_t Rj = this->populations.get_flat_index({j, InfectionState::Recovered});

                size_t SVj = this->populations.get_flat_index({j, InfectionState::SusceptiblePartialImmunity});
                size_t EVj = this->populations.get_flat_index({j, InfectionState::ExposedPartialImmunity});
                size_t CVj = this->populations.get_flat_index({j, InfectionState::CarrierPartialImmunity});
                size_t IVj = this->populations.get_flat_index({j, InfectionState::InfectedPartialImmunity});
                size_t HVj = this->populations.get_flat_index({j, InfectionState::HospitalizedPartialImmunity});
                size_t UVj = this->populations.get_flat_index({j, InfectionState::ICUPartialImmunity});

                size_t CDj = this->populations.get_flat_index({j, InfectionState::CarrierNaiveConfirmed});
                size_t IDj = this->populations.get_flat_index({j, InfectionState::InfectedNaiveConfirmed});

                size_t CDVj = this->populations.get_flat_index({j, InfectionState::CarrierPartialImmunityConfirmed});
                size_t IDVj = this->populations.get_flat_index({j, InfectionState::InfectedPartialImmunityConfirmed});

                size_t EV2j = this->populations.get_flat_index({j, InfectionState::ExposedImprovedImmunity});
                size_t CV2j = this->populations.get_flat_index({j, InfectionState::CarrierImprovedImmunity});
                size_t IV2j = this->populations.get_flat_index({j, InfectionState::InfectedImprovedImmunity});
                size_t HV2j = this->populations.get_flat_index({j, InfectionState::HospitalizedImprovedImmunity});
                size_t UV2j = this->populations.get_flat_index({j, InfectionState::ICUImprovedImmunity});

                size_t CDV2j = this->populations.get_flat_index({j, InfectionState::CarrierImprovedImmunityConfirmed});
                size_t IDV2j = this->populations.get_flat_index({j, InfectionState::InfectedImprovedImmunityConfirmed});

                // effective contact rate by contact rate between groups i and j and damping j
                double season_val =
                    (1 + params.get<Seasonality>() *
                             sin(3.141592653589793 * (std::fmod((params.get<StartDay>() + t), 365.0) / 182.5 + 0.5)));
                double cont_freq_eff =
                    season_val * contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                 static_cast<Eigen::Index>((size_t)j));

                double Nj = pop[Sj] + pop[Ej] + pop[Cj] + pop[Ij] + pop[Hj] + pop[Uj] + pop[Rj] + pop[CDj] + pop[IDj] +
                            pop[SVj] + pop[EVj] + pop[CVj] + pop[IVj] + pop[HVj] + pop[UVj] + pop[CDVj] + pop[IDVj] +
                            pop[EV2j] + pop[CV2j] + pop[IV2j] + pop[HV2j] + pop[UV2j] + pop[CDV2j] +
                            pop[IDV2j]; // without died people

                double divNj = 1.0 / Nj; // precompute 1.0/Nj

                double ext_inf_force_dummy = cont_freq_eff * divNj *
                                             params.template get<InfectionProbabilityFromContact>()[(AgeGroup)i] *
                                             (risk_from_carrier * (pop[Cj] + pop[CVj] + pop[CV2j]) +
                                              risk_from_symptomatic * (pop[Ij] + pop[IVj] + pop[IV2j]));

                double dummy_partial_imm = y[Si] * params.get<RateOfDailyPartialVaccinations>()[i] /
                                           params.get<DaysUntilEffectivePartialImmunity>()[i];
                double dummy_improved_imm = y[SVi] * params.get<RateOfDailyImprovedVaccinations>()[i] /
                                            params.get<DaysUntilEffectivePartialImmunity>()[i];

                double dummy_S = y[Si] * ext_inf_force_dummy;

                double dummy_SV = y[SVi] * exp_fac_part_immune * ext_inf_force_dummy;

                double dummy_R = y[Ri] * exp_fac_impr_immune * ext_inf_force_dummy;

                dydt[Si] -= dummy_S + dummy_partial_imm;
                dydt[Ei] += dummy_S;

                dydt[SVi] -= dummy_SV + dummy_improved_imm;
                dydt[EVi] += dummy_SV;

                dydt[TImm1] += dummy_partial_imm;
                dydt[TImm2] += dummy_improved_imm;

                dydt[Ri] -= dummy_R;
                dydt[EV2i] += dummy_R;

                // waning immunity
                dydt[Si] += 1 / params.get<WainingPartialImmunity>()[i] * y[SVi];
                dydt[SVi] -= 1 / params.get<WainingPartialImmunity>()[i] * y[SVi];
                dydt[SVi] += 1 / params.get<WainingImprovedImmunity>()[i] * y[Ri];
                dydt[Ri] -= 1 / params.get<WainingImprovedImmunity>()[i] * y[Ri];
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
            dydt[Ui] = -((1 - params.get<DeathsPerICU>()[i]) / params.get<ICUToHomeTime>()[i] +
                         params.get<DeathsPerICU>()[i] / params.get<ICUToDeathTime>()[i]) *
                       y[Ui];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[Ui] += prob_hosp2icu / params.get<HospitalizedToICUTime>()[i] * y[Hi];

            /**** path of (partially, i.e., one dose) vaccinated ***/

            dydt[EVi] -= dummy_R2 * y[EVi]; // only exchange of E and C done here
            dydt[CVi] =
                dummy_R2 * y[EVi] - ((inf_fac_part_immune / exp_fac_part_immune) *
                                         (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 +
                                     (1 - (inf_fac_part_immune / exp_fac_part_immune) *
                                              (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                                         (params.get<InfectiousTimeAsymptomatic>()[i] * inf_time_factor_immune)) *
                                        y[CVi];
            dydt[CDVi] = -((inf_fac_part_immune / exp_fac_part_immune) *
                               (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 +
                           (1 - (inf_fac_part_immune / exp_fac_part_immune) *
                                    (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                               (params.get<InfectiousTimeAsymptomatic>()[i] * inf_time_factor_immune)) *
                         y[CDVi];
            dydt[IVi] =
                (inf_fac_part_immune / exp_fac_part_immune) * (1 - params.get<AsymptoticCasesPerInfectious>()[i]) *
                    dummy_R3 * y[CVi] -
                ((1 - hosp_fac_part_immune / inf_fac_part_immune * params.get<HospitalizedCasesPerInfectious>()[i]) /
                     (params.get<InfectiousTimeMild>()[i] * inf_time_factor_immune) +
                 hosp_fac_part_immune / inf_fac_part_immune * params.get<HospitalizedCasesPerInfectious>()[i] /
                     params.get<HomeToHospitalizedTime>()[i]) *
                    y[IVi];
            dydt[IDVi] =
                (inf_fac_part_immune / exp_fac_part_immune) * (1 - params.get<AsymptoticCasesPerInfectious>()[i]) *
                    dummy_R3 * y[CDVi] -
                ((1 - hosp_fac_part_immune / inf_fac_part_immune * params.get<HospitalizedCasesPerInfectious>()[i]) /
                     (params.get<InfectiousTimeMild>()[i] * inf_time_factor_immune) +
                 hosp_fac_part_immune / inf_fac_part_immune * params.get<HospitalizedCasesPerInfectious>()[i] /
                     params.get<HomeToHospitalizedTime>()[i]) *
                    y[IDVi];
            dydt[HVi] = hosp_fac_part_immune / inf_fac_part_immune * params.get<HospitalizedCasesPerInfectious>()[i] /
                            params.get<HomeToHospitalizedTime>()[i] * y[IVi] +
                        hosp_fac_part_immune / inf_fac_part_immune * params.get<HospitalizedCasesPerInfectious>()[i] /
                            params.get<HomeToHospitalizedTime>()[i] * y[IDVi] -
                        ((1 - icu_fac_part_immune / hosp_fac_part_immune * params.get<ICUCasesPerHospitalized>()[i]) /
                             params.get<HospitalizedToHomeTime>()[i] +
                         icu_fac_part_immune / hosp_fac_part_immune * params.get<ICUCasesPerHospitalized>()[i] /
                             params.get<HospitalizedToICUTime>()[i]) *
                            y[HVi];
            dydt[UVi] = -((1 - death_fac_part_immune / icu_fac_part_immune * params.get<DeathsPerICU>()[i]) /
                              params.get<ICUToHomeTime>()[i] +
                          death_fac_part_immune / icu_fac_part_immune * params.get<DeathsPerICU>()[i] /
                              params.get<ICUToDeathTime>()[i]) *
                        y[UVi];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[UVi] += icu_fac_part_immune / hosp_fac_part_immune * prob_hosp2icu /
                         params.get<HospitalizedToICUTime>()[i] * y[HVi];

            /**** path of twice vaccinated, here called immune although reinfection is possible now ***/

            dydt[EV2i] -= dummy_R2 * y[EV2i]; // only exchange of E and C done here

            dydt[CV2i] =
                dummy_R2 * y[EV2i] - ((inf_fac_impr_immune / exp_fac_impr_immune) *
                                          (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 +
                                      (1 - (inf_fac_impr_immune / exp_fac_impr_immune) *
                                               (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                                          (params.get<InfectiousTimeAsymptomatic>()[i] * inf_time_factor_immune)) *
                                         y[CV2i];
            dydt[CDV2i] = -((inf_fac_impr_immune / exp_fac_impr_immune) *
                                (1 - params.get<AsymptoticCasesPerInfectious>()[i]) * dummy_R3 +
                            (1 - (inf_fac_impr_immune / exp_fac_impr_immune) *
                                     (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                                (params.get<InfectiousTimeAsymptomatic>()[i] * inf_time_factor_immune)) *
                          y[CDV2i];

            dydt[IV2i] =
                (inf_fac_impr_immune / exp_fac_impr_immune) * (1 - params.get<AsymptoticCasesPerInfectious>()[i]) *
                    dummy_R3 * y[CV2i] -
                ((1 - hosp_fac_impr_immune / inf_fac_impr_immune * params.get<HospitalizedCasesPerInfectious>()[i]) /
                     (params.get<InfectiousTimeMild>()[i] * inf_time_factor_immune) +
                 hosp_fac_impr_immune / inf_fac_impr_immune * params.get<HospitalizedCasesPerInfectious>()[i] /
                     params.get<HomeToHospitalizedTime>()[i]) *
                    y[IV2i];
            dydt[IDV2i] =
                (inf_fac_impr_immune / exp_fac_impr_immune) * (1 - params.get<AsymptoticCasesPerInfectious>()[i]) *
                    dummy_R3 * y[CDV2i] -
                ((1 - hosp_fac_impr_immune / inf_fac_impr_immune * params.get<HospitalizedCasesPerInfectious>()[i]) /
                     (params.get<InfectiousTimeMild>()[i] * inf_time_factor_immune) +
                 hosp_fac_impr_immune / inf_fac_impr_immune * params.get<HospitalizedCasesPerInfectious>()[i] /
                     params.get<HomeToHospitalizedTime>()[i]) *
                    y[IDV2i];
            dydt[HV2i] = hosp_fac_impr_immune / inf_fac_impr_immune * params.get<HospitalizedCasesPerInfectious>()[i] /
                             params.get<HomeToHospitalizedTime>()[i] * y[IV2i] +
                         hosp_fac_impr_immune / inf_fac_impr_immune * params.get<HospitalizedCasesPerInfectious>()[i] /
                             params.get<HomeToHospitalizedTime>()[i] * y[IDV2i] -
                         ((1 - icu_fac_impr_immune / hosp_fac_impr_immune * params.get<ICUCasesPerHospitalized>()[i]) /
                              params.get<HospitalizedToHomeTime>()[i] +
                          icu_fac_impr_immune / hosp_fac_impr_immune * params.get<ICUCasesPerHospitalized>()[i] /
                              params.get<HospitalizedToICUTime>()[i]) *
                             y[HV2i];
            dydt[UV2i] = -((1 - death_fac_impr_immune / icu_fac_impr_immune * params.get<DeathsPerICU>()[i]) /
                               params.get<ICUToHomeTime>()[i] +
                           death_fac_impr_immune / icu_fac_impr_immune * params.get<DeathsPerICU>()[i] /
                               params.get<ICUToDeathTime>()[i]) *
                         y[UV2i];
            // add flow from hosp to icu according to potentially adjusted probability due to ICU limits
            dydt[UV2i] += icu_fac_impr_immune / hosp_fac_impr_immune * prob_hosp2icu /
                          params.get<HospitalizedToICUTime>()[i] * y[HV2i];

            // compute auxiliary compartment of all past infections
            dydt[ITi] =
                ((1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] +
                 params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i]) *
                    y[Ii] +
                ((1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] +
                 params.get<HospitalizedCasesPerInfectious>()[i] / params.get<HomeToHospitalizedTime>()[i]) *
                    y[IDi] +
                ((1 - hosp_fac_part_immune / inf_fac_part_immune * params.get<HospitalizedCasesPerInfectious>()[i]) /
                     (params.get<InfectiousTimeMild>()[i] * inf_time_factor_immune) +
                 hosp_fac_part_immune / inf_fac_part_immune * params.get<HospitalizedCasesPerInfectious>()[i] /
                     params.get<HomeToHospitalizedTime>()[i]) *
                    y[IDVi] +
                ((1 - hosp_fac_part_immune / inf_fac_part_immune * params.get<HospitalizedCasesPerInfectious>()[i]) /
                     (params.get<InfectiousTimeMild>()[i] * inf_time_factor_immune) +
                 hosp_fac_part_immune / inf_fac_part_immune * params.get<HospitalizedCasesPerInfectious>()[i] /
                     params.get<HomeToHospitalizedTime>()[i]) *
                    y[IVi] +
                ((1 - hosp_fac_impr_immune / inf_fac_impr_immune * params.get<HospitalizedCasesPerInfectious>()[i]) /
                     (params.get<InfectiousTimeMild>()[i] * inf_time_factor_immune) +
                 hosp_fac_impr_immune / inf_fac_impr_immune * params.get<HospitalizedCasesPerInfectious>()[i] /
                     params.get<HomeToHospitalizedTime>()[i]) *
                    y[IDV2i] +
                ((1 - hosp_fac_impr_immune / inf_fac_impr_immune * params.get<HospitalizedCasesPerInfectious>()[i]) /
                     (params.get<InfectiousTimeMild>()[i] * inf_time_factor_immune) +
                 hosp_fac_impr_immune / inf_fac_impr_immune * params.get<HospitalizedCasesPerInfectious>()[i] /
                     params.get<HomeToHospitalizedTime>()[i]) *
                    y[IV2i];

            // recovered and deaths from all paths
            dydt[TImm1] +=
                params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i] * y[Ci] +
                params.get<AsymptoticCasesPerInfectious>()[i] / params.get<InfectiousTimeAsymptomatic>()[i] * y[CDi] +
                (1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] * y[Ii] +
                (1 - params.get<HospitalizedCasesPerInfectious>()[i]) / params.get<InfectiousTimeMild>()[i] * y[IDi] +
                (1 - params.get<ICUCasesPerHospitalized>()[i]) / params.get<HospitalizedToHomeTime>()[i] * y[Hi] +
                (1 - params.get<DeathsPerICU>()[i]) / params.get<ICUToHomeTime>()[i] * y[Ui];

            dydt[TImm2] +=
                (1 -
                 (inf_fac_part_immune / exp_fac_part_immune) * (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                    (params.get<InfectiousTimeAsymptomatic>()[i] * inf_time_factor_immune) * y[CVi] +
                (1 -
                 (inf_fac_part_immune / exp_fac_part_immune) * (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                    (params.get<InfectiousTimeAsymptomatic>()[i] * inf_time_factor_immune) * y[CDVi] +
                (1 - (hosp_fac_part_immune / inf_fac_part_immune) * params.get<HospitalizedCasesPerInfectious>()[i]) /
                    (params.get<InfectiousTimeMild>()[i] * inf_time_factor_immune) * y[IVi] +
                (1 - (hosp_fac_part_immune / inf_fac_part_immune) * params.get<HospitalizedCasesPerInfectious>()[i]) /
                    (params.get<InfectiousTimeMild>()[i] * inf_time_factor_immune) * y[IDVi] +
                (1 - (icu_fac_part_immune / hosp_fac_part_immune) * params.get<ICUCasesPerHospitalized>()[i]) /
                    params.get<HospitalizedToHomeTime>()[i] * y[HVi] +
                (1 - (death_fac_part_immune / icu_fac_part_immune) * params.get<DeathsPerICU>()[i]) /
                    params.get<ICUToHomeTime>()[i] * y[UVi];

            dydt[TImm2] +=
                (1 -
                 (inf_fac_impr_immune / exp_fac_impr_immune) * (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                    (params.get<InfectiousTimeAsymptomatic>()[i] * inf_time_factor_immune) * y[CV2i] +
                (1 -
                 (inf_fac_impr_immune / exp_fac_impr_immune) * (1 - params.get<AsymptoticCasesPerInfectious>()[i])) /
                    (params.get<InfectiousTimeAsymptomatic>()[i] * inf_time_factor_immune) * y[CDV2i] +
                (1 - (hosp_fac_impr_immune / inf_fac_impr_immune) * params.get<HospitalizedCasesPerInfectious>()[i]) /
                    (params.get<InfectiousTimeMild>()[i] * inf_time_factor_immune) * y[IV2i] +
                (1 - (hosp_fac_impr_immune / inf_fac_impr_immune) * params.get<HospitalizedCasesPerInfectious>()[i]) /
                    (params.get<InfectiousTimeMild>()[i] * inf_time_factor_immune) * y[IDV2i] +
                (1 - (icu_fac_impr_immune / hosp_fac_impr_immune) * params.get<ICUCasesPerHospitalized>()[i]) /
                    params.get<HospitalizedToHomeTime>()[i] * y[HV2i] +
                (1 - (death_fac_impr_immune / icu_fac_impr_immune) * params.get<DeathsPerICU>()[i]) /
                    params.get<ICUToHomeTime>()[i] * y[UV2i];

            double dummy_TImm1 = 1 / params.get<ImmunityInterval1>()[i];
            double dummy_TImm2 = 1 / params.get<ImmunityInterval2>()[i];
            dydt[SVi]          = dummy_TImm1 * y[TImm1];
            dydt[Ri]           = dummy_TImm2 * y[TImm2];

            dydt[Di] = params.get<DeathsPerICU>()[i] / params.get<ICUToDeathTime>()[i] * y[Ui] +
                       death_fac_part_immune / icu_fac_part_immune * params.get<DeathsPerICU>()[i] /
                           params.get<ICUToDeathTime>()[i] * y[UVi] +
                       death_fac_impr_immune / icu_fac_impr_immune * params.get<DeathsPerICU>()[i] /
                           params.get<ICUToDeathTime>()[i] * y[UV2i];
            ;
            // add potential, additional deaths due to ICU overflow
            dydt[Di] += prob_hosp2dead / params.get<HospitalizedToICUTime>()[i] * y[Hi];

            dydt[Di] += (death_fac_part_immune / hosp_fac_part_immune) * prob_hosp2dead /
                        params.get<HospitalizedToICUTime>()[i] * y[HVi];

            dydt[Di] += (death_fac_impr_immune / hosp_fac_impr_immune) * prob_hosp2dead /
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

    void apply_variant(double t)
    {
        // TODO: Wachtstumsrate anpassen/ evtl als Parameter
        auto start_day   = this->get_model().parameters.template get<StartDay>();
        auto variant_day = get_day_in_year(Date(2022, 11, 1));
        if (start_day + t > variant_day) {
            auto variant_growth = (start_day - variant_day) * 0.1666667;
            // 2 equal to the share of the delta variant on June 6
            double share_new_variant = std::min(1.0, pow(2, t * 0.1666667 + variant_growth) * 0.01);
            size_t num_groups        = (size_t)this->get_model().parameters.get_num_groups();
            for (size_t i = 0; i < num_groups; ++i) {
                double new_transmission =
                    (1 - share_new_variant) *
                        this->get_model().parameters.template get<BaseInfectiousness>()[(AgeGroup)i] +
                    share_new_variant *
                        this->get_model().parameters.template get<BaseInfectiousnessNewVariant>()[(AgeGroup)i];
                this->get_model().parameters.template get<InfectionProbabilityFromContact>()[(AgeGroup)i] =
                    new_transmission;

                double new_severtiy =
                    (1 - share_new_variant) * this->get_model().parameters.template get<BaseSeverity>()[(AgeGroup)i] +
                    share_new_variant *
                        this->get_model().parameters.template get<BaseSeverityNewVariant>()[(AgeGroup)i];
                this->get_model().parameters.template get<HospitalizedCasesPerInfectious>()[(AgeGroup)i] = new_severtiy;
            }
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
        auto R     = (size_t)InfectionState::Recovered;

        for (size_t i = 0; i < num_groups; ++i) {

            double first_vacc;
            double full_vacc;
            if (t_idx == SimulationDay(0)) {
                first_vacc = params.template get<DailyPartialVaccination>()[{(AgeGroup)i, t_idx}];
                full_vacc  = params.template get<DailyFullVaccination>()[{(AgeGroup)i, t_idx}];
            }
            else {
                first_vacc = params.template get<DailyPartialVaccination>()[{(AgeGroup)i, t_idx}] -
                             params.template get<DailyPartialVaccination>()[{(AgeGroup)i, t_idx - SimulationDay(1)}];
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
                this->apply_variant(t);
            }
            Base::advance(t + dt_eff);
            if (t + 0.5 + dt_eff - std::floor(t + 0.5) >= 1) {
                // this->apply_vaccination(t + 0.5 + dt_eff);
                this->apply_variant(t);
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
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedNaive});
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedNaiveConfirmed});
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedPartialImmunity});
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedImprovedImmunity});
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedPartialImmunityConfirmed});
        sum_inf += sim.get_model().populations.get_from(y, {i, InfectionState::InfectedImprovedImmunityConfirmed});
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
    auto& p_asymp   = params.template get<AsymptoticCasesPerInfectious>().array().template cast<double>();
    auto& p_inf     = params.template get<RiskOfInfectionFromSympomatic>().array().template cast<double>();
    auto& p_inf_max = params.template get<MaxRiskOfInfectionFromSympomatic>().array().template cast<double>();
    //slice of carriers
    auto y_car = slice(y, {Eigen::Index(InfectionState::CarrierNaive), Eigen::Index(size_t(params.get_num_groups())),
                           Eigen::Index(InfectionState::Count)}) +
                 slice(y, {Eigen::Index(InfectionState::CarrierPartialImmunity),
                           Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)}) +
                 slice(y, {Eigen::Index(InfectionState::CarrierImprovedImmunity),
                           Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)});

    //compute isolation, same as infection risk from main model
    auto R3                      = 0.5 / (t_inc - t_ser);
    auto test_and_trace_required = ((1 - p_asymp) * R3 * y_car.array()).sum();
    auto risk_from_symptomatic =
        smoother_cosine(test_and_trace_required, double(params.template get<TestAndTraceCapacity>()),
                        params.template get<TestAndTraceCapacity>() * 5, p_inf.matrix(), p_inf_max.matrix());

    //set factor for infected
    auto factors = Eigen::VectorXd::Ones(y.rows()).eval();
    slice(factors, {Eigen::Index(InfectionState::InfectedNaive), Eigen::Index(size_t(params.get_num_groups())),
                    Eigen::Index(InfectionState::Count)})
        .array() = risk_from_symptomatic;
    slice(factors, {Eigen::Index(InfectionState::InfectedPartialImmunity),
                    Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)})
        .array() = risk_from_symptomatic;
    slice(factors, {Eigen::Index(InfectionState::InfectedImprovedImmunity),
                    Eigen::Index(size_t(params.get_num_groups())), Eigen::Index(InfectionState::Count)})
        .array() = risk_from_symptomatic;
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
        auto Ii  = model.populations.get_flat_index({i, InfectionState::InfectedNaive});
        auto IDi = model.populations.get_flat_index({i, InfectionState::InfectedNaiveConfirmed});
        auto Ci  = model.populations.get_flat_index({i, InfectionState::CarrierNaive});
        auto CDi = model.populations.get_flat_index({i, InfectionState::CarrierNaiveConfirmed});

        auto IV1i  = model.populations.get_flat_index({i, InfectionState::InfectedPartialImmunity});
        auto IDV1i = model.populations.get_flat_index({i, InfectionState::InfectedPartialImmunityConfirmed});
        auto CV1i  = model.populations.get_flat_index({i, InfectionState::CarrierPartialImmunity});
        auto CDV1i = model.populations.get_flat_index({i, InfectionState::CarrierPartialImmunityConfirmed});

        auto IV2i  = model.populations.get_flat_index({i, InfectionState::InfectedImprovedImmunity});
        auto IDV2i = model.populations.get_flat_index({i, InfectionState::InfectedImprovedImmunityConfirmed});
        auto CV2i  = model.populations.get_flat_index({i, InfectionState::CarrierImprovedImmunity});
        auto CDV2i = model.populations.get_flat_index({i, InfectionState::CarrierImprovedImmunityConfirmed});

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

} // namespace osecirvvs
} // namespace mio

#endif //ODESECIRVVS_MODEL_H
