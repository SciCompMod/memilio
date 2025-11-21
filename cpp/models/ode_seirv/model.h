/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#ifndef SEIRV_MODEL_H
#define SEIRV_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/utils/time_series.h"
#include "ode_seirv/infection_state.h"
#include "ode_seirv/parameters.h"

GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wshadow")
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
GCC_CLANG_DIAGNOSTIC(pop)

#include <cmath>
#include <numbers>

namespace mio
{
namespace oseirv
{

// clang-format off
using Flows = TypeList<Flow<InfectionState::Susceptible,            InfectionState::Exposed>,
                       Flow<InfectionState::SusceptibleVaccinated,  InfectionState::ExposedVaccinated>,
                       Flow<InfectionState::Exposed,                InfectionState::Infected>,
                       Flow<InfectionState::ExposedVaccinated,      InfectionState::InfectedVaccinated>,
                       Flow<InfectionState::Infected,               InfectionState::Recovered>,
                       Flow<InfectionState::InfectedVaccinated,     InfectionState::RecoveredVaccinated>>;
// clang-format on

template <typename FP>
class Model
    : public FlowModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>
{
    using Base = FlowModel<FP, InfectionState, mio::Populations<FP, AgeGroup, InfectionState>, Parameters<FP>, Flows>;

public:
    using ParameterSet    = typename Base::ParameterSet;
    using PopulationsType = typename Base::Populations;

    Model(const PopulationsType& pop, const ParameterSet& params)
        : Base(pop, params)
    {
    }
    Model(int num_agegroups)
        : Base(PopulationsType({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        const auto& params      = this->parameters;
        const size_t num_groups = (size_t)params.get_num_groups();

        // Get effective contact matrix effective_contacts  = contacts_healthy + m*contacts_sick
        Eigen::MatrixX<FP> contacts_healthy =
            params.template get<ContactPatternsHealthy<FP>>().get_cont_freq_mat().get_matrix_at(SimulationTime<FP>(t));
        Eigen::MatrixX<FP> contacts_sick =
            params.template get<ContactPatternsSick<FP>>().get_cont_freq_mat().get_matrix_at(SimulationTime<FP>(t));
        const FP m                            = params.template get<SickMixing<FP>>();
        Eigen::MatrixX<FP> effective_contacts = contacts_healthy + m * contacts_sick;

        // Normalization ν(m): prepare susceptibility scaling and next-generation matrix (NG)
        Eigen::VectorX<FP> susceptibility(num_groups);
        for (size_t i = 0; i < num_groups; ++i)
            susceptibility[(Eigen::Index)i] = params.template get<SusceptibilityByAge<FP>>()[AgeGroup(i)] *
                                              params.template get<SusceptibleFraction<FP>>();
        // Sigma is a diagonal matrix applying susceptibility on the receiving (row) side of contacts.
        Eigen::MatrixX<FP> Sigma = susceptibility.asDiagonal();

        // NG = Sigma * effective_contacts is the (unnormalized) next-generation matrix. Its spectral radius will be used
        // below to scale contacts effective_contacts so that structure (mixing pattern) and transmissibility magnitude (BaselineTransmissibility)
        // are cleanly separated.
        Eigen::MatrixX<FP> NG = Sigma * effective_contacts;

        Eigen::ComplexEigenSolver<Eigen::MatrixX<FP>> ces;
        ces.compute(NG);
        FP scale_transmissibility = ces.eigenvalues().cwiseAbs().maxCoeff();
        if (scale_transmissibility > FP(0.)) {
            effective_contacts /= scale_transmissibility;
        }

        const FP transmissibility_baseline =
            params.template get<BaselineTransmissibility<FP>>(); // baseline transmissibility scaling (R-like factor)
        const FP infection_rate   = FP(1) / params.template get<TimeExposed<FP>>(); // 1 / mean latent duration
        const FP recovery_rate    = FP(1) / params.template get<TimeInfected<FP>>(); // 1 / mean infectious duration
        const FP season_amplitude = params.template get<SeasonalityAmplitude<FP>>(); // amplitude of seasonal modulation
        const FP season_shift_per_subtype =
            params.template get<SeasonalityShiftPerSubtype<FP>>(); // base phase shift for seasonality
        const FP season_shift_per_season =
            params.template get<SeasonalityShiftPerSeason<FP>>(); // additional phase shift (e.g. season-specific)
        const FP outside_foi    = params.template get<OutsideFoI<FP>>(); // constant external force of infection
        const FP clustering_exp = params.template get<ClusteringExponent<
            FP>>(); // clustering exponent; clustering_exp<1 dampens, clustering_exp>1 amplifies prevalence effect

        const FP season =
            std::exp(season_amplitude * std::sin(FP(2) * std::numbers::pi_v<FP> *
                                                 (t / FP(52.0) - season_shift_per_subtype + season_shift_per_season)));

        for (auto i : make_index_range((mio::AgeGroup)num_groups)) {
            // Flat indices
            const size_t Si  = this->populations.get_flat_index({i, InfectionState::Susceptible});
            const size_t SVi = this->populations.get_flat_index({i, InfectionState::SusceptibleVaccinated});
            const size_t Ei  = this->populations.get_flat_index({i, InfectionState::Exposed});
            const size_t EVi = this->populations.get_flat_index({i, InfectionState::ExposedVaccinated});
            const size_t Ii  = this->populations.get_flat_index({i, InfectionState::Infected});
            const size_t IVi = this->populations.get_flat_index({i, InfectionState::InfectedVaccinated});

            FP sum = FP(0);
            for (auto j : make_index_range((mio::AgeGroup)num_groups)) {
                const size_t Sj  = this->populations.get_flat_index({j, InfectionState::Susceptible});
                const size_t SVj = this->populations.get_flat_index({j, InfectionState::SusceptibleVaccinated});
                const size_t Ej  = this->populations.get_flat_index({j, InfectionState::Exposed});
                const size_t EVj = this->populations.get_flat_index({j, InfectionState::ExposedVaccinated});
                const size_t Ij  = this->populations.get_flat_index({j, InfectionState::Infected});
                const size_t IVj = this->populations.get_flat_index({j, InfectionState::InfectedVaccinated});
                const size_t Rj  = this->populations.get_flat_index({j, InfectionState::Recovered});
                const size_t RVj = this->populations.get_flat_index({j, InfectionState::RecoveredVaccinated});

                const FP Nj       = pop[Sj] + pop[SVj] + pop[Ej] + pop[EVj] + pop[Ij] + pop[IVj] + pop[Rj] + pop[RVj];
                const FP inf_frac = (Nj > FP(0)) ? ((pop[Ij] + pop[IVj]) / Nj) : FP(0);
                const FP transmission_rate = (clustering_exp == FP(1)) ? inf_frac : std::pow(inf_frac, clustering_exp);
                sum += effective_contacts(i.get(), j.get()) * transmission_rate;
            }
            const FP foi_i = transmissibility_baseline * recovery_rate * season * sum + outside_foi;

            // Flows: S->E, S^V->E^V, E->I, E^V->I^V, I->R, I^V->R^V
            flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Exposed>(i)] =
                foi_i * y[Si];
            flows[Base::template get_flat_flow_index<InfectionState::SusceptibleVaccinated,
                                                     InfectionState::ExposedVaccinated>(i)] = foi_i * y[SVi];

            flows[Base::template get_flat_flow_index<InfectionState::Exposed, InfectionState::Infected>(i)] =
                infection_rate * y[Ei];
            flows[Base::template get_flat_flow_index<InfectionState::ExposedVaccinated,
                                                     InfectionState::InfectedVaccinated>(i)] = infection_rate * y[EVi];

            flows[Base::template get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>(i)] =
                recovery_rate * y[Ii];
            flows[Base::template get_flat_flow_index<InfectionState::InfectedVaccinated,
                                                     InfectionState::RecoveredVaccinated>(i)] = recovery_rate * y[IVi];
        }
    }
};

} // namespace oseirv
} // namespace mio
#endif
