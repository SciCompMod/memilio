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

        // Contact matrices
        const auto cm_h_expr =
            params.template get<ContactPatternsHealthy<FP>>().get_cont_freq_mat().get_matrix_at(SimulationTime<FP>(t));
        const auto cm_s_expr =
            params.template get<ContactPatternsSick<FP>>().get_cont_freq_mat().get_matrix_at(SimulationTime<FP>(t));

        // Get effective contact matrix B = H + m*S
        Eigen::MatrixX<FP> H = cm_h_expr;
        Eigen::MatrixX<FP> S = cm_s_expr;
        const FP m           = params.template get<SickMixingM<FP>>();
        Eigen::MatrixX<FP> B = H + m * S;

        // Normalization ν(m): prepare susceptibility scaling and next-generation matrix (NG)
        // sigma_i = Phi * SigmaByAge[i] (age-specific susceptibility scaled by a global factor)
        Eigen::VectorX<FP> sigma(num_groups);
        for (size_t i = 0; i < num_groups; ++i)
            sigma[(Eigen::Index)i] =
                params.template get<SigmaByAge<FP>>()[AgeGroup(i)] * params.template get<Phi<FP>>();
        // Sigma is a diagonal matrix applying susceptibility on the receiving (row) side of contacts.
        Eigen::MatrixX<FP> Sigma = sigma.asDiagonal();

        // NG = Sigma * B is the (unnormalized) next-generation matrix. Its spectral radius will be used
        // below to scale contacts B so that structure (mixing pattern) and transmissibility magnitude (BaselineTransmissibility)
        // are cleanly separated.
        Eigen::MatrixX<FP> NG = Sigma * B;

        Eigen::ComplexEigenSolver<Eigen::MatrixX<FP>> ces;
        ces.compute(NG);
        FP nu = ces.eigenvalues().cwiseAbs().maxCoeff();
        if (nu > FP(0.)) {
            B /= nu;
        }

        const FP Re =
            params.template get<BaselineTransmissibility<FP>>(); // baseline transmissibility scaling (R-like factor)
        const FP gamma   = params.template get<Gamma<FP>>(); // progression/recovery rate (1 / mean duration)
        const FP delta   = params.template get<SeasonalityAmplitude<FP>>(); // amplitude of seasonal modulation
        const FP tz      = params.template get<ShiftTZ<FP>>(); // base phase shift for seasonality
        const FP ts      = params.template get<ShiftTS<FP>>(); // additional phase shift (e.g. scenario-specific)
        const FP lambda0 = params.template get<OutsideFoI<FP>>(); // constant external force of infection
        const FP rho     = params.template get<
            ClusteringRho<FP>>(); // clustering exponent; rho<1 dampens, rho>1 amplifies prevalence effect

        const FP season = std::exp(delta * std::sin(FP(2) * std::numbers::pi_v<FP> * (t / FP(52.0) - tz + ts)));

        for (auto i : make_index_range((mio::AgeGroup)num_groups)) {
            // Flat indices
            const size_t Si  = this->populations.get_flat_index({i, InfectionState::Susceptible});
            const size_t SVi = this->populations.get_flat_index({i, InfectionState::SusceptibleVaccinated});
            const size_t Ei  = this->populations.get_flat_index({i, InfectionState::Exposed});
            const size_t EVi = this->populations.get_flat_index({i, InfectionState::ExposedVaccinated});
            const size_t Ii  = this->populations.get_flat_index({i, InfectionState::Infected});
            const size_t IVi = this->populations.get_flat_index({i, InfectionState::InfectedVaccinated});

            // λ_i(t) = Re*γ*season * sum_j B_{ij} * ((I_j+I^V_j)/N_j)^ρ + λ0
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
                const FP term     = (rho == FP(1)) ? inf_frac : std::pow(inf_frac, rho);
                sum += B(i.get(), j.get()) * term;
            }
            const FP lambda_i = Re * gamma * season * sum + lambda0;

            // Flows: S->E, S^V->E^V, E->I, E^V->I^V, I->R, I^V->R^V
            flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Exposed>(i)] =
                lambda_i * y[Si];
            flows[Base::template get_flat_flow_index<InfectionState::SusceptibleVaccinated,
                                                     InfectionState::ExposedVaccinated>(i)] = lambda_i * y[SVi];

            flows[Base::template get_flat_flow_index<InfectionState::Exposed, InfectionState::Infected>(i)] =
                gamma * y[Ei];
            flows[Base::template get_flat_flow_index<InfectionState::ExposedVaccinated,
                                                     InfectionState::InfectedVaccinated>(i)] = gamma * y[EVi];

            flows[Base::template get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>(i)] =
                gamma * y[Ii];
            flows[Base::template get_flat_flow_index<InfectionState::InfectedVaccinated,
                                                     InfectionState::RecoveredVaccinated>(i)] = gamma * y[IVi];
        }
    }
};

} // namespace oseirv
} // namespace mio
#endif
