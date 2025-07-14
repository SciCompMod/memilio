/*
* Copyright (C) 2020-2025 MEmilio
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
#ifndef SEIR_MODEL_H
#define SEIR_MODEL_H

#include "memilio/compartments/flow_model.h"
#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/math/interpolation.h"
#include "memilio/utils/time_series.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"

GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wshadow")
#include <Eigen/Dense>
GCC_CLANG_DIAGNOSTIC(pop)

namespace mio
{
namespace oseir
{

/********************
 * define the model *
 ********************/

// clang-format off
using Flows = TypeList<Flow<InfectionState::Susceptible, InfectionState::Exposed>,
                       Flow<InfectionState::Exposed,     InfectionState::Infected>,
                       Flow<InfectionState::Infected,    InfectionState::Recovered>>;
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
        : Base(Populations({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorX<FP>> pop, Eigen::Ref<const Eigen::VectorX<FP>> y, FP t,
                   Eigen::Ref<Eigen::VectorX<FP>> flows) const override
    {
        const Index<AgeGroup> age_groups = reduce_index<Index<AgeGroup>>(this->populations.size());
        const auto& params               = this->parameters;

        for (auto i : make_index_range(age_groups)) {
            const size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            const size_t Ei = this->populations.get_flat_index({i, InfectionState::Exposed});
            const size_t Ii = this->populations.get_flat_index({i, InfectionState::Infected});

            for (auto j : make_index_range(age_groups)) {
                const size_t Sj = this->populations.get_flat_index({j, InfectionState::Susceptible});
                const size_t Ej = this->populations.get_flat_index({j, InfectionState::Exposed});
                const size_t Ij = this->populations.get_flat_index({j, InfectionState::Infected});
                const size_t Rj = this->populations.get_flat_index({j, InfectionState::Recovered});

                const FP Nj        = pop[Sj] + pop[Ej] + pop[Ij] + pop[Rj];
                const FP divNj     = (Nj < Limits<FP>::zero_tolerance()) ? FP(0.0) : FP(1.0 / Nj);
                const FP coeffStoE = params.template get<ContactPatterns<FP>>().get_cont_freq_mat().get_matrix_at(
                                         SimulationTime<FP>(t))(i.get(), j.get()) *
                                     params.template get<TransmissionProbabilityOnContact<FP>>()[i] * divNj;

                flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Exposed>(i)] +=
                    coeffStoE * y[Si] * pop[Ij];
            }
            flows[Base::template get_flat_flow_index<InfectionState::Exposed, InfectionState::Infected>(i)] =
                (1.0 / params.template get<TimeExposed<FP>>()[i]) * y[Ei];
            flows[Base::template get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>(i)] =
                (1.0 / params.template get<TimeInfected<FP>>()[i]) * y[Ii];
        }
    }

    /**
    *@brief Computes the reproduction number at a given index time of the Model output obtained by the Simulation.
    *@param t_idx The index time at which the reproduction number is computed.
    *@param y The TimeSeries obtained from the Model Simulation.
    *@returns The computed reproduction number at the provided index time.
    */
    IOResult<FP> get_reproduction_number(size_t t_idx, const mio::TimeSeries<FP>& y)
    {
        if (!(t_idx < static_cast<size_t>(y.get_num_time_points()))) {
            return mio::failure(mio::StatusCode::OutOfRange, "t_idx is not a valid index for the TimeSeries");
        }

        auto const& params = this->parameters;

        const size_t num_groups                    = (size_t)params.get_num_groups();
        constexpr size_t num_infected_compartments = 2;
        const size_t total_infected_compartments   = num_infected_compartments * num_groups;

        ContactMatrixGroup<FP> const& contact_matrix = params.template get<ContactPatterns<ScalarType>>();

        Eigen::MatrixX<FP> F = Eigen::MatrixX<FP>::Zero(total_infected_compartments, total_infected_compartments);
        Eigen::MatrixX<FP> V = Eigen::MatrixX<FP>::Zero(total_infected_compartments, total_infected_compartments);

        for (auto i = AgeGroup(0); i < AgeGroup(num_groups); i++) {
            size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            for (auto j = AgeGroup(0); j < AgeGroup(num_groups); j++) {

                const ScalarType Nj    = this->populations.get_group_total(j);
                const ScalarType divNj = (Nj < 1e-12) ? 0.0 : 1.0 / Nj;

                FP coeffStoE = contact_matrix.get_matrix_at(SimulationTime<FP>(y.get_time(t_idx)))(i.get(), j.get()) *
                               params.template get<TransmissionProbabilityOnContact<FP>>()[i] * divNj;
                F((size_t)i, (size_t)j + num_groups) = coeffStoE * y.get_value(t_idx)[Si];
            }

            FP T_Ei                                           = params.template get<mio::oseir::TimeExposed<FP>>()[i];
            FP T_Ii                                           = params.template get<mio::oseir::TimeInfected<FP>>()[i];
            V((size_t)i, (size_t)i)                           = 1.0 / T_Ei;
            V((size_t)i + num_groups, (size_t)i)              = -1.0 / T_Ei;
            V((size_t)i + num_groups, (size_t)i + num_groups) = 1.0 / T_Ii;
        }

        V = V.inverse();

        //Compute F*V
        Eigen::MatrixX<FP> NextGenMatrix(total_infected_compartments, total_infected_compartments);
        NextGenMatrix.noalias() = F * V;

        Eigen::MatrixXd NextGenMatrix_dbl = NextGenMatrix.unaryExpr([](const FP& x) {
            return ad::value(x);
        });

        Eigen::ComplexEigenSolver<Eigen::MatrixXd> ces;
        ces.compute(NextGenMatrix_dbl);

        const Eigen::VectorXcd eigvals_complex = ces.eigenvalues();
        Eigen::VectorXd eigvals_abs(eigvals_complex.size());
        for (int i = 0; i < eigvals_complex.size(); ++i) {
            eigvals_abs[i] = std::abs(eigvals_complex[i]);
        }

        FP rho = static_cast<FP>(eigvals_abs.maxCoeff());
        return mio::success(rho);
    }

    /**
    *@brief Computes the reproduction number for all time points of the Model output obtained by the Simulation.
    *@param y The TimeSeries obtained from the Model Simulation.
    *@returns vector containing all reproduction numbers
    */
    Eigen::VectorX<FP> get_reproduction_numbers(const mio::TimeSeries<FP>& y)
    {
        auto num_time_points = y.get_num_time_points();
        Eigen::VectorX<FP> temp(num_time_points);
        for (size_t i = 0; i < static_cast<size_t>(num_time_points); i++) {
            temp[i] = get_reproduction_number(i, y).value();
        }
        return temp;
    }

    /**
    *@brief Computes the reproduction number at a given time point of the Model output obtained by the Simulation. If the particular time point is not inside the output, a linearly interpolated value is returned.
    *@param t_value The time point at which the reproduction number is computed.
    *@param y The TimeSeries obtained from the Model Simulation.
    *@returns The computed reproduction number at the provided time point, potentially using linear interpolation.
    */
    IOResult<FP> get_reproduction_number(FP t_value, const mio::TimeSeries<FP>& y)
    {
        if (t_value < y.get_time(0) || t_value > y.get_last_time()) {
            return mio::failure(mio::StatusCode::OutOfRange,
                                "Cannot interpolate reproduction number outside computed horizon of the TimeSeries");
        }

        if (t_value == y.get_time(0)) {
            return mio::success(get_reproduction_number((size_t)0, y).value());
        }

        auto times = std::vector<FP>(y.get_times().begin(), y.get_times().end());

        auto time_late = std::distance(times.begin(), std::lower_bound(times.begin(), times.end(), t_value));

        FP y1 = get_reproduction_number(static_cast<size_t>(time_late - 1), y).value();
        FP y2 = get_reproduction_number(static_cast<size_t>(time_late), y).value();

        auto result = linear_interpolation(t_value, y.get_time(time_late - 1), y.get_time(time_late), y1, y2);
        return mio::success(static_cast<FP>(result));
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

} // namespace oseir
} // namespace mio

#endif // SEIR_MODEL_H
