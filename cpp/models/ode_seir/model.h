/* 
* Copyright (C) 2020-2024 MEmilio
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

#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/utils/type_list.h"
#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/io/io.h"
#include "memilio/math/interpolation.h"
#include "memilio/utils/time_series.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include <algorithm>
#include <iterator>
#include "memilio/epidemiology/age_group.h"
#include "memilio/compartments/flow_simulation.h"

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

class Model : public FlowModel<InfectionState, Populations<AgeGroup, InfectionState>, Parameters, Flows>
{
    using Base = FlowModel<InfectionState, mio::Populations<AgeGroup, InfectionState>, Parameters, Flows>;

public:
    Model(int num_agegroups)
        : Base(Populations({AgeGroup(num_agegroups), InfectionState::Count}), ParameterSet(AgeGroup(num_agegroups)))
    {
    }

    void get_flows(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                   Eigen::Ref<Eigen::VectorXd> flows) const override
    {
        auto& params                             = this->parameters;
        AgeGroup n_agegroups                     = params.get_num_groups();
        ContactMatrixGroup const& contact_matrix = params.get<ContactPatterns>();

        for (auto i = AgeGroup(0); i < n_agegroups; i++) {
            size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            size_t Ei = this->populations.get_flat_index({i, InfectionState::Exposed});
            size_t Ii = this->populations.get_flat_index({i, InfectionState::Infected});

            for (auto j = AgeGroup(0); j < n_agegroups; j++) {

                size_t Sj = this->populations.get_flat_index({j, InfectionState::Susceptible});
                size_t Ej = this->populations.get_flat_index({j, InfectionState::Exposed});
                size_t Ij = this->populations.get_flat_index({j, InfectionState::Infected});
                size_t Rj = this->populations.get_flat_index({j, InfectionState::Recovered});

                double Nj    = pop[Sj] + pop[Ej] + pop[Ij] + pop[Rj];
                double divNj = 1.0 / Nj;

                double coeffStoE = contact_matrix.get_matrix_at(t)(static_cast<Eigen::Index>((size_t)i),
                                                                   static_cast<Eigen::Index>((size_t)j)) *
                                   params.get<TransmissionProbabilityOnContact>()[i] * divNj;

                double dummy_S = y[Si] * y[Ij] * coeffStoE;

                flows[get_flat_flow_index<InfectionState::Susceptible, InfectionState::Exposed>({i})] += dummy_S;
            }
            flows[get_flat_flow_index<InfectionState::Exposed, InfectionState::Infected>({i})] =
                (1.0 / params.get<TimeExposed>()[i]) * y[Ei];
            flows[get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>({i})] =
                (1.0 / params.get<TimeInfected>()[i]) * y[Ii];
        }
    }

    /**
    *@brief Computes the reproduction number at a given index time of the Model output obtained by the Simulation.
    *@param t_idx The index time at which the reproduction number is computed.
    *@param y The TimeSeries obtained from the Model Simulation.
    *@returns The computed reproduction number at the provided index time.
    */
    IOResult<ScalarType> get_reproduction_number(size_t t_idx, const mio::TimeSeries<ScalarType>& y)
    {
        if (!(t_idx < static_cast<size_t>(y.get_num_time_points()))) {
            return mio::failure(mio::StatusCode::OutOfRange, "t_idx is not a valid index for the TimeSeries");
        }

        auto const& params = this->parameters;

        const size_t num_groups                  = (size_t)params.get_num_groups();
        const size_t num_infected_compartments   = 2;
        const size_t total_infected_compartments = num_infected_compartments * num_groups;

        ContactMatrixGroup const& contact_matrix = params.get<ContactPatterns>();

        Eigen::MatrixXd F = Eigen::MatrixXd::Zero(total_infected_compartments, total_infected_compartments);
        Eigen::MatrixXd V = Eigen::MatrixXd::Zero(total_infected_compartments, total_infected_compartments);

        for (auto i = AgeGroup(0); i < AgeGroup(num_groups); i++) {
            size_t Si = this->populations.get_flat_index({i, InfectionState::Susceptible});
            for (auto j = AgeGroup(0); j < AgeGroup(num_groups); j++) {

                size_t Sj = this->populations.get_flat_index({j, InfectionState::Susceptible});
                size_t Ej = this->populations.get_flat_index({j, InfectionState::Exposed});
                size_t Ij = this->populations.get_flat_index({j, InfectionState::Infected});
                size_t Rj = this->populations.get_flat_index({j, InfectionState::Recovered});

                std::cout << " y = " << y.get_value(t_idx) << std::endl;
                double Nj =
                    y.get_value(t_idx)[Sj] + y.get_value(t_idx)[Ej] + y.get_value(t_idx)[Ij] + y.get_value(t_idx)[Rj];
                double divNj = 1.0 / Nj;

                double coeffStoE = contact_matrix.get_matrix_at(y.get_time(t_idx))(
                                       static_cast<Eigen::Index>((size_t)i), static_cast<Eigen::Index>((size_t)j)) *
                                   params.get<TransmissionProbabilityOnContact>()[i] * divNj;
                F((size_t)i, (size_t)j + num_groups) = coeffStoE * y.get_value(t_idx)[Si];
            }

            double T_Ei                                       = params.get<mio::oseir::TimeExposed>()[i];
            double T_Ii                                       = params.get<mio::oseir::TimeInfected>()[i];
            V((size_t)i, (size_t)i)                           = 1.0 / T_Ei;
            V((size_t)i + num_groups, (size_t)i)              = -1.0 / T_Ei;
            V((size_t)i + num_groups, (size_t)i + num_groups) = 1.0 / T_Ii;
        }

        V = V.inverse();

        Eigen::MatrixXd NextGenMatrix = Eigen::MatrixXd::Zero(total_infected_compartments, total_infected_compartments);
        NextGenMatrix                 = F * V;

        //Compute the largest eigenvalue in absolute value
        Eigen::ComplexEigenSolver<Eigen::MatrixXd> ces;

        ces.compute(NextGenMatrix);
        const Eigen::VectorXcd eigen_vals = ces.eigenvalues();

        Eigen::VectorXd eigen_vals_abs;
        eigen_vals_abs.resize(eigen_vals.size());

        for (int i = 0; i < eigen_vals.size(); i++) {
            eigen_vals_abs[i] = std::abs(eigen_vals[i]);
        }
        return mio::success(eigen_vals_abs.maxCoeff());
    }

    /**
    *@brief Computes the reproduction number for all time points of the Model output obtained by the Simulation.
    *@param y The TimeSeries obtained from the Model Simulation.
    *@returns vector containing all reproduction numbers
    */
    Eigen::VectorXd get_reproduction_numbers(const mio::TimeSeries<ScalarType>& y)
    {
        auto num_time_points = y.get_num_time_points();
        Eigen::VectorXd temp(num_time_points);
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
    IOResult<ScalarType> get_reproduction_number(ScalarType t_value, const mio::TimeSeries<ScalarType>& y)
    {
        if (t_value < y.get_time(0) || t_value > y.get_last_time()) {
            return mio::failure(mio::StatusCode::OutOfRange,
                                "Cannot interpolate reproduction number outside computed horizon of the TimeSeries");
        }

        if (t_value == y.get_time(0)) {
            return mio::success(get_reproduction_number((size_t)0, y).value());
        }

        auto times = std::vector<ScalarType>(y.get_times().begin(), y.get_times().end());

        auto time_late = std::distance(times.begin(), std::lower_bound(times.begin(), times.end(), t_value));

        ScalarType y1 = get_reproduction_number(static_cast<size_t>(time_late - 1), y).value();
        ScalarType y2 = get_reproduction_number(static_cast<size_t>(time_late), y).value();

        auto result = linear_interpolation(t_value, y.get_time(time_late - 1), y.get_time(time_late), y1, y2);
        return mio::success(static_cast<ScalarType>(result));
    }
};

} // namespace oseir
} // namespace mio

#endif // SEIR_MODEL_H
