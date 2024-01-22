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

#include "memilio/compartments/flow_model.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/utils/type_list.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/io/io.h"
#include "memilio/math/interpolation.h"
#include "memilio/utils/time_series.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include <algorithm>
#include <iterator>

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
template<typename FP=double>
class Model : public FlowModel<InfectionState, Populations<InfectionState>, Parameters<FP>, Flows, FP>
{
    using Base = FlowModel<InfectionState, mio::Populations<InfectionState>, Parameters<FP>, Flows, FP>;

public:
    Model()
        : Base(typename Base::Populations({InfectionState::Count}, 0.),typename Base::ParameterSet())
    {
    }

    void get_flows(Eigen::Ref<const Eigen::Matrix<FP,Eigen::Dynamic,1>> pop, Eigen::Ref<const Eigen::Matrix<FP,Eigen::Dynamic,1>> y, FP t,
                   Eigen::Ref<Eigen::Matrix<FP,Eigen::Dynamic,1>> flows) const override
    {
        auto& params     = this->parameters;
        double coeffStoE = params.template get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                           params.template get<TransmissionProbabilityOnContact<FP>>() / this->populations.get_total();

        flows[Base::template get_flat_flow_index<InfectionState::Susceptible, InfectionState::Exposed>()] =
            coeffStoE * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected];
        flows[Base::template get_flat_flow_index<InfectionState::Exposed, InfectionState::Infected>()] =
            (1.0 / params.template get<TimeExposed<FP>>()) * y[(size_t)InfectionState::Exposed];
        flows[Base::template get_flat_flow_index<InfectionState::Infected, InfectionState::Recovered>()] =
            (1.0 / params.template get<TimeInfected<FP>>()) * y[(size_t)InfectionState::Infected];
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

        ScalarType TimeInfected = this->parameters.template get<mio::oseir::TimeInfected>();

        ScalarType coeffStoE = this->parameters.template get<mio::oseir::ContactPatterns>().get_matrix_at(
                                   y.get_time(static_cast<Eigen::Index>(t_idx)))(0, 0) *
                               this->parameters.template get<mio::oseir::TransmissionProbabilityOnContact>() /
                               this->populations.get_total();

        ScalarType result =
            y.get_value(static_cast<Eigen::Index>(t_idx))[(Eigen::Index)mio::oseir::InfectionState::Susceptible] *
            TimeInfected * coeffStoE;

        return mio::success(result);
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
