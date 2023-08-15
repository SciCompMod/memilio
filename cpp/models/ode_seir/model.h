/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#include "memilio/compartments/compartmentalmodel.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/io/io.h"
#include "memilio/math/interpolation.h"
#include "memilio/utils/time_series.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include <math.h>

namespace mio
{
namespace oseir
{

/********************
    * define the model *
    ********************/

class Model : public CompartmentalModel<InfectionState, Populations<InfectionState>, Parameters>
{
    using Base = CompartmentalModel<InfectionState, mio::Populations<InfectionState>, Parameters>;

public:
    Model()
        : Base(Populations({InfectionState::Count}, 0.), ParameterSet())
    {
    }

    void get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop, Eigen::Ref<const Eigen::VectorXd> y, double t,
                         Eigen::Ref<Eigen::VectorXd> dydt) const override
    {
        auto& params     = this->parameters;
        double coeffStoE = params.get<ContactPatterns>().get_matrix_at(t)(0, 0) *
                           params.get<TransmissionProbabilityOnContact>() / populations.get_total();

        dydt[(size_t)InfectionState::Susceptible] =
            -coeffStoE * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected];
        dydt[(size_t)InfectionState::Exposed] =
            coeffStoE * y[(size_t)InfectionState::Susceptible] * pop[(size_t)InfectionState::Infected] -
            (1.0 / params.get<TimeExposed>()) * y[(size_t)InfectionState::Exposed];
        dydt[(size_t)InfectionState::Infected] =
            (1.0 / params.get<TimeExposed>()) * y[(size_t)InfectionState::Exposed] -
            (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected];
        dydt[(size_t)InfectionState::Recovered] =
            (1.0 / params.get<TimeInfected>()) * y[(size_t)InfectionState::Infected];
    }

    /**
    *@brief Computes the reproduction number at a given index time of the Model output obtained by the Simulation.
    *@param t_idx The index time at which the reproduction number is computed.
    *@param y The TimeSeries obtained from the Model Simulation.
    *@returns The computed reproduction number at the provided index time.
    */
    IOResult<double> get_reproduction_number(size_t t_idx, const mio::TimeSeries<ScalarType>& y) noexcept
    {
        if (!(0 <= t_idx && t_idx < static_cast<size_t>(y.get_num_time_points()))) {
            return mio::failure(mio::StatusCode::OutOfRange, "t_idx is not a valid index for the TimeSeries");
        }

        ScalarType TimeInfected = this->parameters.get<mio::oseir::TimeInfected>();

        ScalarType coeffStoE = this->parameters.get<mio::oseir::ContactPatterns>().get_matrix_at(t_idx)(0, 0) *
                               this->parameters.get<mio::oseir::TransmissionProbabilityOnContact>() /
                               this->populations.get_total();

        double result =
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
        for (int i = 0; i < num_time_points; i++) {
            temp[i] = get_reproduction_number(i, y).value();
        }
        return temp;
    }

    /**
    * @brief interpolate reproduction numbers at freely chosen time value that lies in between the time points of the given time series.
    * values are linearly interpolated from their immediate neighbors from the old time points.
    * @param y TimeSeries whose reproduction numbers we want to interpolate
    * @param t_value double time at which we approximate the reproduction number
    * @return interpolated reproduction number at given time value
    */
    IOResult<double> interpolate_reproduction_numbers(const mio::TimeSeries<ScalarType>& y, const double t_value)
    {
        if (t_value < -1 || t_value > y.get_num_time_points()) {
            return mio::failure(
                mio::StatusCode::OutOfRange,
                "Cannot interpolate reproduction number too far outside computed horizon of the TimeSeries");
        }
        if (t_value < 0) {
            return get_reproduction_number(0, y);
        }
        if (t_value > y.get_num_time_points() - 1) {
            return get_reproduction_number(y.get_num_time_points() - 1, y);
        }
        double y1 = get_reproduction_number(floor(t_value), y).value();
        double y2 = get_reproduction_number(ceil(t_value), y).value();

        double result = linear_interpolation(t_value, floor(t_value), ceil(t_value), y1, y2);

        return mio::success(result);
    }
};

} // namespace oseir
} // namespace mio

#endif // SEIR_MODEL_H
