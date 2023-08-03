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
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"

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
    *@brief After the simulation is finished and we get a resulting TimeSeries, this function uses the data to compute the reproduction number
    *at an arbitrary time
    *@param timeindex The time at which we want to compute the reproduction number
    *@param y The TimeSeries. We actually only use the number of Susceptibles at a given time
    *@returns The reproduction number in the seir model
    *@see The same functions in the model.h files of the secir and secirvvs models
    */

    ScalarType get_reproduction_number(Eigen::Index timeindex, const mio::TimeSeries<ScalarType>& y)
    { 
        if(!(timeindex < y.get_num_time_points())){
            mio::log_error("timeindex is out of range of the TimeSeries");
        }

        ScalarType TimeInfected = this->parameters.get<mio::oseir::TimeInfected>();

        ScalarType coeffStoE = this->parameters.get<mio::oseir::ContactPatterns>().get_matrix_at(timeindex)(0,0)*
                                this->parameters.get<mio::oseir::TransmissionProbabilityOnContact>()/
                                this->populations.get_total();

        return y.get_value(timeindex)[(Eigen::Index)mio::oseir::InfectionState::Susceptible] * TimeInfected * coeffStoE;
    }

     /**
    *@brief Compute the reproduction numbers for all time points based on the given TimeSeries.
    *@param y The TimeSeries containing the SEIR model data 
    *@returns vector containing all reproduction numbers
    */

    Eigen::VectorXd get_reproduction_numbers(const mio::TimeSeries<ScalarType>& y)
    {
        auto num_time_points = y.get_num_time_points();
        Eigen::VectorXd temp(num_time_points);
        for (int i = 0; i < num_time_points; i++) {
            temp[i] = get_reproduction_number(i, y);
        }
        return temp;
    }


    
};

} // namespace oseir
} // namespace mio

#endif // SEIR_MODEL_H
