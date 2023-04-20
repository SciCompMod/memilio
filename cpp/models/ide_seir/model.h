/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke
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

#ifndef IDE_SEIR_H
#define IDE_SEIR_H

#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "ide_seir/parameters.h"
#include "ide_seir/infection_state.h"

#include <vector>

namespace mio
{
namespace iseir
{
class Model
{
    using Pa = ParametersBase;

public:
    /**
        * @brief Create an IDE SEIR model.
        *
        * @param[in, out] init TimeSeries with the initial values of the number of susceptibles at associated initial times.
        *   The time steps in this vector should be equidistant and equal to the time step used for the simulation. 
        *   A certain history of time steps and values for susceptibles is needed. 
        *   A warning is displayed if the condition is violated.
        *   Co be more precise, the first time point needs to be smaller than -(k-1)*TimeStep with 
        *       k=ceil((InfectiousTime + LatencyTime)/TimeStep).
        *   The last time point in this vector should be a time 0.
        * @param[in] dt_init The size of the time step used for numerical simulation.
        * @param[in] N_init The population of the considered region. 
        */
    Model(TimeSeries<double>&& init, double dt_init, int N_init, const Pa& Parameterset_init = Pa());

    /**
        * @brief Simulate the evolution of infection numbers with the given IDE SEIR model.
        *
        * The simulation is performed by solving the underlying model equation numerically. 
        * Here, an integro-differential equation is to be solved. The model parameters and the initial data are used.
        *
        * @param[in] t_max Last simulation day. 
        *   If the last point of time of the initial TimeSeries was 0, the simulation will be executed for t_max days.
        * @return The result of the simulation, stored in a TimeSeries with simulation time and 
        *       associated number of susceptibles.
        */
    TimeSeries<double> const& simulate(int t_max);

    /**
        * @brief Calculate the distribution of the population in E, I and, R based on the calculated values for S.
        *
        * The values are calculated using the average latency and infection time, not using model equations. 
        * The simulated values of S are used for this purpose, so the simulate() function should be called beforehand.
        * 
        * @return The result of the calculation stored in an TimeSeries. The TimeSeries contains the simulation time and an
        *   associated Vector with values for S, E, I, and R.
        */
    TimeSeries<double> const& calculate_EIR();

    /**
        * @brief Displays the results of the simulation.
        *
        * You can either output only the simulation times with the simulated values for S, or additionally the 
        * calculated numbers for E, I and R. In any case, the function simulate() should have been called before. 
        * If the values for E, I and R are to be displayed, the function calculate_EIR() must be executed beforehand.
        *
        * @param[in] calculated_SEIR If true, the calculated numbers for E, I and R are displayed 
        *    in addition to the results for S.
        */
    void print_result(bool calculated_SEIR = false) const;

    // Used Parameters for the simulation.
    Pa parameters{};

private:
    /**
        * @brief Density of the generalized beta distribution used for the function f_{beta} of the IDE SEIR model.
        *
        * @param[in] tau evaluation point of the generalized Beta distribution.
        * @param[in] p parameter p of the generalized Beta distribution.
        * @param[in] q parameter q of the generalized Beta distribution.
        * @result Evaluation of the generalized beta distribution at the given evaluation point.
        */
    double generalized_beta_distribution(double tau, double p = 3.0, double q = 10.0) const;

    /**
        * @brief Numerical differentiation of one compartment using a central difference quotient.
        * 
        * @param[in] ts_ide TimeSeries with the time steps already calculated. 
        *       Used as function values in numerical differentiation.
        * @param[in] idx Time index at which the numerical differentiation should be performed.
        * @param[in] compartment Compartment for which the numerical differentiation is to be performed.
        * @return Numerically approximated derivative of the function belonging to the compartment at the point t[idx].
        */
    double central_difference_quotient(TimeSeries<double> const& ts_ide, InfectionState compartment,
                                       Eigen::Index idx) const;

    /**
        * @brief Numerical integration of the inner integral of the integro-differential equation for the group S using
        *    a trapezoidal sum.
        * 
        * @param[in] idx Index of the point of time used in the inner integral.
        * @return Result of the numerical integration.
        */
    double num_integration_inner_integral(Eigen::Index idx) const;

    // TimeSeries containing points of time and the corresponding number of susceptibles.
    TimeSeries<double> m_result;
    // TimeSeries containing points of time and the corresponding number of susceptibles, exposed,
    // infected and recovered.
    TimeSeries<double> m_result_SEIR;

    // Timestep used for simulation.
    double m_dt{0};
    // Population of the considered region.
    int m_N{0};

    // Two Indices used for simulation.
    Eigen::Index m_k{0};
    Eigen::Index m_l{0};
};
} // namespace iseir
} // namespace mio
#endif
