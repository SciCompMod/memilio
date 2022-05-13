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
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/parameter_set.h"

#include <vector>

namespace mio
{
namespace iseir
{

    enum class iSeirInfType
    {
        S,
        E,
        I,
        R,
        Count = 4
    };

    /**************************
 * define some parameters *
 **************************/

    struct LatencyTime {
        using Type = double;
        static constexpr Type get_default()
        {
            return 3.3;
        }
    };
    struct InfectiousTime {
        using Type = double;
        static constexpr Type get_default()
        {
            return 8.2;
        }
    };
    struct TransmissionRisk {
        using Type = double;
        static constexpr Type get_default()
        {
            return 0.1;
        }
    };
    struct ContactFrequency {
        using Type = UncertainContactMatrix;
        static Type get_default()
        {
            ContactMatrixGroup contact_matrix = ContactMatrixGroup(1, 1);
            contact_matrix[0]                 = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10.));
            contact_matrix[0].add_damping(0.7, mio::SimulationTime(10.));
            return Type(contact_matrix);
        }
    };

    using SeirIdeParameters = ParameterSet<TransmissionRisk, LatencyTime, InfectiousTime, ContactFrequency>;

    class IdeModel
    {
        using Pa = SeirIdeParameters;

    public:
        /**
    * Constructor of IDEModel
    * @param init TimeSeries containing the initial (time, quantity of susceptible) values
    *   the initial time steps should be equidistant with distance dt_init; 
    *   the first point of time should be earlier than -(k-1)*TimeStep with k=ceil((InfectiousTime + LatencyTime)/TimeStep)
    *   the last point of time should be 0
    * @param dt_init size of the time step used for numerical integration.
    * @param N_init Population of the considered region. 
    */
        IdeModel(TimeSeries<double>&& init, double dt_init, int N_init, Pa Parameterset_init = Pa());

        /**
        * @brief Simulates the course of infection with the given IDE model.
        * @param t_max last simulation day
        * @return result of the simulation, stored in a TimeSeries with simulation time and associated number of susceptibles (S).
        */
        TimeSeries<double> const& simulate(int t_max);

        /**
        * @brief Calculate the distribution of the population in E,I and R based on the calculated values for S.
        *
        * Average values are calculated using the average latency and infectious time and not the distributions
        * actually used by the model. (because they are not explicitly calculable)
        * @return result of the calculation stored in an TimeSeries. The TimeSeries contains the simulation time and an
        * associated Vector with values for S, E, I, and R
        */
        TimeSeries<double> const& calculate_EIR();

        /**
        * print simulation result
        * @param calculated_SEIR if the values for E,I and R are calculated before with function calculate_EIR(),
        * one may want to print out these values too.
        */
        void print_result(bool calculated_SEIR = false) const;

        Pa m_parameters{};

    private:
        /**
        * Density of the generalized beta distribution used for the function f_{beta} of the IDGL model.
        *
        * @param tau evaluation point of the generalized Beta distribution
        * @param p parameter p of the generalized Beta distribution
        * @param q parameter q of the generalized Beta distribution
        */
        double generalized_beta_distribution(double tau, double p = 3.0, double q = 10.0) const;

        /**
        * Numerical differentiation of one compartment using a central difference quotient.
        * @param idx Time index at which the numerical differentiation should be performed.
        * @param compartment Index of the compartment in m_result for which the numerical differentiation is to be performed.
        * @return compartment'(t[idx]) numerically approximated
        */
        double central_difference_quotient(TimeSeries<double> const& ts_ide, iSeirInfType compartment,
                                           Eigen::Index idx) const;

        /**
        * Numerical integration of the inner integral of the integro-differential equation for the group S using a trapezoidal sum.
        * @param idx Index of the time in the inner integral
        * @return Result of the numerical integration
        */
        double num_integration_inner_integral(Eigen::Index idx) const;

        // vector containing one time Step per entry stored in an Eigen Vector (time, number of Susceptible at time t, R0t)
        TimeSeries<double> m_result;
        TimeSeries<double> m_result_SEIR = TimeSeries<double>(4);

        double m_dt{0};
        int m_N{0};

        Eigen::Index m_k{0};
        Eigen::Index m_l{0};
    };
} // namespace iseir
} // namespace mio
#endif