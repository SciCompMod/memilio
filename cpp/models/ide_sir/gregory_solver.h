/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Anna Wendler
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

#include "ide_sir/parameters.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace isir
{

class GregorySolver
{
    using ParameterSet = Parameters;

public:
    GregorySolver(TimeSeries<ScalarType> pop, Parameters model_parameters, size_t gregory_order);

    ScalarType fixed_point_function(ScalarType s, TimeSeries<ScalarType> populations, ScalarType dt,
                                    Parameters model_parameters, ScalarType N);

    void compute_S(ScalarType s_guess, TimeSeries<ScalarType> populations, ScalarType dt, Parameters model_parameters,
                   ScalarType N, ScalarType tol = 1e-5, size_t max_iterations = 1000);

    void compute_S_deriv();
    void compute_I_and_R();

    TimeSeries<ScalarType> m_pop;

private:
    Parameters m_model_parameters;
    size_t m_gregory_order;
};

} // namespace isir
} // namespace mio
