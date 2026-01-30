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

/* This file conatains ModelMessinaExtendedDetailedInit whioch is the model that will be further investigated and developed.*/

#ifndef IDESIR_MODEL_H
#define IDESIR_MODEL_H

#include "ide_sir_analytical_S_deriv/infection_state.h"
#include "ide_sir_analytical_S_deriv/parameters.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace isir
{
// This class is a further extension of ModelMessinaExtended where use detailed initial conditions, i.e. we assume that
// we have knowledge of the compartment values on a time interval before the simulation start and not just at the simulation start.
class ModelAnalyticalSDeriv
{
    using ParameterSet = Parameters;

public:
    ModelAnalyticalSDeriv(TimeSeries<ScalarType>&& populations_init, size_t gregory_order,
                          size_t finite_difference_order);

    ScalarType get_totalpop() const;

    size_t get_gregory_order() const
    {
        return m_gregory_order;
    }

    size_t get_finite_difference_order() const
    {
        return m_finite_difference_order;
    }

    void set_tol_for_support_max(ScalarType new_tol)
    {
        m_tol = new_tol;
    }

    ScalarType compute_calctime(ScalarType dt, ScalarType tol = 1e-10) const
    {
        return parameters.get<TransitionDistributions>()[(size_t)InfectionTransition::InfectedToRecovered]
            .get_support_max(dt, tol);
    }

    ScalarType sum_part1_weight(size_t n, size_t j);
    ScalarType sum_part2_weight(size_t n, size_t j);

    // Returns the number of iterations needed in fixed point iteration.
    size_t compute_S(ScalarType s_init, ScalarType dt, ScalarType tol = 1e-14, size_t max_iterations = 100);

    ScalarType fixed_point_function(ScalarType s, ScalarType dt);

    void compute_S_deriv(ScalarType dt, size_t time_point_index);
    void compute_S_deriv(ScalarType dt);

    void compute_I_and_R(ScalarType dt, size_t time_point_index);
    void compute_I_and_R(ScalarType dt);

    void set_transitiondistribution_vector(ScalarType dt, ScalarType tmax, size_t t0_index = 0);
    void set_parameter_vectors(ScalarType dt, ScalarType tmax, size_t t0_index = 0);

    // ---- Public parameters. ----
    ParameterSet parameters{}; ///< ParameterSet of Model Parameters.
    TimeSeries<ScalarType> populations; ///< TimeSeries containing points of time and the corresponding number of
    // people in defined #InfectionState%s for every AgeGroup.
    TimeSeries<ScalarType> flows;

private:
    // ---- Private parameters. ----
    ScalarType m_N; ///< Vector containing the total population size of the considered region for every AgeGroup.
    size_t m_gregory_order;
    size_t m_finite_difference_order;
    std::vector<ScalarType> m_transitiondistribution_vector;
    std::vector<ScalarType> m_transmissionproboncontact_vector;
    std::vector<ScalarType> m_riskofinffromsymptomatic_vector;
    ScalarType m_calctime;
    ScalarType m_tol{1e-10};
};

} // namespace isir
} // namespace mio

#endif // IDESIR_MODEL_H
