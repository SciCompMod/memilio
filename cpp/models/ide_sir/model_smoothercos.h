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

#ifndef IDESIR_MODEL_SMOOTHERCOS_H
#define IDESIR_MODEL_SMOOTHERCOS_H

#include "ide_sir/infection_state.h"
#include "ide_sir/parameters.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace isir
{
// This class is a further extension of ModelMessinaExtended where use detailed initial conditions, i.e. we assume that
// we have knowledge of the compartment values on a time interval before the simulation start and not just at the simulation start.
class ModelSmootherCos
{
    using ParameterSet = Parameters;

public:
    ModelSmootherCos(TimeSeries<ScalarType>&& populations_init, size_t fd_order, ScalarType damping_time,
                     ScalarType damping, ScalarType cont_freq, ScalarType smoother_window);

    ScalarType get_totalpop() const;

    size_t get_finite_difference_order() const
    {
        return m_finite_difference_order;
    }

    // Returns the number of iterations needed in fixed point iteration.
    void compute_S();
    void compute_S_deriv(ScalarType dt, size_t j);
    void compute_S_deriv(ScalarType dt);
    void compute_S_deriv_analytical(ScalarType current_time);

    ScalarType smoothstep(ScalarType current_time);
    void compute_S_smoothstep(ScalarType current_time);
    void compute_smoothstep_deriv_analytical(ScalarType current_time);
    void compute_smoothstep_deriv_numerical(ScalarType current_time, ScalarType dt);

    // ---- Public parameters. ----
    ParameterSet parameters{}; ///< ParameterSet of Model Parameters.
    TimeSeries<ScalarType> populations; ///< TimeSeries containing points of time and the corresponding number of
    // people in defined #InfectionState%s for every AgeGroup.

private:
    // ---- Private parameters. ----
    size_t m_finite_difference_order;
    ScalarType m_damping_time;
    ScalarType m_damping;
    ScalarType m_cont_freq;
    ScalarType m_smoother_window;
};

} // namespace isir
} // namespace mio

#endif // IDESIR_MODEL_H
