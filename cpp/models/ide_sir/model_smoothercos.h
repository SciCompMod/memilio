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
    ModelSmootherCos(TimeSeries<ScalarType>&& populations_init, TimeSeries<ScalarType>&& groundtruth_ts,
                     size_t fd_order, ScalarType damping_time, ScalarType damping, ScalarType cont_freq,
                     ScalarType smoother_window, ScalarType k = 1);

    ScalarType get_totalpop() const;

    size_t get_finite_difference_order() const
    {
        return m_finite_difference_order;
    }

    // Function that compute numerical derivative of some function.
    ScalarType compute_deriv_numerical(ScalarType dt, ScalarType current_time,
                                       std::function<ScalarType(ScalarType)> smoother_func);

    // Functions for computing/approximating smoothercosine.
    ScalarType smoothercos_via_contacts(ScalarType current_time);
    ScalarType smoothercos(ScalarType current_time);
    ScalarType smoothercos_deriv(ScalarType current_time);
    void approximate_smoothercos(ScalarType dt, ScalarType current_time);

    // Functions for computing/approximating smoothstep.
    ScalarType smoothstep(ScalarType current_time);
    ScalarType smoothstep_deriv(ScalarType current_time);
    void approximate_smoothstep(ScalarType dt, ScalarType current_time);

    // Functions for computing/approximating smoothstep where the resulting function is C2.
    ScalarType smoothstep_c2(ScalarType current_time);
    ScalarType smoothstep_c2_deriv(ScalarType current_time);
    void approximate_smoothstep_c2(ScalarType dt, ScalarType current_time);

    // Functions for computing/approximating smoothstep where the resulting function is C2.
    ScalarType sigmoid(ScalarType current_time);
    ScalarType sigmoid_smoother(ScalarType current_time);
    ScalarType sigmoid_smoother_deriv(ScalarType current_time);
    void approximate_sigmoid_smoother(ScalarType dt, ScalarType current_time);

    // Set groundtruth.
    void set_groundtruth(ScalarType current_time, std::string smoother_func_str);

    // ---- Public parameters. ----
    ParameterSet parameters{}; ///< ParameterSet of Model Parameters.
    TimeSeries<ScalarType> populations; ///< TimeSeries containing points of time and the corresponding number of
    // people in defined #InfectionState%s for every AgeGroup.
    TimeSeries<ScalarType> groundtruth;

private:
    // ---- Private parameters. ----
    size_t m_finite_difference_order;
    ScalarType m_damping_time;
    ScalarType m_damping;
    ScalarType m_cont_freq;
    ScalarType m_smoother_window;
    ScalarType m_k{1.}; /// Parameter of sigmoid function.
};

} // namespace isir
} // namespace mio

#endif // IDESIR_MODEL_H
