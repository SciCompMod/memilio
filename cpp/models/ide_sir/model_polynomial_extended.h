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

#ifndef IDESIR_MODEL_POLY_H
#define IDESIR_MODEL_POLY_H

#include "ide_sir/infection_state.h"
#include "ide_sir/parameters.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"

namespace mio
{
namespace isir
{
class ModelPolynomialExtended
{
    // using ParameterSet = Parameters;

public:
    ModelPolynomialExtended(TimeSeries<ScalarType>&& populations_init, size_t gregory_order, ScalarType N);

    ScalarType get_totalpop() const;

    size_t get_gregory_order() const
    {
        return m_gregory_order;
    }

    ScalarType sum_part1_weight(size_t n, size_t j);
    ScalarType sum_part2_weight(size_t n, size_t j);

    void approx_single_integral(ScalarType dt, size_t t0_index, ScalarType alpha = 0.);
    void approx_double_integral(ScalarType dt, size_t t0_index, ScalarType alpha = 0.);

    // ---- Public parameters. ----
    TimeSeries<ScalarType> populations; ///< TimeSeries containing points of time and the corresponding number of
    // people in defined #InfectionState%s for every AgeGroup.
    TimeSeries<ScalarType> flows;

private:
    // ---- Private parameters. ----
    size_t m_gregory_order;
    ScalarType m_N;
};

} // namespace isir
} // namespace mio

#endif // IDESIR_MODEL_POLY_H
