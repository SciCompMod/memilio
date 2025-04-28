/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding
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
#ifndef MIO_MATH_TIME_SERIES_FUNCTOR_H
#define MIO_MATH_TIME_SERIES_FUNCTOR_H

#include "memilio/io/default_serialize.h"
#include "memilio/math/interpolation.h"
#include "memilio/utils/time_series.h"

#include <cassert>
#include <vector>

namespace mio
{

/**
 * @brief Type of a TimeSeriesFunctor.
 * The available types are:
 * - LinearInterpolation:
 *   - Requires at least one time point with exactly one value each. Time must be strictly monotic increasing.
 *   - Linearly interpolates between data points. Stays constant outside of provided data with first/last value.
 */
enum class TimeSeriesFunctorType
{
    LinearInterpolation,
};

template <class FP>
class TimeSeriesFunctor
{
public:
    /**
     * @brief Creates a functor using the given data.
     * Note the data requirements of the given functor type.
     * @param type The type of the functor.
     * @param table A list of time points, passed to the TimeSeries constructor.
     */
    TimeSeriesFunctor(TimeSeriesFunctorType type, const TimeSeries<FP>& data)
        : m_type(type)
        , m_data(data)
    {
        // data shape checks and preprocessing
        switch (m_type) {
        case TimeSeriesFunctorType::LinearInterpolation:
            // make sure data has the correct shape, i.e. a list of (time, value) pairs
            assert(m_data.get_num_time_points() > 0 && "Need at least one time point for LinearInterpolation.");
            assert(m_data.get_num_elements() == 1 && "LinearInterpolation requires exactly one value per time point.");
            assert(m_data.is_strictly_monotonic());
            break;
        default:
            assert(false && "Unhandled TimeSeriesFunctorType!");
            break;
        }
    }

    /**
     * @brief Creates a functor using the given table.
     * Note the data requirements of the given functor type.
     * @param type The type of the functor.
     * @param table A list of time points, passed to the TimeSeries constructor.
     */
    TimeSeriesFunctor(TimeSeriesFunctorType type, std::vector<std::vector<FP>>&& table)
        : TimeSeriesFunctor(type, TimeSeries<FP>{table})
    {
    }

    /**
     * @brief Creates a Zero functor.
     */
    TimeSeriesFunctor()
        : TimeSeriesFunctor(TimeSeriesFunctorType::LinearInterpolation, {{FP(0.0), FP(0.0)}})
    {
    }

    /**
     * @brief Function returning a scalar value.
     * @param time A scalar value.
     * @return A scalar value computed from data, depending on the functor's type.
     */
    FP operator()(FP time) const
    {
        switch (m_type) {
        case TimeSeriesFunctorType::LinearInterpolation:
            return linear_interpolation(time, m_data)[0];
        default:
            assert(false && "Unhandled TimeSeriesFunctorType!");
            return FP();
        }
    }

    /// This method is used by the default serialization feature.
    auto default_serialize()
    {
        return Members("TimeSeriesFunctor").add("type", m_type).add("data", m_data);
    }

private:
    TimeSeriesFunctorType m_type; ///< Determines what kind of functor this is, e.g. linear interpolation.
    TimeSeries<FP> m_data; ///< Data used by the functor to compute its values. Its shape depends on type.
};

} // namespace mio

#endif
