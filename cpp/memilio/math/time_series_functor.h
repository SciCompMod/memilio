/* 
* Copyright (C) 2020-2024 MEmilio
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

#include "memilio/io/auto_serialize.h"
#include "memilio/math/interpolation.h"
#include "memilio/utils/time_series.h"

#include <algorithm>
#include <cassert>
#include <vector>

namespace mio
{

template <class FP>
class TimeSeriesFunctor
{
public:
    enum Type
    {
        Zero,
        LinearInterpolation,
    };

    TimeSeriesFunctor(Type type, const TimeSeries<FP>& data)
        : m_type(type)
        , m_data(data)
    {
        // data preprocessing
        switch (m_type) {
        case Type::Zero:
            // no preprocessing needed
            break;
        case Type::LinearInterpolation:
            // make sure data has the correct shape, i.e. a list of (time, value) pairs
            assert(m_data.get_num_time_points() > 0 && "Need at least one time point for LinearInterpolation.");
            assert(m_data.get_num_elements() == 1 && "LinearInterpolation requires exactly one value per time point.");
            assert(m_data.is_sorted());
        }
    }

    TimeSeriesFunctor(Type type, std::vector<std::vector<FP>>&& table)
        : TimeSeriesFunctor(type, TimeSeries<FP>{table})
    {
    }

    TimeSeriesFunctor()
        : TimeSeriesFunctor(Type::Zero, TimeSeries<FP>{0})
    {
    }

    FP operator()(FP time) const
    {
        FP value = 0.0;
        switch (m_type) {
        case Type::Zero:
            // value is explicitly zero-initialized
            break;
        case Type::LinearInterpolation:
            // find next time point in m_data (strictly) after time
            auto tp_range      = m_data.get_times();
            const auto next_tp = std::upper_bound(tp_range.begin(), tp_range.end(), time, [](auto&& t, auto&& tp) {
                return t < tp;
            });
            if (next_tp == tp_range.begin()) { // time is before first data point
                value = m_data.get_value(0)[0];
            }
            else if (next_tp == tp_range.end()) { // time is past last data point
                value = m_data.get_last_value()[0];
            }
            else { // time is in between data points
                const auto i = next_tp - tp_range.begin();
                value        = linear_interpolation(time, m_data.get_time(i - 1), m_data.get_time(i),
                                                    m_data.get_value(i - 1)[0], m_data.get_value(i)[0]);
            }
            break;
        }
        return value;
    }

    /// This method is used by the auto-serialization feature.
    auto auto_serialize()
    {
        return make_auto_serialization("TimeSeriesFunctor", NVP("type", m_type), NVP("data", m_data));
    }

private:
    Type m_type;
    TimeSeries<FP> m_data;
};

} // namespace mio

#endif
