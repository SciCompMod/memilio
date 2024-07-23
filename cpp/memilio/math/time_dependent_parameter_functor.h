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
#ifndef MIO_MATH_TIME_DEPENDENT_PARAMETER_FUNCTOR_H
#define MIO_MATH_TIME_DEPENDENT_PARAMETER_FUNCTOR_H

#include "memilio/config.h"
#include "memilio/io/auto_serialize.h"
#include "memilio/math/interpolation.h"

#include <algorithm>
#include <vector>

namespace mio
{

class TimeDependentParameterFunctor
{
public:
    enum class Type
    {
        Zero,
        LinearInterpolation,
    };

    using DataType = std::vector<std::vector<ScalarType>>;
    TimeDependentParameterFunctor(Type type, const DataType& data)
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
            assert(m_data.size() > 0);
            assert(std::all_of(m_data.begin(), m_data.end(), [](auto&& a) {
                return a.size() == 2;
            }));
            // sort by time
            std::sort(m_data.begin(), m_data.end(), [](auto&& a, auto&& b) {
                return a[0] < b[0];
            });
        }
    }

    TimeDependentParameterFunctor()
        : TimeDependentParameterFunctor(Type::Zero, {})
    {
    }

    ScalarType operator()(ScalarType time) const
    {
        ScalarType value = 0.0;
        switch (m_type) {
        case Type::Zero:
            // value is explicitly zero-initialized
            break;
        case Type::LinearInterpolation:
            // find next time point in m_data (strictly) after time
            const auto next_tp = std::upper_bound(m_data.begin(), m_data.end(), time, [](auto&& t, auto&& tp) {
                return t < tp[0];
            });
            if (next_tp == m_data.begin()) { // time is before first data point
                return m_data.front()[1];
            }
            if (next_tp == m_data.end()) { // time is past last data point
                return m_data.back()[1];
            }
            const auto tp = next_tp - 1;
            value         = linear_interpolation(time, (*tp)[0], (*next_tp)[0], (*tp)[1], (*next_tp)[1]);
            break;
        }
        return value;
    }

    /// This method is used by the auto-serialization feature.
    auto auto_serialize()
    {
        return make_auto_serialization("TimeDependentParameterFunctor", NVP("type", m_type), NVP("data", m_data));
    }

private:
    Type m_type;
    DataType m_data;
};

} // namespace mio

#endif
