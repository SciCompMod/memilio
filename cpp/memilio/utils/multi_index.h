/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: René Schmieding
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
#ifndef MIO_UTILS_MULTI_INDEX_H_
#define MIO_UTILS_MULTI_INDEX_H_

#include "memilio/utils/index.h"

namespace mio
{

template <class MultiIndex>
class MultiIndexRange
{
public:
    MultiIndexRange(MultiIndex dimensions)
        : m_dimensions(dimensions)
        , m_end(make_end_index())
    {
    }

    class MultiIndexIterator
    {
    public:
        using value_type = MultiIndex;
        using reference  = const MultiIndex&;

        explicit MultiIndexIterator(value_type index, reference dimensions)
            : m_index(index)
            , m_dims(dimensions)
        {
        }

        reference operator*() const
        {
            return m_index;
        }

        MultiIndexIterator& operator++()
        {
            if constexpr (MultiIndex::size > 0) {
                increment_index();
            }
            return *this;
        }

        MultiIndexIterator operator++(int)
        {
            auto retval = *this;
            ++(*this);
            return retval;
        }

        bool operator==(MultiIndexIterator& other) const
        {
            return m_index == other.m_index;
        }

        bool operator!=(MultiIndexIterator& other) const
        {
            return !(*this == other);
        }

    private:
        template <size_t I = MultiIndex::size - 1>
        inline void increment_index()
        {
            assert(mio::get<I>(m_dims).get() > 0 && "All dimensions must be positive.");
            assert(mio::get<I>(m_index) < mio::get<I>(m_dims) && "Index out of bounds.");
            if constexpr (I > 0) {
                ++mio::get<I>(m_index);
                if (mio::get<I>(m_index) == mio::get<I>(m_dims)) {
                    mio::get<I>(m_index) = mio::get<I>(m_index).Zero();
                    increment_index<I - 1>();
                }
            }
            else {
                ++mio::get<0>(m_index);
            }
        }
        value_type m_index;
        reference m_dims;
    };

    MultiIndexIterator begin() const
    {
        return MultiIndexIterator(MultiIndex::Zero(), m_dimensions);
    }
    MultiIndexIterator end() const
    {
        return MultiIndexIterator(m_end, m_dimensions);
    }

private:
    MultiIndex make_end_index()
    {
        auto end         = MultiIndex::Zero();
        mio::get<0>(end) = mio::get<0>(m_dimensions);
        return end;
    }
    const MultiIndex m_dimensions; ///< Upper bound for each Index in MultiIndex.
    const MultiIndex m_end; ///< Index returned by
};

} // namespace mio

#endif