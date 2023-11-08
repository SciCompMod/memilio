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
#ifndef MIO_UTILS_MULTI_INDEX_RANGE_H_
#define MIO_UTILS_MULTI_INDEX_RANGE_H_

#include "memilio/utils/index.h"

namespace mio
{

/**
 * @brief A Range that can be used to iterate over a MultiIndex.
 * @tparam MultiIndex A type mio::Index<Categories...> for some set of categories.
 */
template <class MultiIndex>
class MultiIndexRange
{
public:
    /**
     * @brief Construct a Range that can be used to iterate over all MultiIndices in the given dimensions.
     * The range for each Index i in the MultiIndex is determined by [0, d_i), where d_i is the dimension of i.
     * @param[in] dimensions A MultIndex that contains the dimension for each Category.
     */
    MultiIndexRange(const MultiIndex& dimensions)
        : m_dimensions(dimensions)
        , m_end([](const MultiIndex& dims) {
            // set end to the first invalid index that is reached by increments of 1,
            // i.e. 0 everywhere except for the most significant index (position 0),
            // which is set to its dimension
            auto end         = MultiIndex::Zero();
            mio::get<0>(end) = mio::get<0>(dims);
            return end;
        }(dimensions))
    {
    }

    /**
     * @brief Iterator for MultiIndices.
     */
    class MultiIndexIterator
    {
    public:
        using value_type = MultiIndex;
        using reference  = const MultiIndex&;

        /**
         * @brief Iterator for MultiIndices.
         * @param index Initial value for the iterator position.
         * @param dimensions A reference to the dimensions of the MultiIndex.
         */
        explicit MultiIndexIterator(value_type index, reference dimensions)
            : m_index(index)
            , m_dims(dimensions)
        {
        }

        /// Dereference operator.
        reference operator*() const
        {
            return m_index;
        }

        /// Pre-increment operator.
        MultiIndexIterator& operator++()
        {
            if constexpr (MultiIndex::size > 0) {
                increment_index();
            }
            return *this;
        }

        /// Post-increment operator.
        MultiIndexIterator operator++(int)
        {
            auto retval = *this;
            ++(*this);
            return retval;
        }

        /// Equality operator.
        bool operator==(MultiIndexIterator& other) const
        {
            return m_index == other.m_index;
        }

        /// Inequality operator.
        bool operator!=(MultiIndexIterator& other) const
        {
            return !(*this == other);
        }

    private:
        /**
         * @brief Implementation of ++. Increments m_index.
         * @tparam I position in the MultiIndex.
         */
        template <size_t I = MultiIndex::size - 1>
        inline void increment_index()
        {
            assert(mio::get<I>(m_dims).get() > 0 && "All dimensions must be positive.");
            assert(mio::get<I>(m_index) < mio::get<I>(m_dims) && "Index out of bounds.");

            if constexpr (I > 0) {
                // increment first, then do a carry-over if necessary
                ++mio::get<I>(m_index);
                if (mio::get<I>(m_index) == mio::get<I>(m_dims)) {
                    mio::get<I>(m_index) = mio::get<I>(m_index).Zero();
                    increment_index<I - 1>();
                }
            }
            else {
                ++mio::get<0>(m_index);
                // no carry check for the most significant index
            }
        }
        value_type m_index; ///< Index used for iteration.
        reference m_dims; ///< Read only reference to range dimensions.
    };

    /**
     * @brief STL iterator for a MultiIndexRange.
     * @return Returns the first index for the given dimensions, i.e. 0.
     */
    MultiIndexIterator begin() const
    {
        return MultiIndexIterator(MultiIndex::Zero(), m_dimensions);
    }

    /**
     * @brief STL iterator for a MultiIndexRange.
     * @return Returns the first index outside of the given dimensions.
     */
    MultiIndexIterator end() const
    {
        return MultiIndexIterator(m_end, m_dimensions);
    }

private:
    MultiIndex m_dimensions; ///< Contains strict upper bounds for each Index in MultiIndex.
    MultiIndex m_end; ///< Index to stop at when iterating using increments.
};

} // namespace mio

#endif // MIO_UTILS_MULTI_INDEX_RANGE_H_