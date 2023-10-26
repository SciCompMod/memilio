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
#ifndef MIO_UTILS_TYPE_LIST_H_
#define MIO_UTILS_TYPE_LIST_H_

#include "memilio/io/io.h"
#include "memilio/utils/index.h"

#include <tuple>
#include <type_traits>
#include <utility>

namespace mio
{

/// @brief Collection of types. Each type is mapped to an index of type size_t.
template <class... Types>
struct TypeList {
    /**
     * @brief Get index of a given type.
     * @tparam Type A type contained in TypeList.
     * @return Position of Type within TypeList.
     */
    template <class Type>
    static constexpr size_t index_of()
    {
        return index_of_v<Type, Types...>;
    }

    /// @brief returns the number of Types in TypeList
    static constexpr size_t size()
    {
        return sizeof...(Types);
    }
};

/// Specialization of type_at for TypeList. @see type_at.
template <size_t Index, class... Types>
struct type_at<Index, TypeList<Types...>> : public type_at<Index, Types...> {
};

/// Specialization of index_of for TypeList. @see index_of.
template <class Type, class... Types>
struct index_of<Type, TypeList<Types...>> : public index_of<Type, Types...> {
};

} // namespace mio

#endif // MIO_UTILS_TYPE_LIST_H_