/* 
* Copyright (C) 2020-2025 MEmilio
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
#ifndef FLOW_H_
#define FLOW_H_

#include <type_traits>

namespace mio
{

/**
 * @brief A Flow defines a possible transition between two Compartments in a FlowModel.
 * Use in a TypeList to define the "Flows" parameter of a FlowModel.
 */
template <auto Source, auto Target>
struct Flow {
    using Type = decltype(Source);
    static_assert(std::is_same<Type, decltype(Target)>::value && std::is_enum<Type>::value,
                  "The Source and Target parameters of a Flow must have the same enum type.");
    const static Type source = Source;
    const static Type target = Target;
};

} // namespace mio

#endif // FLOW_H_