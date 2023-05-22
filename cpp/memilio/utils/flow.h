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
#ifndef FLOW_H_
#define FLOW_H_

namespace mio
{

/// @brief A Flow defines a transition between two Compartments in a CompartmentalModel. Use with TypeChart
template <class Status, Status Source, Status Target>
struct Flow {
    using Type                 = Status;
    const static Status source = Source;
    const static Status target = Target;
};

} // namespace mio

#endif // FLOW_H_