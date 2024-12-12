/* 
* Copyright (C) 2020-2024 German Aerospace Center (DLR-SC)
*
* Authors: René Schmieding, Julia Bicker
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

#ifndef MIO_D_ABM_MODEL_H
#define MIO_D_ABM_MODEL_H

namespace mio
{

namespace dabm
{
template <class Implementation>
class Model : public Implementation
{
public:
    using Implementation::Implementation;

    using Implementation::adopt;
    using Implementation::adoption_rate;
    using Implementation::get_rng;
    using Implementation::move;
    using Implementation::time_point;

    using Status = typename Implementation::Status;
    using Agent  = typename Implementation::Agent;

    inline constexpr void check_constraints() const
    {
    }
};

} // namespace dabm
} // namespace mio

#endif
