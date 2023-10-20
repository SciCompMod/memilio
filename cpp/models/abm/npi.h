/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: David Kerkmann
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
#ifndef EPI_ABM_NPI_H
#define EPI_ABM_NPI_H

#include <functional>
#include "abm/mask_type.h"
#include "abm/time.h"
#include "abm/testing_strategy.h"
#include "memilio/utils/random_number_generator.h"

namespace mio
{
namespace abm
{

/**
 * TODO: Template everything with variable arguments.
 * TODO: Probably remove std::function for better efficiency and access to final NPI data members.
*/
class Location;
class Person;
using LocationNPIStrategy =
    std::function<bool(Person::RandomNumberGenerator&, const Person&, const Location&, TimePoint)>;

/**
 * @brief Base NPI to always allow entrance.
*/
struct AcceptBase {
    virtual ~AcceptBase() = default;
    virtual bool operator()(Person::RandomNumberGenerator&, const Person&, const Location&, TimePoint)
    {
        return true;
    }
};

} // namespace abm
} // namespace mio

#endif //EPI_ABM_NPI_H
