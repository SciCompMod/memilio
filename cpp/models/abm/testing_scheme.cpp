/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Elisabeth Kluth
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
#include "abm/testing_scheme.h"
#include "abm/world.h"
#include "abm/location.h"
#include "abm/parameters.h"
#include "memilio/utils/random_number_generator.h"

namespace mio
{
namespace abm
{

TestingScheme::TestingScheme(TimeSpan interval, double probability)
    : m_time_interval(interval)
    , m_probability(probability)
{
}

TestingScheme::TestingScheme()
    : TestingScheme(seconds(std::numeric_limits<int>::max()), 1)
{
}

bool TestingScheme::run_scheme(Person& person, const GlobalTestingParameters& params) const
{
    if (person.get_time_since_negative_test() > m_time_interval) {
        double random = UniformDistribution<double>::get_instance()();
        if (random < m_probability) {
            return !person.get_tested(params.get<AntigenTest>());
        }
    }
    return true;
}

} // namespace abm
} // namespace mio
