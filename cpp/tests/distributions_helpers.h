/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele
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
#ifndef DISTRIBUTIONS_HELPERS_H
#define DISTRIBUTIONS_HELPERS_H

#include "memilio/utils/parameter_distributions.h"
#include <gtest/gtest.h>
#include <gmock/gmock.h>

void check_distribution(const mio::ParameterDistribution& dist, const mio::ParameterDistribution& dist_read);

//Mocks are not copyable because they need to store the call counters etc.
//ParameterDistribution must be copyable (for clone)
//so our mock consists of two classes
//the first class contains the actual mock logic but does not derive from anything
class MockParameterDistribution
{
public:
    MOCK_METHOD(double, get_rand_sample, (), ());
};
//the second class is clonable etc. and forwards calls to a stable instance of the first class
//templated so it allows StrictMock, NiceMock, etc.
//inherits from ParameterDistributionNormal instead of the base class so it is visitable
//all copies of an instance of this class will share the same mock and call counters etc.
template <class Mock = MockParameterDistribution>
class MockParameterDistributionRef : public mio::ParameterDistributionNormal
{
public:
    using mio::ParameterDistributionNormal::ParameterDistributionNormal;

    double get_rand_sample(mio::RandomNumberGenerator& /*rng*/) override
    {
        return mock->get_rand_sample();
    }

    mio::ParameterDistribution* clone() const override
    {
        return new MockParameterDistributionRef(*this);
    }

    Mock& get_mock()
    {
        return *mock;
    }

private:
    std::shared_ptr<Mock> mock = std::make_shared<Mock>();
};

#endif
