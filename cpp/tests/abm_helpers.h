/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, David Kerkmann, Sascha Korf, Martin J. Kuehn, Khoa Nguyen
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
#ifndef ABM_HELPERS_H
#define ABM_HELPERS_H

#include "abm/model.h"
#include "abm/person_id.h"

#include "gmock/gmock.h"

// Assign the name to general age group.
const size_t num_age_groups   = 6;
const auto age_group_0_to_4   = mio::AgeGroup(0);
const auto age_group_5_to_14  = mio::AgeGroup(1);
const auto age_group_15_to_34 = mio::AgeGroup(2);
const auto age_group_35_to_59 = mio::AgeGroup(3);
const auto age_group_60_to_79 = mio::AgeGroup(4);
const auto age_group_80_plus  = mio::AgeGroup(5);

/**
 * mock of the generator function of DistributionAdapter<DistT>.
 * can't be used directly as a generator function because it is not copyable.
 * see MockDistributionRef
 */
template <class DistT>
struct MockDistribution {
    using Distribution = DistT;
    // using invoke() instead of operator() because operators cant be mocked in the GMock framework */
    MOCK_METHOD(typename Distribution::ResultType, invoke, (const typename Distribution::ParamType&), ());
};

/**
 * reference wrapper of a MockDistribution object.
 * Mocks are not copyable but the generator function of a distribution must be copyable.
 * This wrapper is copyable and all copies redirect invocations to a shared underlying mock
 * so it can be used as a generator function.
 */
template <class MockDistribution>
struct MockDistributionRef {
    using Distribution = typename MockDistribution::Distribution;
    typename Distribution::ResultType operator()(const typename Distribution::ParamType& p)
    {
        return mock->invoke(p);
    }
    std::shared_ptr<MockDistribution> mock = std::make_shared<MockDistribution>();
};

/**
 * Replaces the generator function in the static instance of DistributionAdapter with a mock.
 * On construction sets the generator and on destruction restores the previous generator.
 */
template <class MockDistribution>
struct ScopedMockDistribution {
    using Distribution = typename MockDistribution::Distribution;
    /**
     * constructor replaces the generator function with a mock
     */
    ScopedMockDistribution()
    {
        old = Distribution::get_instance().get_generator();
        Distribution::get_instance().set_generator(mock_ref);
    }
    ~ScopedMockDistribution()
    {
        Distribution::get_instance().set_generator(old);
    }
    MockDistribution& get_mock()
    {
        return *mock_ref.mock;
    }

    MockDistributionRef<MockDistribution> mock_ref;
    typename Distribution::GeneratorFunction old;
};

/**
 * @brief Create a Person without a Model object. Intended for simple use in tests.
 */
mio::abm::Person make_test_person(mio::RandomNumberGenerator& rng, mio::abm::Location& location,
                                  mio::AgeGroup age                        = age_group_15_to_34,
                                  mio::abm::InfectionState infection_state = mio::abm::InfectionState::Susceptible,
                                  mio::abm::TimePoint t                    = mio::abm::TimePoint(0),
                                  mio::abm::Parameters params              = mio::abm::Parameters(num_age_groups),
                                  mio::abm::PersonId id                    = mio::abm::PersonId(0));

/**
 * @brief Add a Person to the Model. Intended for simple use in tests.
 */
mio::abm::PersonId add_test_person(mio::abm::Model& model, mio::abm::LocationId loc_id,
                                   mio::AgeGroup age                        = age_group_15_to_34,
                                   mio::abm::InfectionState infection_state = mio::abm::InfectionState::Susceptible,
                                   mio::abm::TimePoint t                    = mio::abm::TimePoint(0));

/// @brief Calls mio::abm::interact, but it computes the correct exposures for you.
void interact_testing(mio::abm::PersonalRandomNumberGenerator& personal_rng, mio::abm::Person& person,
                      const mio::abm::Location& location, const std::vector<mio::abm::Person>& local_population,
                      const mio::abm::TimePoint t, const mio::abm::TimeSpan dt,
                      const mio::abm::Parameters& global_parameters);

#endif //ABM_HELPERS_H
