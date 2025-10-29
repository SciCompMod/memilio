/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Jan Kleinert
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
#include "memilio/epidemiology/populations.h"
#include "memilio/utils/logging.h"
#include <gtest/gtest.h>
#include <array>

// Three categories, one defined by an enum, one by an enum class and one by a struct.
enum class TestInfectionState
{
    S,
    E,
    C,
    I,
    H,
    U,
    R,
    D,
    Count
};

struct TestAgeGroup {
};

enum TestContinent
{
    Europe,
    Asia,
    NorthAmerica,
    SouthAmerica,
    Australia,
    Antarctica,
    Africa,
    Count
};

TEST(TestPopulations, sizes)
{

    mio::Index<TestInfectionState> num_infType(TestInfectionState::Count);
    mio::Index<TestAgeGroup> num_ageGroup(7);
    mio::Index<TestContinent> num_continents(TestContinent::Count);

    size_t num_compartments = (size_t)num_infType * (size_t)num_ageGroup * (size_t)num_continents;
    ASSERT_EQ(7 * 7 * 8, num_compartments);

    mio::Populations<double, TestInfectionState, TestAgeGroup, TestContinent> m(
        {num_infType, num_ageGroup, num_continents});

    ASSERT_EQ(num_compartments, m.get_num_compartments());
    ASSERT_EQ(num_compartments, m.numel());
    ASSERT_EQ(num_compartments, (size_t)m.get_compartments().size());

    ASSERT_EQ(m.size<TestInfectionState>(), num_infType);
    ASSERT_EQ(m.size<TestAgeGroup>(), num_ageGroup);
    ASSERT_EQ(m.size<TestContinent>(), num_continents);
    ASSERT_EQ(0, m.get_total());
}

TEST(TestPopulations, set_population)
{
    mio::Index<TestInfectionState> num_infType(TestInfectionState::Count);
    mio::Index<TestAgeGroup> num_ageGroup(7);
    mio::Index<TestContinent> num_continents(TestContinent::Count);

    mio::Populations<double, TestInfectionState, TestAgeGroup, TestContinent> m(
        {num_infType, num_ageGroup, num_continents});

    m.set_total(1.);
    size_t num_compartments = (size_t)num_infType * (size_t)num_ageGroup * (size_t)num_continents;

    for (auto i = mio::Index<TestInfectionState>(0); i < m.size<TestInfectionState>(); ++i) {
        for (auto j = mio::Index<TestAgeGroup>(0); j < m.size<TestAgeGroup>(); ++j) {
            for (auto k = mio::Index<TestContinent>(0); k < m.size<TestContinent>(); ++k) {
                ASSERT_NEAR(1. / num_compartments, (m[{i, j, k}]), 1e-12);
            }
        }
    }
    ASSERT_NEAR(1., m.get_total(), 1e-12);
}

TEST(TestPopulations, group_population)
{
    mio::Index<TestInfectionState> num_infType(TestInfectionState::Count);
    mio::Index<TestAgeGroup> num_ageGroup(7);
    mio::Index<TestContinent> num_continents(TestContinent::Count);

    mio::Populations<double, TestInfectionState, TestAgeGroup, TestContinent> m(
        {num_infType, num_ageGroup, num_continents});

    m.set_total(1.);
    size_t num_compartments = (size_t)num_infType * (size_t)num_ageGroup * (size_t)num_continents;

    mio::Index<TestAgeGroup> fortyToFifty(5);
    ASSERT_NEAR(1. / static_cast<size_t>(num_ageGroup), m.get_group_total(mio::Index<TestAgeGroup>(5)), 1e-12);
    m.set_group_total(mio::Index<TestAgeGroup>(5), 1.);
    ASSERT_NEAR(1., m.get_group_total(mio::Index<TestAgeGroup>(5)), 1e-12);
    ASSERT_NEAR(2 - 1. / static_cast<size_t>(num_ageGroup), m.get_total(), 1e-12);

    Eigen::Matrix<double, Eigen::Dynamic, 1> y_tmp = m.get_compartments();
    Eigen::VectorXd y                              = y_tmp.template cast<double>();
    size_t idx                                     = 0;
    for (auto i = mio::Index<TestInfectionState>(0); i < m.size<TestInfectionState>(); ++i) {
        for (auto j = mio::Index<TestAgeGroup>(0); j < m.size<TestAgeGroup>(); ++j) {
            for (auto k = mio::Index<TestContinent>(0); k < m.size<TestContinent>(); ++k) {
                ASSERT_EQ(idx, m.get_flat_index({i, j, k}));

                if (j == fortyToFifty) {
                    ASSERT_NEAR(y[idx], 1. / (static_cast<size_t>(num_infType) * static_cast<size_t>(num_continents)),
                                1e-12);
                }
                else {
                    ASSERT_NEAR(y[idx], 1. / num_compartments, 1e-12);
                }
                idx++;
            }
        }
    }
}

TEST(TestPopulations, set_difference_from_total)
{
    mio::Index<TestInfectionState> num_infType(TestInfectionState::Count);
    mio::Index<TestAgeGroup> num_ageGroup(7);
    mio::Index<TestContinent> num_continents(TestContinent::Count);

    using Po = mio::Populations<double, TestInfectionState, TestAgeGroup, TestContinent>;
    Po m({num_infType, num_ageGroup, num_continents});

    Po::Index S_2_Africa = {mio::Index<TestInfectionState>(TestInfectionState::S), mio::Index<TestAgeGroup>(2),
                            mio::Index<TestContinent>(Africa)};

    Po::Index E_2_Africa = {mio::Index<TestInfectionState>(TestInfectionState::E), mio::Index<TestAgeGroup>(2),
                            mio::Index<TestContinent>(Africa)};

    m[S_2_Africa] = 100;

    m.set_difference_from_total(E_2_Africa, 1000);
    ASSERT_NEAR(1000, m.get_total(), 1e-12);
    ASSERT_NEAR(900, (m[E_2_Africa]), 1e-12);

    m.set_difference_from_total(E_2_Africa, 2000);
    ASSERT_NEAR(2000, m.get_total(), 1e-12);
    ASSERT_NEAR(1900, (m[E_2_Africa]), 1e-12);

    for (auto i = mio::Index<TestInfectionState>(0); i < m.size<TestInfectionState>(); ++i) {
        for (auto j = mio::Index<TestAgeGroup>(0); j < m.size<TestAgeGroup>(); ++j) {
            for (auto k = mio::Index<TestContinent>(0); k < m.size<TestContinent>(); ++k) {
                auto current = Po::Index(i, j, k);
                if (current == S_2_Africa) {
                    ASSERT_NEAR(100, (m[current]), 1e-12);
                }
                else if (current == E_2_Africa) {
                    ASSERT_NEAR(1900, (m[current]), 1e-12);
                }
                else {
                    ASSERT_NEAR(0, (m[current]), 1e-12);
                }
            }
        }
    }
}

TEST(TestPopulations, set_difference_from_group_total)
{
    mio::Index<TestInfectionState> num_infType(TestInfectionState::Count);
    mio::Index<TestAgeGroup> num_ageGroup(7);
    mio::Index<TestContinent> num_continents(TestContinent::Count);

    using Po = mio::Populations<double, TestInfectionState, TestAgeGroup, TestContinent>;
    Po m({num_infType, num_ageGroup, num_continents});

    Po::Index S_2_Africa = {mio::Index<TestInfectionState>(TestInfectionState::S), mio::Index<TestAgeGroup>(2),
                            mio::Index<TestContinent>(Africa)};

    Po::Index E_2_Africa = {mio::Index<TestInfectionState>(TestInfectionState::E), mio::Index<TestAgeGroup>(2),
                            mio::Index<TestContinent>(Africa)};

    Po::Index S_2_Europe = {mio::Index<TestInfectionState>(TestInfectionState::E), mio::Index<TestAgeGroup>(2),
                            mio::Index<TestContinent>(Europe)};

    m[S_2_Africa] = 100;
    m[S_2_Europe] = 200;

    m.set_difference_from_group_total<TestContinent>(E_2_Africa, 1000);
    ASSERT_NEAR(1000, m.get_group_total(mio::Index<TestContinent>(Africa)), 1e-12);
    ASSERT_NEAR(900, (m[E_2_Africa]), 1e-12);
    ASSERT_NEAR(1200, m.get_total(), 1e-12);

    m.set_difference_from_group_total<TestContinent>(E_2_Africa, 2000);
    ASSERT_NEAR(2000, m.get_group_total(mio::Index<TestContinent>(Africa)), 1e-12);
    ASSERT_NEAR(1900, (m[E_2_Africa]), 1e-12);
    ASSERT_NEAR(2200, m.get_total(), 1e-12);
    for (auto i = mio::Index<TestInfectionState>(0); i < m.size<TestInfectionState>(); ++i) {
        for (auto j = mio::Index<TestAgeGroup>(0); j < m.size<TestAgeGroup>(); ++j) {
            for (auto k = mio::Index<TestContinent>(0); k < m.size<TestContinent>(); ++k) {
                auto current = Po::Index(i, j, k);
                if (current == S_2_Africa) {
                    ASSERT_NEAR(100, (m[current]), 1e-12);
                }
                else if (current == E_2_Africa) {
                    ASSERT_NEAR(1900, (m[current]), 1e-12);
                }
                else if (current == S_2_Europe) {
                    ASSERT_NEAR(200, (m[current]), 1e-12);
                }
                else {
                    ASSERT_NEAR(0, (m[current]), 1e-12);
                }
            }
        }
    }
}

TEST(TestPopulations, populations_constraints)
{
    // Verify that the functions check_constraints() and apply_constraints() work as expected.
    mio::Populations<double, TestInfectionState, TestAgeGroup> pop(
        {mio::Index<TestInfectionState>(TestInfectionState::Count), mio::Index<TestAgeGroup>(7)});

    // Check that with valid inputs, the output of both functions is false and
    // that apply_constraints does not change anything.
    pop.set_group_total<TestAgeGroup>(mio::Index<TestAgeGroup>(5), 10);
    EXPECT_FALSE(pop.check_constraints());
    EXPECT_FALSE(pop.apply_constraints());
    EXPECT_NEAR(10, pop.get_total(), 1e-12);

    // Deactivate temporarily log output for next tests.
    mio::set_log_level(mio::LogLevel::off);

    // Check that with invalid inputs, the output of both functions is true and
    // that apply_constraints corrects all wrong values.
    // Set two entries to negative values.
    pop[{mio::Index<TestInfectionState>(TestInfectionState::E), mio::Index<TestAgeGroup>(5)}] = -100;
    pop[{mio::Index<TestInfectionState>(TestInfectionState::H), mio::Index<TestAgeGroup>(3)}] = -10;
    EXPECT_TRUE(pop.check_constraints());
    EXPECT_TRUE(pop.apply_constraints());
    // Negative values should be corrected to zero.
    EXPECT_NEAR(0., (pop[{mio::Index<TestInfectionState>(TestInfectionState::E), mio::Index<TestAgeGroup>(5)}]), 1e-12);
    EXPECT_NEAR(0., (pop[{mio::Index<TestInfectionState>(TestInfectionState::H), mio::Index<TestAgeGroup>(3)}]), 1e-12);

    // Reactive log output.
    mio::set_log_level(mio::LogLevel::warn);
}
