#include "epidemiology/model/populations.h"
#include <gtest/gtest.h>
#include <array>

// Three categories, one defined by an enum, one by an enum class and one by a struct.
enum class InfectionState
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

struct AgeGroup{};

enum Continent
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

    epi::Index<InfectionState> num_infType(InfectionState::Count);
    epi::Index<AgeGroup>      num_ageGroup(7);
    epi::Index<Continent>     num_continents(Continent::Count);

    size_t num_compartments = (size_t)num_infType * (size_t)num_ageGroup * (size_t)num_continents;
    ASSERT_EQ(7 * 7 * 8, num_compartments);

    epi::Populations<InfectionState, AgeGroup, Continent> m({num_infType, num_ageGroup, num_continents});

    ASSERT_EQ(num_compartments, m.get_num_compartments());
    ASSERT_EQ(num_compartments, m.numel());
    ASSERT_EQ(num_compartments, (size_t)m.get_compartments().size());

    ASSERT_EQ(m.size<InfectionState>(), num_infType);
    ASSERT_EQ(m.size<AgeGroup>(), num_ageGroup);
    ASSERT_EQ(m.size<Continent>(), num_continents);
    ASSERT_EQ(0, m.get_total());
}

TEST(TestPopulations, set_population)
{
    epi::Index<InfectionState> num_infType(InfectionState::Count);
    epi::Index<AgeGroup>      num_ageGroup(7);
    epi::Index<Continent>     num_continents(Continent::Count);

    epi::Populations<InfectionState, AgeGroup, Continent> m({num_infType, num_ageGroup, num_continents});

    m.set_total(1.);
    size_t num_compartments = (size_t)num_infType * (size_t)num_ageGroup * (size_t)num_continents;

    for (auto i = epi::Index<InfectionState>(0); i < m.size<InfectionState>(); ++i) {
        for (auto j = epi::Index<AgeGroup>(0); j < m.size<AgeGroup>(); ++j) {
            for (auto k = epi::Index<Continent>(0); k < m.size<Continent>(); ++k) {
                ASSERT_NEAR(1. / num_compartments, (m[{i, j ,k}]), 1e-12);
            }
        }
    }
    ASSERT_NEAR(1., m.get_total(), 1e-12);
}

TEST(TestPopulations, group_population)
{
    epi::Index<InfectionState> num_infType(InfectionState::Count);
    epi::Index<AgeGroup>      num_ageGroup(7);
    epi::Index<Continent>     num_continents(Continent::Count);

    epi::Populations<InfectionState, AgeGroup, Continent> m({num_infType, num_ageGroup, num_continents});

    m.set_total(1.);
    size_t num_compartments = (size_t)num_infType * (size_t)num_ageGroup * (size_t)num_continents;

    epi::Index<AgeGroup> fortyToFifty(5);
    ASSERT_NEAR(1. / static_cast<size_t>(num_ageGroup), m.get_group_total(epi::Index<AgeGroup>(5)), 1e-12);
    m.set_group_total(epi::Index<AgeGroup>(5), 1.);
    ASSERT_NEAR(1., m.get_group_total(epi::Index<AgeGroup>(5)), 1e-12);
    ASSERT_NEAR(2 - 1. / static_cast<size_t>(num_ageGroup), m.get_total(), 1e-12);

    Eigen::VectorXd y = m.get_compartments();
    size_t idx        = 0;
    for (auto i = epi::Index<InfectionState>(0); i < m.size<InfectionState>(); ++i) {
        for (auto j = epi::Index<AgeGroup>(0); j < m.size<AgeGroup>(); ++j) {
            for (auto k = epi::Index<Continent>(0); k < m.size<Continent>(); ++k) {
                ASSERT_EQ(idx, m.get_flat_index({i, j , k}));

                if (j == fortyToFifty) {
                    ASSERT_NEAR(y[idx],
                                1. /
                                    (static_cast<size_t>(num_infType) * static_cast<size_t>(num_continents)),
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
    epi::Index<InfectionState> num_infType(InfectionState::Count);
    epi::Index<AgeGroup>      num_ageGroup(7);
    epi::Index<Continent>     num_continents(Continent::Count);

    using Po = epi::Populations<InfectionState, AgeGroup, Continent>;
    Po m({num_infType, num_ageGroup, num_continents});


    Po::MultiIndex S_2_Africa = {epi::Index<InfectionState>(InfectionState::S),
                                 epi::Index<AgeGroup>(2),
                                 epi::Index<Continent>(Africa)};

    Po::MultiIndex E_2_Africa = {epi::Index<InfectionState>(InfectionState::E),
                                 epi::Index<AgeGroup>(2),
                                 epi::Index<Continent>(Africa)};

    m[S_2_Africa] = 100;

    m.set_difference_from_total(E_2_Africa, 1000);
    ASSERT_NEAR(1000, m.get_total(), 1e-12);
    ASSERT_NEAR(900, (m[E_2_Africa]), 1e-12);

    m.set_difference_from_total(E_2_Africa, 2000);
    ASSERT_NEAR(2000, m.get_total(), 1e-12);
    ASSERT_NEAR(1900, (m[E_2_Africa]), 1e-12);

    for (auto i = epi::Index<InfectionState>(0); i < m.size<InfectionState>(); ++i) {
        for (auto j = epi::Index<AgeGroup>(0); j < m.size<AgeGroup>(); ++j) {
            for (auto k = epi::Index<Continent>(0); k < m.size<Continent>(); ++k) {
                auto current = Po::MultiIndex(i, j, k);
                if ( current == S_2_Africa ) {
                    ASSERT_NEAR(100, (m[current]), 1e-12);
                }
                else if ( current == E_2_Africa ) {
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
    epi::Index<InfectionState> num_infType(InfectionState::Count);
    epi::Index<AgeGroup>      num_ageGroup(7);
    epi::Index<Continent>     num_continents(Continent::Count);

    using Po = epi::Populations<InfectionState, AgeGroup, Continent>;
    Po m({num_infType, num_ageGroup, num_continents});

    Po::MultiIndex S_2_Africa = {epi::Index<InfectionState>(InfectionState::S),
                                 epi::Index<AgeGroup>(2),
                                 epi::Index<Continent>(Africa)};

    Po::MultiIndex E_2_Africa = {epi::Index<InfectionState>(InfectionState::E),
                                 epi::Index<AgeGroup>(2),
                                 epi::Index<Continent>(Africa)};

    Po::MultiIndex S_2_Europe = {epi::Index<InfectionState>(InfectionState::E),
                                 epi::Index<AgeGroup>(2),
                                 epi::Index<Continent>(Europe)};

    m[S_2_Africa] = 100;
    m[S_2_Europe] = 200;

    m.set_difference_from_group_total<Continent>(E_2_Africa, 1000);
    ASSERT_NEAR(1000, m.get_group_total(epi::Index<Continent>(Africa)), 1e-12);
    ASSERT_NEAR(900, (m[E_2_Africa]), 1e-12);
    ASSERT_NEAR(1200, m.get_total(), 1e-12);

    m.set_difference_from_group_total<Continent>(E_2_Africa, 2000);
    ASSERT_NEAR(2000, m.get_group_total(epi::Index<Continent>(Africa)), 1e-12);
    ASSERT_NEAR(1900, (m[E_2_Africa]), 1e-12);
    ASSERT_NEAR(2200, m.get_total(), 1e-12);
    for (auto i = epi::Index<InfectionState>(0); i < m.size<InfectionState>(); ++i) {
        for (auto j = epi::Index<AgeGroup>(0); j < m.size<AgeGroup>(); ++j) {
            for (auto k = epi::Index<Continent>(0); k < m.size<Continent>(); ++k) {
                auto current = Po::MultiIndex(i, j, k);
                if (current  == S_2_Africa) {
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
