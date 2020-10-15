#include "epidemiology/model/populations.h"
#include <gtest/gtest.h>
#include <array>

enum class InfectionType
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

enum class AgeGroup
{
    TwelveAndYounger,
    TwelveToTwenty,
    TwentyToThirtyfive,
    ThirtyfiveToFourtyFive,
    FourtyFiveToSixtyFive,
    SixtyFiveToSeventyFive,
    SeventyFiveAndOlder,
    Count
};

enum class Continent
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

    int num_compartments =
        static_cast<int>(InfectionType::Count) * static_cast<int>(AgeGroup::Count) * static_cast<int>(Continent::Count);
    ASSERT_EQ(7 * 7 * 8, num_compartments);

    epi::Populations<InfectionType, AgeGroup, Continent> m;

    ASSERT_EQ(num_compartments, m.get_num_compartments());
    ASSERT_EQ(num_compartments, m.get_compartments().size());
    auto category_sizes = m.dimensions;
    ASSERT_EQ(category_sizes[0], static_cast<size_t>(InfectionType::Count));
    ASSERT_EQ(category_sizes[1], static_cast<size_t>(AgeGroup::Count));
    ASSERT_EQ(category_sizes[2], static_cast<size_t>(Continent::Count));
    ASSERT_EQ(0, m.get_total());
}

TEST(TestPopulations, set_population)
{
    epi::Populations<InfectionType, AgeGroup, Continent> m;

    m.set_total(1.);
    int num_compartments =
        static_cast<int>(InfectionType::Count) * static_cast<int>(AgeGroup::Count) * static_cast<int>(Continent::Count);

    for (size_t i = 0; i < static_cast<size_t>(InfectionType::Count); ++i) {
        for (size_t j = 0; j < static_cast<size_t>(AgeGroup::Count); ++j) {
            for (size_t k = 0; k < static_cast<size_t>(Continent::Count); ++k) {
                ASSERT_NEAR(1. / num_compartments,
                            m.get(static_cast<InfectionType>(i), static_cast<AgeGroup>(j), static_cast<Continent>(k)),
                            1e-12);
            }
        }
    }
    ASSERT_NEAR(1., m.get_total(), 1e-12);
}

TEST(TestPopulations, group_population)
{
    epi::Populations<InfectionType, AgeGroup, Continent> m;

    m.set_total(1.);
    int num_compartments =
        static_cast<int>(InfectionType::Count) * static_cast<int>(AgeGroup::Count) * static_cast<int>(Continent::Count);

    ASSERT_NEAR(1. / static_cast<size_t>(AgeGroup::Count), m.get_group_total<AgeGroup>(AgeGroup::FourtyFiveToSixtyFive),
                1e-12);
    m.set_group_total<AgeGroup>(1., AgeGroup::FourtyFiveToSixtyFive);
    ASSERT_NEAR(1., m.get_group_total<AgeGroup>(AgeGroup::FourtyFiveToSixtyFive), 1e-12);
    ASSERT_NEAR(2 - 1. / static_cast<size_t>(AgeGroup::Count), m.get_total(), 1e-12);

    Eigen::VectorXd y = m.get_compartments();
    size_t idx        = 0;
    for (size_t i = 0; i < static_cast<size_t>(InfectionType::Count); ++i) {
        for (size_t j = 0; j < static_cast<size_t>(AgeGroup::Count); ++j) {
            for (size_t k = 0; k < static_cast<size_t>(Continent::Count); ++k) {
                ASSERT_EQ(idx, m.get_flat_index(static_cast<InfectionType>(i), static_cast<AgeGroup>(j),
                                                static_cast<Continent>(k)));

                if (j == static_cast<size_t>(AgeGroup::FourtyFiveToSixtyFive)) {
                    ASSERT_NEAR(y[idx],
                                1. /
                                    (static_cast<size_t>(InfectionType::Count) * static_cast<size_t>(Continent::Count)),
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
    epi::Populations<InfectionType, AgeGroup, Continent> m;

    m.set(100, InfectionType::S, AgeGroup::TwentyToThirtyfive, Continent::Africa);

    m.set_difference_from_total(1000, InfectionType::E, AgeGroup::TwentyToThirtyfive, Continent::Africa);
    ASSERT_NEAR(1000, m.get_total(), 1e-12);
    ASSERT_NEAR(900, m.get(InfectionType::E, AgeGroup::TwentyToThirtyfive, Continent::Africa), 1e-12);

    m.set_difference_from_total(2000, InfectionType::E, AgeGroup::TwentyToThirtyfive, Continent::Africa);
    ASSERT_NEAR(2000, m.get_total(), 1e-12);
    ASSERT_NEAR(1900, m.get(InfectionType::E, AgeGroup::TwentyToThirtyfive, Continent::Africa), 1e-12);

    for (size_t i = 0; i < static_cast<size_t>(InfectionType::Count); ++i) {
        for (size_t j = 0; j < static_cast<size_t>(AgeGroup::Count); ++j) {
            for (size_t k = 0; k < static_cast<size_t>(Continent::Count); ++k) {
                if (i == static_cast<size_t>(InfectionType::S) &&
                    j == static_cast<size_t>(AgeGroup::TwentyToThirtyfive) &&
                    k == static_cast<size_t>(Continent::Africa)) {
                    ASSERT_NEAR(
                        100, m.get(static_cast<InfectionType>(i), static_cast<AgeGroup>(j), static_cast<Continent>(k)),
                        1e-12);
                }
                else if (i == static_cast<size_t>(InfectionType::E) &&
                         j == static_cast<size_t>(AgeGroup::TwentyToThirtyfive) &&
                         k == static_cast<size_t>(Continent::Africa)) {
                    ASSERT_NEAR(
                        1900, m.get(static_cast<InfectionType>(i), static_cast<AgeGroup>(j), static_cast<Continent>(k)),
                        1e-12);
                }
                else {
                    ASSERT_NEAR(
                        0, m.get(static_cast<InfectionType>(i), static_cast<AgeGroup>(j), static_cast<Continent>(k)),
                        1e-12);
                }
            }
        }
    }
}

TEST(TestPopulations, set_difference_from_group_total)
{
    epi::Populations<InfectionType, AgeGroup, Continent> m;

    m.set(100, InfectionType::S, AgeGroup::TwentyToThirtyfive, Continent::Africa);
    m.set(200, InfectionType::S, AgeGroup::TwentyToThirtyfive, Continent::Europe);

    m.set_difference_from_group_total<Continent>(1000, Continent::Africa, InfectionType::E,
                                                 AgeGroup::TwentyToThirtyfive, Continent::Africa);
    ASSERT_NEAR(1000, m.get_group_total<Continent>(Continent::Africa), 1e-12);
    ASSERT_NEAR(900, m.get(InfectionType::E, AgeGroup::TwentyToThirtyfive, Continent::Africa), 1e-12);
    ASSERT_NEAR(1200, m.get_total(), 1e-12);

    m.set_difference_from_group_total<Continent>(2000, Continent::Africa, InfectionType::E,
                                                 AgeGroup::TwentyToThirtyfive, Continent::Africa);
    ASSERT_NEAR(2000, m.get_group_total<Continent>(Continent::Africa), 1e-12);
    ASSERT_NEAR(1900, m.get(InfectionType::E, AgeGroup::TwentyToThirtyfive, Continent::Africa), 1e-12);
    ASSERT_NEAR(2200, m.get_total(), 1e-12);

    for (size_t i = 0; i < static_cast<size_t>(InfectionType::Count); ++i) {
        for (size_t j = 0; j < static_cast<size_t>(AgeGroup::Count); ++j) {
            for (size_t k = 0; k < static_cast<size_t>(Continent::Count); ++k) {
                if (i == static_cast<size_t>(InfectionType::S) &&
                    j == static_cast<size_t>(AgeGroup::TwentyToThirtyfive) &&
                    k == static_cast<size_t>(Continent::Africa)) {
                    ASSERT_NEAR(
                        100, m.get(static_cast<InfectionType>(i), static_cast<AgeGroup>(j), static_cast<Continent>(k)),
                        1e-12);
                }
                else if (i == static_cast<size_t>(InfectionType::E) &&
                         j == static_cast<size_t>(AgeGroup::TwentyToThirtyfive) &&
                         k == static_cast<size_t>(Continent::Africa)) {
                    ASSERT_NEAR(
                        1900, m.get(static_cast<InfectionType>(i), static_cast<AgeGroup>(j), static_cast<Continent>(k)),
                        1e-12);
                }
                else if (i == static_cast<size_t>(InfectionType::S) &&
                         j == static_cast<size_t>(AgeGroup::TwentyToThirtyfive) &&
                         k == static_cast<size_t>(Continent::Europe)) {
                    ASSERT_NEAR(
                        200, m.get(static_cast<InfectionType>(i), static_cast<AgeGroup>(j), static_cast<Continent>(k)),
                        1e-12);
                }
                else {
                    ASSERT_NEAR(
                        0, m.get(static_cast<InfectionType>(i), static_cast<AgeGroup>(j), static_cast<Continent>(k)),
                        1e-12);
                }
            }
        }
    }
}
