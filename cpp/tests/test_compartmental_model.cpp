#include "epidemiology/compartmental_model.h"
#include <gtest/gtest.h>

TEST(TestCompartmentalModel, set_population)
{
    enum InfectionType
    {
        S,
        E,
        C,
        I,
        H,
        U,
        R,
        D,
        InfectionTypeCount
    };

    enum AgeGroup
    {
        TwelveAndYounger,
        TwelveToTwenty,
        TwentyToThirtyfive,
        ThirtyfiveToFourtyFive,
        FourtyFiveToSixtyFive,
        SixtyFiveToSeventyFive,
        SeventyFiveAndOlder,
        AgeGroupCount
    };

    enum Continent
    {
        Europe,
        Asia,
        NorthAmerica,
        SouthAmerica,
        Australia,
        Antarctica,
        Africa,
        ContinentCount
    };

    int num_compartments = InfectionTypeCount * AgeGroupCount * ContinentCount;
    ASSERT_EQ(7 * 7 * 8, num_compartments);

    epi::CompartmentalModel m({InfectionTypeCount, AgeGroupCount, ContinentCount});

    ASSERT_EQ(num_compartments, m.get_num_compartments());
    ASSERT_EQ(num_compartments, m.get_compartments().size());
    auto category_sizes = m.get_category_sizes();
    ASSERT_EQ(category_sizes[0], InfectionTypeCount);
    ASSERT_EQ(category_sizes[1], AgeGroupCount);
    ASSERT_EQ(category_sizes[2], ContinentCount);
    ASSERT_EQ(0, m.get_total_population());

    m.set_total_populaton(1.);
    for (size_t i = 0; i < InfectionTypeCount; ++i) {
        for (size_t j = 0; j < AgeGroupCount; ++j) {
            for (size_t k = 0; k < ContinentCount; ++k) {
                ASSERT_NEAR(1. / num_compartments, m.get({i, j, k}), 1e-12);
            }
        }
    }
    ASSERT_NEAR(1., m.get_total_population(), 1e-12);
}
