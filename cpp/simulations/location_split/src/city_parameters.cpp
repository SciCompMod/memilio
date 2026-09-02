#include "city_parameters.h"

#include <algorithm>
#include <cmath>
#include <numeric>

namespace CityParameters
{

std::string region_to_string(Region region)
{
    switch (region) {
    case Region::Germany:
        return "germany";
    case Region::France:
        return "france";
    case Region::UnitedStates:
        return "usa";
    }
    return "unknown";
}

bool region_from_string(const std::string& name, Region& region)
{
    if (name == "germany" || name == "ger" || name == "de") {
        region = Region::Germany;
        return true;
    }
    if (name == "france" || name == "fra" || name == "fr") {
        region = Region::France;
        return true;
    }
    if (name == "usa" || name == "us") {
        region = Region::UnitedStates;
        return true;
    }
    return false;
}

/*
 * Sources:
 * - Destatis 2024: Population: Germany, reference date, age. Code 12411-0005.
 * - Destatis 2024: Household statistics. Code 12421-0100.
 * - Destatis 2024: Erwerbstaetigenquoten, age group 15 to under 65.
 * - Destatis: Schueler: Deutschland, Schuljahr, Geschlecht, Schulart, Jahrgangsstufen. Code 21111-0002.
 * - HDE / Statista 2023: number of retail stores in Germany.
 */
RegionParameters germany()
{
    return RegionParameters{Region::Germany,
                            {
                                0.044, // 0-4 years
                                0.094, // 5-14 years
                                0.222, // 15-34 years
                                0.334, // 35-59 years
                                0.233, // 60-79 years
                                0.073 // 80+ years
                            },
                            {
                                0.415, // 1-person households
                                0.342, // 2-person households
                                0.118, // 3-person households
                                0.091, // 4-person households
                                0.034 // 5+ person households
                            },
                            InfrastructureRatios{
                                /*employment_rate=*/0.775,
                                /*people_per_workplace=*/10.0,
                                /*school_rate=*/0.105,
                                /*max_students_per_elementary_school=*/200.0,
                                /*max_students_per_secondary_school=*/300.0,
                                /*ratio_elementary_to_secondary_school=*/1.65,
                                /*retail_stores_per_1000_people=*/3.7,
                                /*people_per_event=*/15.0,
                            }};
}

/*
 * Sources:
 * - INSEE 2024: population by age.
 * - INSEE 2020: Composition des menages.
 * - INSEE 2023: Taux d'emploi, ages 15-64.
 * - Ministere de l'Education nationale: enrolment statistics.
 * - FCD (Federation du Commerce et de la Distribution) 2023.
 */
RegionParameters france()
{
    return RegionParameters{Region::France,
                            {
                                0.052, // 0-4 years
                                0.122, // 5-14 years
                                0.235, // 15-34 years
                                0.319, // 35-59 years
                                0.211, // 60-79 years
                                0.061 // 80+ years
                            },
                            {
                                0.384, // 1-person households
                                0.324, // 2-person households
                                0.130, // 3-person households
                                0.108, // 4-person households
                                0.054 // 5+ person households
                            },
                            InfrastructureRatios{
                                /*employment_rate=*/0.744,
                                /*people_per_workplace=*/10.0,
                                /*school_rate=*/0.151,
                                /*max_students_per_elementary_school=*/200.0,
                                /*max_students_per_secondary_school=*/300.0,
                                /*ratio_elementary_to_secondary_school=*/1.65,
                                /*retail_stores_per_1000_people=*/4.2,
                                /*people_per_event=*/15.0,
                            }};
}

/*
 * Sources:
 * - U.S. Census Bureau: Annual Estimates of the Resident Population 2024.
 * - U.S. Census Bureau: Household Size and Type.
 * - Bureau of Labor Statistics: Employment Situation 2024 (ages 16-64).
 * - National Center for Education Statistics 2023.
 * - U.S. Census Bureau: Annual Retail Trade Survey 2022.
 */
RegionParameters united_states()
{
    return RegionParameters{Region::UnitedStates,
                            {
                                0.057, // 0-4 years
                                0.126, // 5-14 years
                                0.276, // 15-34 years
                                0.318, // 35-59 years
                                0.191, // 60-79 years
                                0.032 // 80+ years
                            },
                            {
                                0.288, // 1-person households
                                0.341, // 2-person households
                                0.152, // 3-person households
                                0.164, // 4-person households
                                0.055 // 5+ person households
                            },
                            InfrastructureRatios{
                                /*employment_rate=*/0.594,
                                /*people_per_workplace=*/10.0,
                                /*school_rate=*/0.128,
                                /*max_students_per_elementary_school=*/200.0,
                                /*max_students_per_secondary_school=*/300.0,
                                /*ratio_elementary_to_secondary_school=*/1.65,
                                /*retail_stores_per_1000_people=*/3.3,
                                /*people_per_event=*/15.0,
                            }};
}

RegionParameters get_region_parameters(Region region)
{
    switch (region) {
    case Region::France:
        return france();
    case Region::UnitedStates:
        return united_states();
    case Region::Germany:
        break;
    }
    return germany();
}

std::vector<int> age_vector(int population, const RegionParameters& params)
{
    std::vector<int> ages(params.age_distribution.size());
    for (size_t i = 0; i < params.age_distribution.size(); ++i) {
        ages[i] = static_cast<int>(population * params.age_distribution[i]);
    }
    return ages;
}

std::vector<int> household_sizes(int population, const RegionParameters& params)
{
    const int n_bins = static_cast<int>(params.household_size_distribution.size());

    // 1. Estimate the household count from the average household size.
    double avg_household_size = 0.0;
    for (int i = 0; i < n_bins; ++i) {
        avg_household_size += (i + 1) * params.household_size_distribution[i];
    }
    const int estimated_total_households = static_cast<int>(std::round(population / avg_household_size));

    // 2. Distribute the households over the size bins and round.
    std::vector<int> households_by_size(n_bins);
    for (int i = 0; i < n_bins; ++i) {
        households_by_size[i] =
            static_cast<int>(std::round(estimated_total_households * params.household_size_distribution[i]));
    }

    // 3. Fix the rounding error by adding or removing single and five person households.
    int total_people = 0;
    for (int i = 0; i < n_bins; ++i) {
        total_people += households_by_size[i] * (i + 1);
    }
    int diff = population - total_people;
    while (diff != 0) {
        const int idx = (diff > 0) ? n_bins - 1 : 0;
        households_by_size[idx] += (diff > 0) ? 1 : -1;
        diff += (diff > 0) ? -(idx + 1) : (idx + 1);
    }

    return households_by_size;
}

namespace
{

std::pair<int, int> calc_num_workplaces_and_worker(int population, const RegionParameters& params)
{
    const auto ages = age_vector(population, params);

    // NOTE: (4 / 20) is integer division and evaluates to 0, so the 60-79 age group does not
    // contribute any potential workers. This is kept as it was on the panvXabmSim branch so that
    // results stay comparable; the shortfall is corrected in CityBuilder, which assigns the
    // remaining workplaces to the 60+ age groups.
    const int n_potential_worker = ages[2] + ages[3] + ((4 / 20) * ages[4]);
    const int n_worker = static_cast<int>(std::round(n_potential_worker * params.ratios.employment_rate));
    const int n_workplaces = static_cast<int>(std::round(n_worker / params.ratios.people_per_workplace));
    return {n_workplaces, n_worker};
}

std::tuple<int, int, int, int> calc_num_elem_and_sec_schools(int population, const RegionParameters& params)
{
    const int total_students = static_cast<int>(population * params.ratios.school_rate);
    const int num_elementary_students =
        static_cast<int>(std::round(total_students / (1.0 + params.ratios.ratio_elementary_to_secondary_school)));
    const int num_secondary_students = total_students - num_elementary_students;
    return {static_cast<int>(std::ceil(num_elementary_students / params.ratios.max_students_per_elementary_school)),
            static_cast<int>(std::ceil(num_secondary_students / params.ratios.max_students_per_secondary_school)),
            num_elementary_students, num_secondary_students};
}

} // namespace

CityInfrastructure CityInfrastructure::calculate(int population, const RegionParameters& params)
{
    CityInfrastructure infra;

    infra.num_households_hh_size = household_sizes(population, params);

    const auto [num_workplaces, num_worker] = calc_num_workplaces_and_worker(population, params);
    infra.num_workplaces                    = std::max(1, num_workplaces);
    infra.num_worker                        = num_worker;

    const auto [num_elem_schools, num_sec_schools, num_elem_students, num_sec_students] =
        calc_num_elem_and_sec_schools(population, params);
    infra.num_elementary_schools         = std::max(1, num_elem_schools);
    infra.num_secondary_schools          = std::max(1, num_sec_schools);
    infra.num_persons_elementary_schools = num_elem_students;
    infra.num_persons_secondary_schools  = num_sec_students;

    infra.num_hospitals = 1;
    infra.num_icus      = 1;

    infra.num_stores = std::max(
        1, static_cast<int>(std::ceil(population * params.ratios.retail_stores_per_1000_people / 1000.0)));
    infra.num_events = std::max(1, static_cast<int>(std::ceil(population / params.ratios.people_per_event)));

    return infra;
}

} // namespace CityParameters
