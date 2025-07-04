#include "../include/city_builder.h"
#include "../include/constants.h"
#include "../include/parameter_setter.h"
#include "../include/city_parameters.h"
#include <random>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>
#include <iostream>
#include <iomanip>

mio::IOResult<mio::abm::World> CityBuilder::build_world(const CityConfig& config)
{
    auto world = mio::abm::World(num_age_groups);
    set_parameters(world.parameters);

    // Get infrastructure configuration based on German demographic data
    auto infra = config.infrastructure();

    // Healthcare facilities
    std::vector<mio::abm::LocationId> hospitals;
    std::vector<mio::abm::LocationId> icus;

    for (int i = 0; i < infra.num_hospitals; ++i) {
        hospitals.push_back(world.add_location(mio::abm::LocationType::Hospital));
    }
    for (int i = 0; i < infra.num_icus; ++i) {
        icus.push_back(world.add_location(mio::abm::LocationType::ICU));
    }

    // Create all location types based on German infrastructure ratios
    auto households = create_households(world, infra.num_households);
    auto workplaces = create_workplaces(world, infra.num_workplaces);
    auto schools    = create_schools(world, infra.num_elementary_schools + infra.num_secondary_schools);
    auto shops      = create_shops(world, infra.num_grocery_stores, infra.num_pharmacies, infra.num_general_stores);
    auto events =
        create_events(world, infra.num_large_events, infra.num_small_events + infra.num_restaurants + infra.num_bars);

    // Assign people to locations using German demographic distribution
    BOOST_OUTCOME_TRY(assign_people_to_locations(world, households, workplaces, schools, shops, events, hospitals, icus,
                                                 config.total_population));

    set_local_parameters(world);

    return mio::success(std::move(world));
}

std::vector<mio::abm::LocationId> CityBuilder::create_households(mio::abm::World& world, int num_households)
{
    std::vector<mio::abm::LocationId> household_ids;
    household_ids.reserve(num_households);

    for (int i = 0; i < num_households; ++i) {
        auto household_id = world.add_location(mio::abm::LocationType::Home);
        household_ids.push_back(household_id);
    }

    return household_ids;
}

std::vector<mio::abm::LocationId> CityBuilder::create_workplaces(mio::abm::World& world, int num_workplaces)
{
    std::vector<mio::abm::LocationId> workplace_ids;
    workplace_ids.reserve(num_workplaces);

    for (int i = 0; i < num_workplaces; ++i) {
        auto workplace_id = world.add_location(mio::abm::LocationType::Work);
        workplace_ids.push_back(workplace_id);
    }

    return workplace_ids;
}

std::vector<mio::abm::LocationId> CityBuilder::create_schools(mio::abm::World& world, int num_schools)
{
    std::vector<mio::abm::LocationId> school_ids;
    school_ids.reserve(num_schools);

    for (int i = 0; i < num_schools; ++i) {
        auto school_id = world.add_location(mio::abm::LocationType::School);
        school_ids.push_back(school_id);
    }

    return school_ids;
}

std::vector<mio::abm::LocationId> CityBuilder::create_events(mio::abm::World& world, int large, int small)
{
    std::vector<mio::abm::LocationId> event_ids;
    event_ids.reserve(large + small);

    // Create large events (concert halls, stadiums)
    for (int i = 0; i < large; ++i) {
        auto event_id = world.add_location(mio::abm::LocationType::SocialEvent);
        event_ids.push_back(event_id);
    }

    // Create small events (restaurants, bars, community centers)
    for (int i = 0; i < small; ++i) {
        auto event_id = world.add_location(mio::abm::LocationType::SocialEvent);
        event_ids.push_back(event_id);
    }

    return event_ids;
}

std::vector<mio::abm::LocationId> CityBuilder::create_shops(mio::abm::World& world, int grocery_stores, int pharmacies,
                                                            int general_stores)
{
    std::vector<mio::abm::LocationId> shop_ids;
    shop_ids.reserve(grocery_stores + pharmacies + general_stores);

    // Create grocery stores
    for (int i = 0; i < grocery_stores; ++i) {
        auto shop_id = world.add_location(mio::abm::LocationType::BasicsShop);
        shop_ids.push_back(shop_id);
    }

    // Create pharmacies
    for (int i = 0; i < pharmacies; ++i) {
        auto shop_id = world.add_location(mio::abm::LocationType::BasicsShop);
        shop_ids.push_back(shop_id);
    }

    // Create general stores
    for (int i = 0; i < general_stores; ++i) {
        auto shop_id = world.add_location(mio::abm::LocationType::BasicsShop);
        shop_ids.push_back(shop_id);
    }

    return shop_ids;
}

mio::IOResult<void> CityBuilder::assign_people_to_locations(
    mio::abm::World& world, const std::vector<mio::abm::LocationId>& households,
    const std::vector<mio::abm::LocationId>& workplaces, const std::vector<mio::abm::LocationId>& schools,
    const std::vector<mio::abm::LocationId>& shops, const std::vector<mio::abm::LocationId>& events,
    const std::vector<mio::abm::LocationId>& hospitals, const std::vector<mio::abm::LocationId>& icus,
    int total_population)
{
    auto world_rng = world.get_rng();

    // assign everyone to a hospital and ICU randomly
    std::uniform_int_distribution<> hospital_dist(0, hospitals.size() - 1);
    std::uniform_int_distribution<> icu_dist(0, icus.size() - 1);

    for (int i = 0; i < total_population; ++i) {
        auto age_group = assign_age_group_from_demographics(world_rng);

        // Assign to household using German household distribution
        auto household_id = assign_household(households, world_rng);
        auto& person      = world.add_person(household_id, age_group);
        person.set_assigned_location(household_id);

        // Assign healthcare locations
        person.set_assigned_location(hospitals[hospital_dist(world_rng)]);
        person.set_assigned_location(icus[icu_dist(world_rng)]);

        // Assign to school based on German school attendance rates
        if (should_attend_school(age_group, world_rng)) {
            auto school_id = schools[i % schools.size()];
            person.set_assigned_location(school_id);
        }

        // Assign to workplace based on German employment rates
        if (should_be_employed(age_group, world_rng) && !workplaces.empty()) {
            auto workplace_id = workplaces[i % workplaces.size()];
            person.set_assigned_location(workplace_id);
        }

        // Assign to shops (everyone needs basic supplies)
        if (!shops.empty()) {
            auto shop_id = shops[i % shops.size()];
            person.set_assigned_location(shop_id);
        }

        // Assign to events (social activities)
        if (!events.empty()) {
            auto event_id = events[i % events.size()];
            person.set_assigned_location(event_id);
        }
    }

    return mio::success();
}

mio::AgeGroup CityBuilder::assign_age_group_from_demographics(mio::RandomNumberGenerator& gen)
{
    // Use German age distribution from 2023 census data
    std::discrete_distribution<> age_dist(CityParameters::GERMAN_AGE_DISTRIBUTION.begin(),
                                          CityParameters::GERMAN_AGE_DISTRIBUTION.end());

    int age_group_index = age_dist(gen);
    return mio::AgeGroup(age_group_index);
}

mio::abm::LocationId CityBuilder::assign_household(const std::vector<mio::abm::LocationId>& households,
                                                   mio::RandomNumberGenerator& gen)
{
    // Simple assignment for now - could be enhanced with household size distribution
    std::uniform_int_distribution<> household_dist(0, households.size() - 1);
    return households[household_dist(gen)];
}

bool CityBuilder::should_attend_school(const mio::AgeGroup& age_group, mio::RandomNumberGenerator& gen)
{
    auto it = CityParameters::SCHOOL_ATTENDANCE_RATES.find(age_group.get());
    if (it == CityParameters::SCHOOL_ATTENDANCE_RATES.end()) {
        return false;
    }

    std::uniform_real_distribution<> rate_dist(0.0, 1.0);
    return rate_dist(gen) < it->second;
}

bool CityBuilder::should_be_employed(const mio::AgeGroup& age_group, mio::RandomNumberGenerator& gen)
{
    auto it = CityParameters::EMPLOYMENT_RATES.find(age_group.get());
    if (it == CityParameters::EMPLOYMENT_RATES.end()) {
        return false;
    }

    std::uniform_real_distribution<> rate_dist(0.0, 1.0);
    return rate_dist(gen) < it->second;
}

void CityBuilder::print_city_summary(const CityConfig& config)
{
    auto infra = config.infrastructure();

    std::cout << "\n" << std::string(50, '=') << "\n";
    std::cout << "     GERMAN CITY SIMULATION SUMMARY\n";
    std::cout << std::string(50, '=') << "\n\n";

    // Population Overview
    std::cout << "📊 POPULATION OVERVIEW\n";
    std::cout << "Total Population: " << config.total_population << "\n";
    std::cout << "Expected Age Distribution (based on German 2023 census):\n";
    std::cout << "  • 0-4 years:   " << std::fixed << std::setprecision(1)
              << (CityParameters::GERMAN_AGE_DISTRIBUTION[0] * 100) << "% ("
              << static_cast<int>(config.total_population * CityParameters::GERMAN_AGE_DISTRIBUTION[0]) << " people)\n";
    std::cout << "  • 5-14 years:  " << (CityParameters::GERMAN_AGE_DISTRIBUTION[1] * 100) << "% ("
              << static_cast<int>(config.total_population * CityParameters::GERMAN_AGE_DISTRIBUTION[1]) << " people)\n";
    std::cout << "  • 15-34 years: " << (CityParameters::GERMAN_AGE_DISTRIBUTION[2] * 100) << "% ("
              << static_cast<int>(config.total_population * CityParameters::GERMAN_AGE_DISTRIBUTION[2]) << " people)\n";
    std::cout << "  • 35-59 years: " << (CityParameters::GERMAN_AGE_DISTRIBUTION[3] * 100) << "% ("
              << static_cast<int>(config.total_population * CityParameters::GERMAN_AGE_DISTRIBUTION[3]) << " people)\n";
    std::cout << "  • 60-79 years: " << (CityParameters::GERMAN_AGE_DISTRIBUTION[4] * 100) << "% ("
              << static_cast<int>(config.total_population * CityParameters::GERMAN_AGE_DISTRIBUTION[4]) << " people)\n";
    std::cout << "  • 80+ years:   " << (CityParameters::GERMAN_AGE_DISTRIBUTION[5] * 100) << "% ("
              << static_cast<int>(config.total_population * CityParameters::GERMAN_AGE_DISTRIBUTION[5])
              << " people)\n\n";

    // Housing
    std::cout << "🏠 HOUSING\n";
    std::cout << "Households: " << infra.num_households << "\n";
    std::cout << "Average household size: " << std::setprecision(2)
              << static_cast<double>(config.total_population) / infra.num_households << " people\n";
    std::cout << "German average: " << CityParameters::AVERAGE_HOUSEHOLD_SIZE << " people\n\n";

    // Employment
    std::cout << "💼 EMPLOYMENT\n";
    std::cout << "Workplaces: " << infra.num_workplaces << "\n";
    std::cout << "Average employees per workplace: " << std::setprecision(1)
              << static_cast<double>(config.total_population) * CityParameters::InfrastructureRatios::EMPLOYMENT_RATE /
                     infra.num_workplaces
              << "\n";
    std::cout << "Expected employment rate: " << std::setprecision(1)
              << (CityParameters::InfrastructureRatios::EMPLOYMENT_RATE * 100) << "%\n\n";

    // Education
    std::cout << "🎓 EDUCATION\n";
    std::cout << "Elementary schools: " << infra.num_elementary_schools << "\n";
    std::cout << "Secondary schools: " << infra.num_secondary_schools << "\n";
    std::cout << "Total schools: " << (infra.num_elementary_schools + infra.num_secondary_schools) << "\n";
    std::cout << "Students per elementary school: ~" << std::setprecision(0)
              << CityParameters::InfrastructureRatios::STUDENTS_PER_ELEMENTARY_SCHOOL << "\n";
    std::cout << "Students per secondary school: ~" << std::setprecision(0)
              << CityParameters::InfrastructureRatios::STUDENTS_PER_SECONDARY_SCHOOL << "\n\n";

    // Healthcare
    std::cout << "🏥 HEALTHCARE\n";
    std::cout << "Hospitals: " << infra.num_hospitals << "\n";
    std::cout << "ICU units: " << infra.num_icus << "\n";
    std::cout << "People per hospital: " << std::setprecision(0)
              << static_cast<double>(config.total_population) / infra.num_hospitals << "\n";
    std::cout << "People per ICU: " << std::setprecision(0)
              << static_cast<double>(config.total_population) / infra.num_icus << "\n\n";

    // Retail & Services
    std::cout << "🛒 RETAIL & SERVICES\n";
    std::cout << "Grocery stores: " << infra.num_grocery_stores << "\n";
    std::cout << "Pharmacies: " << infra.num_pharmacies << "\n";
    std::cout << "General stores: " << infra.num_general_stores << "\n";
    std::cout << "Total retail locations: "
              << (infra.num_grocery_stores + infra.num_pharmacies + infra.num_general_stores) << "\n";
    std::cout << "People per grocery store: " << std::setprecision(0)
              << static_cast<double>(config.total_population) / infra.num_grocery_stores << "\n\n";

    // Social & Entertainment
    std::cout << "🎉 SOCIAL & ENTERTAINMENT\n";
    std::cout << "Restaurants: " << infra.num_restaurants << "\n";
    std::cout << "Bars/Pubs: " << infra.num_bars << "\n";
    std::cout << "Large event venues: " << infra.num_large_events << " (stadiums, concert halls)\n";
    std::cout << "Small event venues: " << infra.num_small_events << " (community centers, clubs)\n";
    std::cout << "People per restaurant: " << std::setprecision(0)
              << static_cast<double>(config.total_population) / infra.num_restaurants << "\n\n";

    // Infrastructure Summary
    std::cout << "📈 INFRASTRUCTURE DENSITY\n";
    std::cout << "Total locations: "
              << (infra.num_households + infra.num_workplaces + infra.num_elementary_schools +
                  infra.num_secondary_schools + infra.num_hospitals + infra.num_icus + infra.num_grocery_stores +
                  infra.num_pharmacies + infra.num_general_stores + infra.num_restaurants + infra.num_bars +
                  infra.num_large_events + infra.num_small_events)
              << "\n";
    std::cout << "Locations per 1000 people: " << std::setprecision(1)
              << static_cast<double>(infra.num_households + infra.num_workplaces + infra.num_elementary_schools +
                                     infra.num_secondary_schools + infra.num_hospitals + infra.num_icus +
                                     infra.num_grocery_stores + infra.num_pharmacies + infra.num_general_stores +
                                     infra.num_restaurants + infra.num_bars + infra.num_large_events +
                                     infra.num_small_events) *
                     1000.0 / config.total_population
              << "\n\n";

    std::cout << "📍 DATA SOURCES\n";
    std::cout << "• German Federal Statistical Office (Destatis) 2023\n";
    std::cout << "• German Trade Association (HDE) 2023\n";
    std::cout << "• German Hotel and Restaurant Association (DEHOGA) 2023\n";

    std::cout << "\n" << std::string(50, '=') << "\n\n";
}