#include "../include/city_builder.h"
#include "../include/constants.h"
#include "../include/parameter_setter.h"
#include <random>
#include <algorithm>

mio::IOResult<mio::abm::World> CityBuilder::build_city(const CityConfig& config)
{
    auto world = mio::abm::World(num_age_groups);
    set_parameters(world.parameters);

    // Create all location types
    auto households = create_households(world, config.num_households);
    auto workplaces = create_workplaces(world, config.num_workplaces);
    auto schools    = create_schools(world, config.num_schools);
    auto events     = create_events(world, config.num_large_events, config.num_small_events);

    // Assign people to locations
    BOOST_OUTCOME_TRY(assign_people_to_locations(world, households, workplaces, schools, config.total_population));

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

    // Create large events
    for (int i = 0; i < large; ++i) {
        auto event_id = world.add_location(mio::abm::LocationType::SocialEvent);
        event_ids.push_back(event_id);
    }

    // Create small events
    for (int i = 0; i < small; ++i) {
        auto event_id = world.add_location(mio::abm::LocationType::BasicsShop);
        event_ids.push_back(event_id);
    }

    return event_ids;
}

mio::IOResult<void> CityBuilder::assign_people_to_locations(mio::abm::World& world,
                                                            const std::vector<mio::abm::LocationId>& households,
                                                            const std::vector<mio::abm::LocationId>& workplaces,
                                                            const std::vector<mio::abm::LocationId>& schools,
                                                            int total_population)
{
    std::random_device rd;
    std::mt19937 gen(rd());

    for (int i = 0; i < total_population; ++i) {
        auto person_id = mio::abm::PersonId(i);
        auto age_group = assign_age_group_from_demographics(i);

        // Assign to household
        auto household_id = households[i % households.size()];
        auto person =
            mio::abm::Person(world.get_individualized_location(mio::abm::LocationType::Home, person_id), age_group);

        // TODO: Assign work/school based on age group and census data
        // For now, simple assignment
        if (age_group == age_group_5_to_14 && !schools.empty()) {
            // School age - assign to school
            auto school_id = schools[i % schools.size()];
            person.set_assigned_location(mio::abm::LocationType::School);
        }
        else if ((age_group == age_group_15_to_34 || age_group == age_group_35_to_59) && !workplaces.empty()) {
            // Working age - assign to workplace
            auto workplace_id = workplaces[i % workplaces.size()];
            person.set_assigned_location(mio::abm::LocationType::Work);
        }

        world.add_person(person_id, std::move(person));
    }

    return mio::success();
}

mio::AgeGroup CityBuilder::assign_age_group_from_demographics(int person_index)
{
    // TODO: Use real census data to assign age groups
    // For now, simple distribution
    std::vector<double> age_distribution = {0.05, 0.15, 0.25, 0.35, 0.15, 0.05}; // 0-4, 5-14, 15-34, 35-59, 60-79, 80+

    int age_group_index = person_index % 6;
    return mio::AgeGroup(age_group_index);
}