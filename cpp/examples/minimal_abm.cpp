/*
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Khoa Nguyen
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
#include "abm/abm.h"
#include "abm/household.h"
#include <cstdio>

/**
 * Determine the infection state of a person at the beginning of the simulation.
 * The infection states are chosen randomly. They are distributed according to the probabilites set in the example.
 * @return random infection state
 */
mio::abm::InfectionState determine_infection_state(double exposed, double infected, double carrier, double recovered)
{
    double susceptible          = 1 - exposed - infected - carrier - recovered;
    std::vector<double> weights = {susceptible,  exposed,       carrier,      infected / 2,
                                   infected / 2, recovered / 2, recovered / 2};
    uint32_t state              = (uint32_t)mio::DiscreteDistribution<size_t>::get_instance()(weights);
    return (mio::abm::InfectionState)state;
}

void create_world(mio::abm::World& world)
{

    // There are 3 households for each hosehold group.
    int n_households = 3;

    // One-person household group with one from each age group.
    auto onePersonHousehold_group = mio::abm::HouseholdGroup();
    auto onePersonHousehold_member = mio::abm::HouseholdMember();
    onePersonHousehold_member.set_age_weight(mio::abm::AgeGroup::Age15to34, 1);
    onePersonHousehold_member.set_age_weight(mio::abm::AgeGroup::Age35to59, 1);
    onePersonHousehold_member.set_age_weight(mio::abm::AgeGroup::Age60to79, 1);
    for (int counter = 0; counter < n_households; counter++) {
        auto onePersonHousehold = mio::abm::Household();
        onePersonHousehold.add_members(onePersonHousehold_member, 1);
        onePersonHousehold_group.add_households(onePersonHousehold, 1);
    }

    // For more than 1 family households we need families. These are parents and children and randoms (which are distributed like the data we have for these households).
    auto child = mio::abm::HouseholdMember(); // A child is 50/50% 0-4 or 5-14.
    child.set_age_weight(mio::abm::AgeGroup::Age0to4, 1);
    child.set_age_weight(mio::abm::AgeGroup::Age5to14, 1);

    auto parent = mio::abm::HouseholdMember(); // A parent is 50/50% 15-34, 35-59 or 60-79.
    parent.set_age_weight(mio::abm::AgeGroup::Age15to34, 2);
    parent.set_age_weight(mio::abm::AgeGroup::Age35to59, 2);
    
    
    // Two-person household with one parent and one child. 
    auto twoPersonHousehold_group = mio::abm::HouseholdGroup();
    for (int counter = 0; counter < n_households; counter++) {
        auto twoPersonHousehold_full = mio::abm::Household();
        twoPersonHousehold_full.add_members(child, 1);
        twoPersonHousehold_full.add_members(parent, 1);
        twoPersonHousehold_group.add_households(twoPersonHousehold_full, 2);  
    }
    add_household_group_to_world(world, twoPersonHousehold_group);

    // Three-person household with two parent and one child. 
    auto threePersonHousehold_group = mio::abm::HouseholdGroup();
    for (int counter = 0; counter < n_households; counter++) {
        auto threePersonHousehold_full = mio::abm::Household();
        threePersonHousehold_full.add_members(child, 1);
        threePersonHousehold_full.add_members(parent, 2);
        threePersonHousehold_group.add_households(threePersonHousehold_full, 2);
    }
   add_household_group_to_world(world, threePersonHousehold_group);
}

/**
 * Add locations to the world and assign locations to the people.
 */
void create_assign_locations(mio::abm::World& world)
{
    // Add one social event with 5 maximum contacts.
    // Maximum contacs limit the number of people that a person can infect while being at this location.
    // People have to get tested in the 2 days before the event
    auto event = world.add_location(mio::abm::LocationType::SocialEvent);
    world.get_individualized_location(event).get_infection_parameters().set<mio::abm::MaximumContacts>(5);

    std::vector<mio::abm::LocationType> test_at_social_event = {mio::abm::LocationType::SocialEvent};
    auto testing_criteria =
        std::vector<mio::abm::TestingCriteria>{mio::abm::TestingCriteria({}, test_at_social_event, {})};
    auto testing_min_time = mio::abm::days(1);
    auto start_date       = mio::abm::TimePoint(0);
    auto end_date         = mio::abm::TimePoint(0) + mio::abm::days(20);
    auto probability      = 1.0;
    auto test_type        = mio::abm::AntigenTest();

    auto testing_scheme =
        mio::abm::TestingScheme(testing_criteria, testing_min_time, start_date, end_date, test_type, probability);

    world.get_testing_strategy().add_testing_scheme(testing_scheme);

    // Add hospital and ICU with 5 maximum contacs.
    auto hospital = world.add_location(mio::abm::LocationType::Hospital);
    world.get_individualized_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto icu = world.add_location(mio::abm::LocationType::ICU);
    world.get_individualized_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);

    // Add schools, workplaces and shops.
    // At every school, the maximum contacts are 20.
    // At every workplace, maximum contacts are 10.
    // People can get tested at work (and do this with 0.5 probability).
    // Add one supermarked, maximum constacts are assumed to be 20.
    auto shop = world.add_location(mio::abm::LocationType::BasicsShop);
    world.get_individualized_location(shop).get_infection_parameters().set<mio::abm::MaximumContacts>(20);

    auto school = world.add_location(mio::abm::LocationType::School);
    world.get_individualized_location(school).get_infection_parameters().set<mio::abm::MaximumContacts>(20);

    auto work = world.add_location(mio::abm::LocationType::Work);
    world.get_individualized_location(work).get_infection_parameters().set<mio::abm::MaximumContacts>(10);

    // Assign locations to the people
    auto persons = world.get_persons();
    for (auto& person : persons) {
        //assign shop and event
        person.set_assigned_location(event);
        person.set_assigned_location(shop);
        //assign hospital and ICU
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
        //assign work/school to people depending on their age
        if (person.get_age() == mio::abm::AgeGroup::Age5to14) {
            person.set_assigned_location(school);
        }
        if (person.get_age() == mio::abm::AgeGroup::Age15to34 || person.get_age() == mio::abm::AgeGroup::Age35to59) {
            person.set_assigned_location(work);
        }
    }

    // add the testing schemes for work
    auto test_at_work = std::vector<mio::abm::LocationType>{mio::abm::LocationType::Work};
    auto testing_criteria_work =
        std::vector<mio::abm::TestingCriteria>{mio::abm::TestingCriteria({}, test_at_work, {})};
    testing_min_time = mio::abm::days(1);
    probability      = 0.5;
    auto testing_scheme_work =
        mio::abm::TestingScheme(testing_criteria_work, testing_min_time, start_date, end_date, test_type, probability);
    world.get_testing_strategy().add_testing_scheme(testing_scheme_work);
}

/**
 * Assign an infection state to each person.
 */

void assign_infection_state(mio::abm::World& world, double exposed_pct, double infected_pct, double carrier_pct,
                            double recovered_pct)
{
    auto persons = world.get_persons();
    for (auto& person : persons) {
        world.set_infection_state(person,
                                  determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct));
    }
}

int main()
{
    // Assumed percentage of infection state at the beginning of the simulation.
    double exposed_pct = 0.01, infected_pct = 0.08, carrier_pct = 0.05, recovered_pct = 0.01;

    // Set global infection parameters (similar to infection parameters in SECIR model) and initialize the world
    mio::abm::GlobalInfectionParameters infection_params;

    // Set same parameter for all age groups
    infection_params.get<mio::abm::IncubationPeriod>()               = 4.;
    infection_params.get<mio::abm::SusceptibleToExposedByCarrier>()  = 0.02;
    infection_params.get<mio::abm::SusceptibleToExposedByInfected>() = 0.02;
    infection_params.get<mio::abm::CarrierToInfected>()              = 0.15;
    infection_params.get<mio::abm::CarrierToRecovered>()             = 0.15;
    infection_params.get<mio::abm::InfectedToRecovered>()            = 0.2;
    infection_params.get<mio::abm::InfectedToSevere>()               = 0.03;
    infection_params.get<mio::abm::SevereToRecovered>()              = 0.1;
    infection_params.get<mio::abm::SevereToCritical>()               = 0.1;
    infection_params.get<mio::abm::CriticalToRecovered>()            = 0.02;
    infection_params.get<mio::abm::CriticalToDead>()                 = 0.06;
    infection_params.get<mio::abm::RecoveredToSusceptible>()         = 0.1;

    // Set parameters for vaccinated people of all age groups
    infection_params.get<mio::abm::IncubationPeriod>().slice(mio::abm::VaccinationState::Vaccinated) = 4.;
    infection_params.get<mio::abm::SusceptibleToExposedByCarrier>().slice(mio::abm::VaccinationState::Vaccinated) =
        0.02;
    infection_params.get<mio::abm::SusceptibleToExposedByInfected>().slice(mio::abm::VaccinationState::Vaccinated) =
        0.02;
    infection_params.get<mio::abm::CarrierToRecovered>().slice(mio::abm::VaccinationState::Vaccinated)     = 0.15;
    infection_params.get<mio::abm::InfectedToRecovered>().slice(mio::abm::VaccinationState::Vaccinated)    = 0.15;
    infection_params.get<mio::abm::InfectedToSevere>().slice(mio::abm::VaccinationState::Vaccinated)       = 0.05;
    infection_params.get<mio::abm::SevereToRecovered>().slice(mio::abm::VaccinationState::Vaccinated)      = 0.05;
    infection_params.get<mio::abm::SevereToCritical>().slice(mio::abm::VaccinationState::Vaccinated)       = 0.005;
    infection_params.get<mio::abm::CriticalToRecovered>().slice(mio::abm::VaccinationState::Vaccinated)    = 0.5;
    infection_params.get<mio::abm::CriticalToDead>().slice(mio::abm::VaccinationState::Vaccinated)         = 0.005;
    infection_params.get<mio::abm::RecoveredToSusceptible>().slice(mio::abm::VaccinationState::Vaccinated) = 0.05;

    auto world = mio::abm::World(infection_params);

    // Create the world object.
    create_world(world);

    // Assign an infection state to each person.
    assign_infection_state(world, exposed_pct, infected_pct, carrier_pct, recovered_pct);

    // Add locations and assign locations to the people.
    create_assign_locations(world);

    auto t0         = mio::abm::TimePoint(0);
    auto t_lockdown = mio::abm::TimePoint(0) + mio::abm::days(10);
    auto tmax       = mio::abm::TimePoint(0) + mio::abm::days(30);

    // During the lockdown, 60% of people work from home and schools are closed for 90% of students.
    // Social events are very rare.
    mio::abm::set_home_office(t_lockdown, 0.25, world.get_migration_parameters());
    mio::abm::set_school_closure(t_lockdown, 0.9, world.get_migration_parameters());
    mio::abm::close_social_events(t_lockdown, 0.9, world.get_migration_parameters());

    auto sim = mio::abm::Simulation(t0, std::move(world));

    sim.advance(tmax);

    // The results are saved in a table with 9 rows.
    // The first row is t = time, the others correspond to the number of people with a certain infection state at this time:
    // S = Susceptible, E = Exposed, C= Carrier, I= Infected, I_s = Infected_Severe,
    // I_c = Infected_Critical, R_C = Recovered_Carrier, R_I = Recovered_Infected, D = Dead
    // E.g. the following gnuplot skrips plots detected infections and deaths.
    // plot "abm.txt" using 1:5 with lines title "infected (detected)", "abm.txt" using 1:10 with lines title "dead"
    // set xlabel "days"
    // set ylabel "number of people"
    // set title "ABM Example"
    // set output "abm.png"
    // set terminal png
    // replot
    auto f_abm = fopen("minimal_abm.txt", "w");
    fprintf(f_abm, "# t S E C I I_s I_c R_C R_I D\n");
    for (auto i = 0; i < sim.get_result().get_num_time_points(); ++i) {
        fprintf(f_abm, "%f ", sim.get_result().get_time(i));
        auto v = sim.get_result().get_value(i);
        for (auto j = 0; j < v.size(); ++j) {
            fprintf(f_abm, "%f", v[j]);
            if (j < v.size() - 1) {
                fprintf(f_abm, " ");
            }
        }
        if (i < sim.get_result().get_num_time_points() - 1) {
            fprintf(f_abm, "\n");
        }
    }
    fclose(f_abm);
}