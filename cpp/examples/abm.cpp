/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
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
#include "epidemiology/abm/abm.h"
/**
 * Determine the infection state of a person at the beginning of the simulation.
 * The infection states are chosen randomly. They are distributed according to the probabilites set in the example.
 * @return random infection state
 */
epi::InfectionState determine_infection_state(double exposed, double infected, double carrier, double recovered)
{
    double susceptible          = 1 - exposed - infected - carrier - recovered;
    std::vector<double> weights = {susceptible,  exposed,       carrier,      infected / 2,
                                   infected / 2, recovered / 2, recovered / 2};
    uint32_t state              = epi::DiscreteDistribution<size_t>::get_instance()(weights);
    return (epi::InfectionState)state;
}

/**
 * Add people to the world.
 * The age structur is equal to the age structure in Germany, 2019.
 * The distribution of household sizes approximates the distribution of household sizes in 2019.
 */
void create_people(epi::World& world, double num_total_people, double exposed_pct, double infected_pct,
                   double carrier_pct, double recovered_pct)
{

    // Create a household with arbitrary many people.
    // The vector age_groups contains the age groups of the people that are part of the household.
    auto create_household = [&](const std::vector<epi::AbmAgeGroup>& age_groups) {
        // Create home
        auto home = world.add_location(epi::LocationType::Home);
        // Genenerate people with defined age groups and assign them their home
        for (auto age_group : age_groups) {
            auto& p = world.add_person(
                home, determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct), age_group);
            p.set_assigned_location(home);
        }
    };

    // Create Children
    // Assumption: Children live in families with one or two children. 50% of the families have 1 child, 50% have two children.
    // 4% of the population is younger than 5
    int num_people_age0to4 = num_total_people * 0.04;
    // THe increment is 3 since in each loop 3 children are added.
    for (auto i = 0; i < num_people_age0to4; i = i + 3) {
        //Families with 2 children younger than 5
        create_household({epi::AbmAgeGroup::Age0to4, epi::AbmAgeGroup::Age0to4, epi::AbmAgeGroup::Age35to59,
                          epi::AbmAgeGroup::Age35to59});
        //Families with 1 child younger than 5
        create_household({epi::AbmAgeGroup::Age0to4, epi::AbmAgeGroup::Age15to34, epi::AbmAgeGroup::Age15to34});
    }

    // 9% of population are childreen between 5 and 14 years.
    int num_people_age5to14 = num_total_people * 0.09;
    // The ncrement is 3 since in each loop 3 children are added.
    for (auto i = 0; i < num_people_age5to14; i = i + 3) {
        //Families with 2 children older than 5
        create_household({epi::AbmAgeGroup::Age5to14, epi::AbmAgeGroup::Age5to14, epi::AbmAgeGroup::Age35to59,
                          epi::AbmAgeGroup::Age35to59});
        //Families with 1 child older than 5
        create_household({epi::AbmAgeGroup::Age5to14, epi::AbmAgeGroup::Age15to34, epi::AbmAgeGroup::Age35to59});
    }

    //Add adults
    // 23% of the population is between 15 and 34 years old.
    int num_people_age15to24 = num_total_people * 0.23;
    // Assumptions: Approx. 25% of people between 15 and 34 live with their children (they have been already added above)
    // 8% lives with their parents, 36% lives alone and 31 % with their partner
    for (auto i = 0; i < num_people_age15to24 * 0.08; ++i) {
        //Person lives with parents
        create_household({epi::AbmAgeGroup::Age15to34, epi::AbmAgeGroup::Age35to59, epi::AbmAgeGroup::Age35to59});
    }
    for (auto i = 0; i < num_people_age15to24 * 0.36; ++i) {
        //Person lives alone
        create_household({epi::AbmAgeGroup::Age15to34});
    }
    // Increment is 2 since in every loop two people between 15 and 34 years are added.
    for (auto i = 0; i < num_people_age15to24 * 0.36; i = i + 2) {
        //Person lives with partner
        create_household({epi::AbmAgeGroup::Age15to34, epi::AbmAgeGroup::Age15to34});
    }

    // 35% of the population is between 35 and 59 years old.
    int num_people_age35to59 = num_total_people * 0.35;
    //Assumption: Approx. 33% of people between 35 and 59 lives with their younger children and approx. 11% lives with their teenager (have been already added above)
    //39% lives with a partner and 17% lives alone
    for (auto i = 0; i < num_people_age35to59 * 0.17; ++i) {
        //Person lives alone
        create_household({epi::AbmAgeGroup::Age35to59});
    }
    // Increment increases by 2 since in each loop 2 people between 35 and 39 are added.
    for (auto i = 0; i < num_people_age35to59 * 0.39; i = i + 2) {
        //Person lives with partner
        create_household({epi::AbmAgeGroup::Age35to59, epi::AbmAgeGroup::Age35to59});
    }

    // Add senior citizen
    // 22% of the population is between 60 and 79 years
    int num_people_age60to79 = num_total_people * 0.22;
    // Assumptions: 66 procent of people between 60 and 79 live with a partner, 33% lives alone
    // The increment is 3 since in every loop 3 people between 60 and 79 are added.
    for (auto i = 0; i < num_people_age60to79; i = i + 3) {
        //Two persons
        create_household({epi::AbmAgeGroup::Age60to79, epi::AbmAgeGroup::Age60to79});
        //Single person
        create_household({epi::AbmAgeGroup::Age60to79});
    }

    // 7% of the population is older than 80 years
    int num_people_age80plus = num_total_people * 0.07;
    // Assumptions: 15% of people of age 80+ live in a home for senior citizens, 24% live with a partner, 51% live alone.
    // In average, in an home for senior citizens live around 80 people.
    // In the home for senior citizens, effective contacs are 10.
    int counter;
    int num_added_persons = 0;
    while (num_added_persons < num_people_age80plus * 0.15) {
        // create home for senior citizens
        counter   = 0;
        auto home = world.add_location(epi::LocationType::Home);
        world.get_individualized_location(home).get_infection_parameters().set<epi::EffectiveContacts>(10);
        // add 80 persons that live in the home
        while (counter < 80 && num_added_persons < num_people_age80plus * 0.15) {
            auto& adult =
                world.add_person(home, determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct),
                                 epi::AbmAgeGroup::Age80plus);
            adult.set_assigned_location(home);
            counter++;
            num_added_persons++;
        }
    }
    // Live alone
    for (auto i = 0; i < num_people_age80plus * 0.51; ++i) {
        create_household({epi::AbmAgeGroup::Age80plus});
    }
    // Live with a partner - per loop two people are added
    for (auto i = 0; i < num_people_age80plus * 0.24; i = i + 2) {
        create_household({epi::AbmAgeGroup::Age80plus, epi::AbmAgeGroup::Age80plus});
    }
}

/**
 * Add locations to the world and assign locations to the people.
 */
void create_assign_locations(epi::World& world)
{
    // Add one social event with 100 effective contacts.
    // Effective contacs limit the number of people that a person can infect while being at this location.
    // People have to get tested in the 2 days before the event
    auto event = world.add_location(epi::LocationType::SocialEvent);
    world.get_individualized_location(event).get_infection_parameters().set<epi::EffectiveContacts>(100);
    world.get_individualized_location(event).set_testing_scheme(epi::days(2), 1);

    // Add hospital and ICU with 5 effective contacs.
    auto hospital = world.add_location(epi::LocationType::Hospital);
    world.get_individualized_location(hospital).get_infection_parameters().set<epi::EffectiveContacts>(5);
    auto icu = world.add_location(epi::LocationType::ICU);
    world.get_individualized_location(icu).get_infection_parameters().set<epi::EffectiveContacts>(5);

    // Add schools, workplaces and shops.
    // At every school there are 600 students. The effective contacs are 40.
    // Students have to get tested once a week.
    // At every workplace work 100 people (needs to be varified), effective contacts are 40.
    // People can get tested at work (and do this with 0.5 probability).
    // Add one supermarked per 15.000 people, effective constacts are assumed to be 20.
    auto shop = world.add_location(epi::LocationType::BasicsShop);
    world.get_individualized_location(shop).get_infection_parameters().set<epi::EffectiveContacts>(20);

    auto school = world.add_location(epi::LocationType::School);
    world.get_individualized_location(school).get_infection_parameters().set<epi::EffectiveContacts>(40);
    world.get_individualized_location(school).set_testing_scheme(epi::days(7), 1);


    auto work = world.add_location(epi::LocationType::Work);
    world.get_individualized_location(work).get_infection_parameters().set<epi::EffectiveContacts>(40);
    world.get_individualized_location(work).set_testing_scheme(epi::days(7), 0.5);
    int counter_school = 0;
    int counter_work   = 0;
    int counter_shop   = 0;
    //Assign locations to the people
    auto persons = world.get_persons();
    for (auto& person : persons) {
        //assign shop and event
        person.set_assigned_location(event);
        person.set_assigned_location(shop);
        counter_shop++;
        //assign hospital and ICU
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
        //assign work/school to people depending on their age
        if (person.get_age() == epi::AbmAgeGroup::Age5to14) {
            person.set_assigned_location(school);
            counter_school++;
        }
        if (person.get_age() == epi::AbmAgeGroup::Age15to34 || person.get_age() == epi::AbmAgeGroup::Age35to59) {
            person.set_assigned_location(work);
            counter_work++;
        }
        //add new school/work/shop if needed
        if (counter_school == 600) {
            counter_school = 0;
            school         = world.add_location(epi::LocationType::School);
            world.get_individualized_location(school).get_infection_parameters().set<epi::EffectiveContacts>(40);
        }
        if (counter_work == 100) {
            counter_work = 0;
            work         = world.add_location(epi::LocationType::Work);
            world.get_individualized_location(work).get_infection_parameters().set<epi::EffectiveContacts>(40);
        }
        if (counter_shop == 15000) {
            counter_shop = 0;
            shop         = world.add_location(epi::LocationType::BasicsShop);
            world.get_individualized_location(shop).get_infection_parameters().set<epi::EffectiveContacts>(20);
        }
    }
}

int main()
{
    //epi::set_log_level(epi::LogLevel::warn);

    // Set seeds of previous run for debugging:
    // epi::thread_local_rng().seed({...});

    //Print seeds to be able to use them again for debugging:
    //printf("Seeds: ");
    //for (auto s : epi::thread_local_rng().get_seeds()) {
    //    printf("%u, ", s);
    //}
    //printf("\n");

    //Parameters
    //total number of people
    double num_total_people = 50000;

    //assumed percentage of infection state at the beginning of the simulation
    double exposed_pct = 0.01, infected_pct = 0.008, carrier_pct = 0.005, recovered_pct = 0.001;

    //Set global infection parameters (similar to infection parameters in SECIR model) and initialize the world
    epi::GlobalInfectionParameters abm_params;
    abm_params.set<epi::IncubationPeriod>({{epi::AbmAgeGroup::Count}, 4.});
    abm_params.set<epi::SusceptibleToExposedByCarrier>({{epi::AbmAgeGroup::Count}, 0.02});
    abm_params.set<epi::SusceptibleToExposedByInfected>({{epi::AbmAgeGroup::Count}, 0.02});
    abm_params.set<epi::CarrierToInfected>({{epi::AbmAgeGroup::Count}, 0.15});
    abm_params.set<epi::CarrierToRecovered>({{epi::AbmAgeGroup::Count}, 0.15});
    abm_params.set<epi::InfectedToRecovered>({{epi::AbmAgeGroup::Count}, 0.2});
    abm_params.set<epi::InfectedToSevere>({{epi::AbmAgeGroup::Count}, 0.03});
    abm_params.set<epi::SevereToRecovered>({{epi::AbmAgeGroup::Count}, 0.1});
    abm_params.set<epi::SevereToCritical>({{epi::AbmAgeGroup::Count}, 0.1});
    abm_params.set<epi::CriticalToRecovered>({{epi::AbmAgeGroup::Count}, 0.02});
    abm_params.set<epi::CriticalToDead>({{epi::AbmAgeGroup::Count}, 0.06});
    abm_params.set<epi::RecoveredToSusceptible>({{epi::AbmAgeGroup::Count}, 0.});

    auto world = epi::World(abm_params);

    //Add people to the world.
    create_people(world, num_total_people, exposed_pct, infected_pct, carrier_pct, recovered_pct);

    //Add locations and assign locations to the people.
    create_assign_locations(world);

    auto t0         = epi::TimePoint(0);
    auto t_lockdown = epi::TimePoint(0) + epi::days(20);
    auto tmax       = epi::TimePoint(0) + epi::days(60);

    // During the lockdown, 60% of people work from home and schools are closed for 90% of students.
    // Social events are very rare.
    epi::set_home_office(t_lockdown, 0.25, world.get_migration_parameters());
    epi::set_school_closure(t_lockdown, 0.9, world.get_migration_parameters());
    epi::close_social_events(t_lockdown, 0.9, world.get_migration_parameters());

    auto sim = epi::AbmSimulation(t0, std::move(world));

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
    auto f_abm = fopen("abm.txt", "w");
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
