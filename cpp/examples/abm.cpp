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
#include "abm/abm.h"
#include "abm/household.h"
#include "abm/world_builder.h"
#include <cstdio>
/**
 * Determine the infection state of a person at the beginning of the simulation.
 * The infection states are chosen randomly. They are distributed according to the probabilites set in the example.
 * @return random infection state
 */
mio::InfectionState determine_infection_state(double exposed, double infected, double carrier, double recovered)
{
    double susceptible          = 1 - exposed - infected - carrier - recovered;
    std::vector<double> weights = {susceptible,  exposed,       carrier,      infected / 2,
                                   infected / 2, recovered / 2, recovered / 2};
    uint32_t state              = mio::DiscreteDistribution<size_t>::get_instance()(weights);
    return (mio::InfectionState)state;
}

/**
 * Calculates a vector in which each entry describes the amount of people living in the corresponding household.
 * This is done with equal distribution and if the number of people is not divisible by number of households the last one gets the rest. E.g. number_of_people = 10, number_of_households = 3. Then the vector household_sizes = {3,3,4}.
 * @param number_of_people The total amount of people to be distributed.
 * @param number_of_households The total amount of households.
 * @return A vector with the size of each household.
 */
std::vector<int> last_household_gets_the_rest(int number_of_people, int number_of_households)
{
    std::vector<int> household_sizes(number_of_households, 0);
    int avarage_household_size_round_down = number_of_people / number_of_households; //int rounds down.
    int people_left                       = number_of_people -
                      avarage_household_size_round_down *
                          number_of_households; // People left if everyone got the same rounded down amount of people.
    for (auto i = 0; i < number_of_households - 1; i++) {
        household_sizes.at(i) = avarage_household_size_round_down;
    }
    household_sizes.at(number_of_households - 1) =
        avarage_household_size_round_down + people_left; // Last one gets the people which would've been left out.
    return household_sizes;
}

/**
 * Constructs a household group which has a single member to represent them all, e.g. all people have the same age distribution.
 * @param age_dist A vector with the amount of people in each age group
 * @param number_of_people The total amount of people living in this household group.
 * @param number_of_hh The number of households in this household group.
 * @return householdGroup A Class Household Group.
 */
mio::HouseholdGroup make_uniform_households(const mio::HouseholdMember& member, int number_of_people, int number_of_hh)
{

    // The size of each household is calculated in a vector household_size_list.
    auto households_size_list = last_household_gets_the_rest(number_of_people, number_of_hh);

    auto householdGroup = mio::HouseholdGroup();
    for (auto& household_size : households_size_list) {
        auto household = mio::Household();
        household.add_members(member, household_size); // Add members according to the amount of people in the list.
        householdGroup.add_households(household, 1); // Add the household to the household group.
    }
    return householdGroup;
}

/**
 * Constructs a household group with families.
 * @param child Child Household Member.
 * @param parent Parent Household Member.
 * @param random Random Household Member. This is for the rest Group where no exact age distribution can be found.
 * @param number_of_persons_in_household Amount of people in this household
 * @param number_of_full_familes Amount of full families, e.g. two parents and (number_of_persons_in_household - 2) children.
 * @param number_of_half_familes Amount of half families, e.g. one parent and (number_of_persons_in_household - 1) children.
 * @param number_of_other_familes number_of_persons_in_household random persons.
 * @return A Household group.
 */
mio::HouseholdGroup make_homes_with_families(const mio::HouseholdMember& child, const mio::HouseholdMember& parent,
                                             const mio::HouseholdMember& random, int number_of_persons_in_household,
                                             int number_of_full_familes, int number_of_half_familes,
                                             int number_of_other_familes)
{

    auto private_household_group = mio::HouseholdGroup();

    // Add full families.
    auto household_full = mio::Household();
    household_full.add_members(child, number_of_persons_in_household - 2);
    household_full.add_members(parent, 2);
    private_household_group.add_households(household_full, number_of_full_familes);

    // Add half families.
    auto household_half = mio::Household();
    household_half.add_members(child, number_of_persons_in_household - 1);
    household_half.add_members(parent, 1);
    private_household_group.add_households(household_half, number_of_half_familes);

    // Add other families.
    if (number_of_persons_in_household < 5) {
        auto household_others = mio::Household();
        household_others.add_members(random, number_of_persons_in_household);
        private_household_group.add_households(household_others, number_of_other_familes);
    }
    else if (number_of_persons_in_household == 5) {
        // For 5 and more people in one household we have to distribute the rest onto the left over households.
        int people_left_size5 = 545;

        auto households_size_list = last_household_gets_the_rest(people_left_size5, number_of_other_familes);

        auto household_rest = mio::HouseholdGroup();
        for (auto& household_size : households_size_list) {
            auto household = mio::Household();
            household.add_members(random, household_size); // Add members according to the amount of people in the list.
            household_rest.add_households(household, 1); // Add the household to the household group.
        }
    }
    return private_household_group;
}

void create_world_from_statistical_data(mio::World& world)
{

    /** The data is taken from
     * https://www-genesis.destatis.de/genesis/online?operation=statistic&levelindex=0&levelid=1627908577036&code=12211#abreadcrumb
     * All numbers are in 1000.
     * Destatis divides the Households into community households and private households.
     * Community Households are: Refugee, Disabled, Retirement and Others. We have an explicit age distribution, amount of households and amount of people for them but not the exact amount of people in each household.
     * The private Households are divided with respect to the amount of people living in each household. For a one person household we have the exact age distribution. For the rest we have data about which kind of family lives in them. The different kinds of families are: A family with two parents and the rest are children, a family with one parent and the rest are children and  "other" families with no exact data about their age.
    */

    // Refugee
    auto refugee = mio::HouseholdMember();
    refugee.set_age_weight(mio::AbmAgeGroup::Age0to4, 25);
    refugee.set_age_weight(mio::AbmAgeGroup::Age5to14, 12);
    refugee.set_age_weight(mio::AbmAgeGroup::Age15to34, 25);
    refugee.set_age_weight(mio::AbmAgeGroup::Age35to59, 9);
    refugee.set_age_weight(mio::AbmAgeGroup::Age60to79, 1);
    refugee.set_age_weight(mio::AbmAgeGroup::Age80plus, 1);
    int refugee_number_of_people     = 74;
    int refugee_number_of_households = 12;
    auto refugeeGroup = make_uniform_households(refugee, refugee_number_of_people, refugee_number_of_households);

    add_household_group_to_world(world, refugeeGroup);

    // Disabled
    auto disabled = mio::HouseholdMember();
    disabled.set_age_weight(mio::AbmAgeGroup::Age0to4, 2);
    disabled.set_age_weight(mio::AbmAgeGroup::Age5to14, 6);
    disabled.set_age_weight(mio::AbmAgeGroup::Age15to34, 13);
    disabled.set_age_weight(mio::AbmAgeGroup::Age35to59, 42);
    disabled.set_age_weight(mio::AbmAgeGroup::Age60to79, 97);
    disabled.set_age_weight(mio::AbmAgeGroup::Age80plus, 32);
    int disabled_number_of_people     = 194;
    int disabled_number_of_households = 8;

    auto disabledGroup = make_uniform_households(disabled, disabled_number_of_people, disabled_number_of_households);

    add_household_group_to_world(world, disabledGroup);

    // Retirement
    auto retired = mio::HouseholdMember();
    retired.set_age_weight(mio::AbmAgeGroup::Age15to34, 1);
    retired.set_age_weight(mio::AbmAgeGroup::Age35to59, 30);
    retired.set_age_weight(mio::AbmAgeGroup::Age60to79, 185);
    retired.set_age_weight(mio::AbmAgeGroup::Age80plus, 530);
    int retirement_number_of_people     = 744;
    int retirement_number_of_households = 16;

    auto retirementGroup =
        make_uniform_households(retired, retirement_number_of_people, retirement_number_of_households);

    add_household_group_to_world(world, retirementGroup);

    // Others
    auto other = mio::HouseholdMember();
    other.set_age_weight(mio::AbmAgeGroup::Age0to4, 30);
    other.set_age_weight(mio::AbmAgeGroup::Age5to14, 40);
    other.set_age_weight(mio::AbmAgeGroup::Age15to34, 72);
    other.set_age_weight(mio::AbmAgeGroup::Age35to59, 40);
    other.set_age_weight(mio::AbmAgeGroup::Age60to79, 30);
    other.set_age_weight(mio::AbmAgeGroup::Age80plus, 10);
    int others_number_of_people     = 222;
    int others_number_of_households = 20;

    auto otherGroup = make_uniform_households(other, others_number_of_people, others_number_of_households);

    add_household_group_to_world(world, otherGroup);

    // One Person Household (we have exact age data about this)
    auto one_person_household_member = mio::HouseholdMember();
    one_person_household_member.set_age_weight(mio::AbmAgeGroup::Age15to34, 4364);
    one_person_household_member.set_age_weight(mio::AbmAgeGroup::Age35to59, 7283);
    one_person_household_member.set_age_weight(mio::AbmAgeGroup::Age60to79, 4100);
    one_person_household_member.set_age_weight(mio::AbmAgeGroup::Age80plus, 1800);
    int one_person_number_of_people     = 15387;
    int one_person_number_of_households = 15387;

    auto onePersonGroup = make_uniform_households(one_person_household_member, one_person_number_of_people,
                                                  one_person_number_of_households);

    add_household_group_to_world(world, onePersonGroup);

    // For more than 1 family households we need families. These are parents and children and randoms (which are distributed like the data we have for these households).
    auto child = mio::HouseholdMember(); // A child is 50/50% 0-4 or 5-14.
    child.set_age_weight(mio::AbmAgeGroup::Age0to4, 1);
    child.set_age_weight(mio::AbmAgeGroup::Age5to14, 1);

    auto parent = mio::HouseholdMember(); // A child is 40/40/20% 15-34, 35-59 or 60-79.
    parent.set_age_weight(mio::AbmAgeGroup::Age15to34, 2);
    parent.set_age_weight(mio::AbmAgeGroup::Age35to59, 2);
    parent.set_age_weight(mio::AbmAgeGroup::Age60to79, 1);

    auto random = mio::HouseholdMember(); // Randoms are distributed according to the left over persons.
    random.set_age_weight(mio::AbmAgeGroup::Age0to4, 5000);
    random.set_age_weight(mio::AbmAgeGroup::Age5to14, 6000);
    random.set_age_weight(mio::AbmAgeGroup::Age15to34, 14943);
    random.set_age_weight(mio::AbmAgeGroup::Age35to59, 22259);
    random.set_age_weight(mio::AbmAgeGroup::Age60to79, 11998);
    random.set_age_weight(mio::AbmAgeGroup::Age80plus, 5038);

    // Two person households
    int two_person_full_families  = 11850;
    int two_person_half_families  = 1765;
    int two_person_other_families = 166;
    auto twoPersonHouseholds      = make_homes_with_families(child, parent, random, 2, two_person_full_families,
                                                        two_person_half_families, two_person_other_families);
    add_household_group_to_world(world, twoPersonHouseholds);

    // Three person households
    int three_person_full_families  = 4155;
    int three_person_half_families  = 662;
    int three_person_other_families = 175;
    auto threePersonHouseholds      = make_homes_with_families(child, parent, random, 3, three_person_full_families,
                                                          three_person_half_families, three_person_other_families);
    add_household_group_to_world(world, threePersonHouseholds);

    // Four person households
    int four_person_full_families  = 3551;
    int four_person_half_families  = 110;
    int four_person_other_families = 122;
    auto fourPersonHouseholds      = make_homes_with_families(child, parent, random, 4, four_person_full_families,
                                                         four_person_half_families, four_person_other_families);
    add_household_group_to_world(world, fourPersonHouseholds);

    // Five plus person households
    int fiveplus_person_full_families  = 1245;
    int fiveplus_person_half_families  = 80;
    int fiveplus_person_other_families = 82;
    auto fivePlusPersonHouseholds =
        make_homes_with_families(child, parent, random, 5, fiveplus_person_full_families, fiveplus_person_half_families,
                                 fiveplus_person_other_families);
    add_household_group_to_world(world, fivePlusPersonHouseholds);
}

/**
 * Add locations to the world and assign locations to the people.
 */
void create_assign_locations(mio::World& world)
{
    // Add one social event with 100 effective contacts.
    // Effective contacs limit the number of people that a person can infect while being at this location.
    // People have to get tested in the 2 days before the event
    auto event = world.add_location(mio::LocationType::SocialEvent);
    world.get_individualized_location(event).get_infection_parameters().set<mio::EffectiveContacts>(100);
    world.get_individualized_location(event).set_testing_scheme(mio::days(2), 1);

    // Add hospital and ICU with 5 effective contacs.
    auto hospital = world.add_location(mio::LocationType::Hospital);
    world.get_individualized_location(hospital).get_infection_parameters().set<mio::EffectiveContacts>(5);
    auto icu = world.add_location(mio::LocationType::ICU);
    world.get_individualized_location(icu).get_infection_parameters().set<mio::EffectiveContacts>(5);

    // Add schools, workplaces and shops.
    // At every school there are 600 students. The effective contacs are 40.
    // Students have to get tested once a week.
    // At every workplace work 100 people (needs to be varified), effective contacts are 40.
    // People can get tested at work (and do this with 0.5 probability).
    // Add one supermarked per 15.000 people, effective constacts are assumed to be 20.
    auto shop = world.add_location(mio::LocationType::BasicsShop);
    world.get_individualized_location(shop).get_infection_parameters().set<mio::EffectiveContacts>(20);

    auto school = world.add_location(mio::LocationType::School);
    world.get_individualized_location(school).get_infection_parameters().set<mio::EffectiveContacts>(40);
    world.get_individualized_location(school).set_testing_scheme(mio::days(7), 1);

    auto work = world.add_location(mio::LocationType::Work);
    world.get_individualized_location(work).get_infection_parameters().set<mio::EffectiveContacts>(40);
    world.get_individualized_location(work).set_testing_scheme(mio::days(7), 0.5);
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
        if (person.get_age() == mio::AbmAgeGroup::Age5to14) {
            person.set_assigned_location(school);
            counter_school++;
        }
        if (person.get_age() == mio::AbmAgeGroup::Age15to34 || person.get_age() == mio::AbmAgeGroup::Age35to59) {
            person.set_assigned_location(work);
            counter_work++;
        }
        //add new school/work/shop if needed
        if (counter_school == 600) {
            counter_school = 0;
            school         = world.add_location(mio::LocationType::School);
            world.get_individualized_location(school).get_infection_parameters().set<mio::EffectiveContacts>(40);
        }
        if (counter_work == 100) {
            counter_work = 0;
            work         = world.add_location(mio::LocationType::Work);
            world.get_individualized_location(work).get_infection_parameters().set<mio::EffectiveContacts>(40);
        }
        if (counter_shop == 15000) {
            counter_shop = 0;
            shop         = world.add_location(mio::LocationType::BasicsShop);
            world.get_individualized_location(shop).get_infection_parameters().set<mio::EffectiveContacts>(20);
        }
    }
}

/**
 * Assign an infection state to each person.
 */

void assign_infection_state(mio::World& world, double exposed_pct, double infected_pct, double carrier_pct,
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
    /*
    Eigen::VectorXd people(6);
    people << 10, 18, 30, 18, 18, 6;
    Eigen::MatrixXd contact_matrix(6, 6);
    contact_matrix << 2.8, 4.3, 1.7, 4.4, 4.3, 1.5, 2.3888888, 3.777777, 6.222222, 3.8333333, 2.277777, 0.5, 0.5666666,
        3.733333, 10.4, 3.23333, 1.03333, 0.03333, 2.44444, 3.833333, 5.388888, 4, 2.72222, 0.611111, 2.388888,
        2.277777, 1.72222, 2.722222, 6.66666, 3.22222, 2.5, 1.5, 0.166666, 1.833333, 9.66666, 3.33333;
    int num_locs              = 5;
    Eigen::VectorXd size_locs = Eigen::VectorXd::Constant(5, 20);
    */

    /*
    Eigen::VectorXd people(2);
    people << 40, 24;
    Eigen::MatrixXd contact_matrix (2,2);
    contact_matrix << 24, 7, 11.6666666, 19.3333333;
    int num_locs = 2;
    Eigen::VectorXd size_locs = Eigen::VectorXd::Constant(2, 32);
     */

    Eigen::VectorXd people(6);
    people << 4091, 6983, 24611, 35921, 10238, 8319;

    Eigen::MatrixXd contact_matrix(6, 6);
    contact_matrix << 0.5170, 0.3997, 0.7957, 0.9958, 0.3239, 0.0428, 0.0632, 0.9121, 0.3254, 0.4731, 0.2355, 0.0148,
        0.0336, 0.1604, 1.7529, 0.8622, 0.1440, 0.0077, 0.0204, 0.1444, 0.5738, 1.2127, 0.3433, 0.0178, 0.0371, 0.0393,
        0.4171, 0.9666, 0.7495, 0.0257, 0.0791, 0.0800, 0.3480, 0.5588, 0.2769, 0.0180;
    int num_locs              = 10;
    Eigen::VectorXd size_locs = Eigen::VectorXd::Constant(num_locs, people.sum() / num_locs);

    Eigen::VectorXd x_sol = mio::find_optimal_locations(people, num_locs, contact_matrix, size_locs);
    //std::cout << x_sol << std::endl;
    Eigen::MatrixXd average_contacs(6, 6);
    mio::compute_average_contact_matrix(average_contacs, x_sol, num_locs);
    std::cout << " Average Matrix: " << average_contacs << std::endl;
    std::cout << "Differenz zur gegebenen Matrix: " << (average_contacs - contact_matrix).norm() << std::endl;

    /*
    //mio::set_log_level(mio::LogLevel::warn);

    // Set seeds of previous run for debugging:
    // mio::thread_local_rng().seed({...});

    //Print seeds to be able to use them again for debugging:
    //printf("Seeds: ");
    //for (auto s : mio::thread_local_rng().get_seeds()) {
    //    printf("%u, ", s);
    //}
    //printf("\n");

    // Assumed percentage of infection state at the beginning of the simulation.
    double exposed_pct = 0.01, infected_pct = 0.08, carrier_pct = 0.05, recovered_pct = 0.01;

    //Set global infection parameters (similar to infection parameters in SECIR model) and initialize the world
    mio::GlobalInfectionParameters abm_params;

    abm_params.set<mio::IncubationPeriod>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 4.});
    abm_params.set<mio::SusceptibleToExposedByCarrier>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.02});
    abm_params.set<mio::SusceptibleToExposedByInfected>(
        {{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.02});
    abm_params.set<mio::CarrierToInfected>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.15});
    abm_params.set<mio::CarrierToRecovered>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.15});
    abm_params.set<mio::InfectedToRecovered>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.2});
    abm_params.set<mio::InfectedToSevere>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.03});
    abm_params.set<mio::SevereToRecovered>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.1});
    abm_params.set<mio::SevereToCritical>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.1});
    abm_params.set<mio::CriticalToRecovered>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.02});
    abm_params.set<mio::CriticalToDead>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.06});
    abm_params.set<mio::RecoveredToSusceptible>({{mio::AbmAgeGroup::Count, mio::VaccinationState::Count}, 0.1});

    // Set each parameter for vaccinated people
    for (auto age = mio::Index<mio::AbmAgeGroup>(0); age < mio::AbmAgeGroup::Count; ++age) {
        abm_params.get<mio::IncubationPeriod>()[{age, mio::VaccinationState::Vaccinated}]               = 4.;
        abm_params.get<mio::SusceptibleToExposedByCarrier>()[{age, mio::VaccinationState::Vaccinated}]  = 0.02;
        abm_params.get<mio::SusceptibleToExposedByInfected>()[{age, mio::VaccinationState::Vaccinated}] = 0.02;
        abm_params.get<mio::CarrierToRecovered>()[{age, mio::VaccinationState::Vaccinated}]             = 0.15;
        abm_params.get<mio::InfectedToRecovered>()[{age, mio::VaccinationState::Vaccinated}]            = 0.15;
        abm_params.get<mio::InfectedToSevere>()[{age, mio::VaccinationState::Vaccinated}]               = 0.05;
        abm_params.get<mio::SevereToRecovered>()[{age, mio::VaccinationState::Vaccinated}]              = 0.05;
        abm_params.get<mio::SevereToCritical>()[{age, mio::VaccinationState::Vaccinated}]               = 0.005;
        abm_params.get<mio::CriticalToRecovered>()[{age, mio::VaccinationState::Vaccinated}]            = 0.05;
        abm_params.get<mio::CriticalToDead>()[{age, mio::VaccinationState::Vaccinated}]                 = 0.005;
        abm_params.get<mio::RecoveredToSusceptible>()[{age, mio::VaccinationState::Vaccinated}]         = 0.05;
    }

    auto world = mio::World(abm_params);

    // Create the world object from statistical data.
    create_world_from_statistical_data(world);

    // Assign an infection state to each person.
    assign_infection_state(world, exposed_pct, infected_pct, carrier_pct, recovered_pct);

    // Add locations and assign locations to the people.
    create_assign_locations(world);

    auto t0         = mio::TimePoint(0);
    auto t_lockdown = mio::TimePoint(0) + mio::days(20);
    auto tmax       = mio::TimePoint(0) + mio::days(60);

    // During the lockdown, 60% of people work from home and schools are closed for 90% of students.
    // Social events are very rare.
    mio::set_home_office(t_lockdown, 0.25, world.get_migration_parameters());
    mio::set_school_closure(t_lockdown, 0.9, world.get_migration_parameters());
    mio::close_social_events(t_lockdown, 0.9, world.get_migration_parameters());

    auto sim = mio::AbmSimulation(t0, std::move(world));

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
     */
}
