//
//  world_test.cpp
//  
//
//  Created by Sascha Korf on 04.08.21.
//

#include "epidemiology/abm/abm.h"

struct infectionStateParameters{
    double exposed, infected, carrier, recovered;
};

struct poolToChooseFrom{
    std::vector<int> pool_ages;
    int total_pool_size;
};

struct privateHHParameters{
    poolToChooseFrom children_pool, parents_pool, rest_pool;
    int index_for_extras_in_5plus_households;
};

struct communityHhParameters{
    poolToChooseFrom age_pool;
    int number_of_people, number_of_households;
    int mean_rounded_down, deviation_from_mean;
    int number_of_extra_persons_left, index_extra_persons;
};

/*
 * Struct for defining the demografic parameters important for the building of world with these parameters.
 */
struct demograficParameters
{
    int total_population;
    
    // people in different hh.
    int people_in_priv_hh, people_in_comm_hh;
    int people_in_retirement_hh, people_in_disabled_hh, people_in_refugee_hh, people_in_other_hh;
    std::vector<int> people_in_priv_hh_each_size; // First entry is the number of people in 1 person HH, second for 2 person HH, etc.
    
    // Age distributions.
    std::vector<int> age_distr_priv_hh, age_distr_comm_hh; // General distribution into the age distribution
    std::vector<int> age_distr_priv_hh_size1; // For Household with 1 person the age distribution is available
    std::vector<int> age_distr_comm_hh_retirement, age_distr_comm_hh_disabled, age_distr_comm_hh_refugee, age_distr_comm_hh_other; // For Community Households the age distribution is available.
    
    // Number of Household's of each Type.
    int number_priv_hh, number_comm_hh;
    std::vector<int> number_priv_hh_each_size; // First entry is the number of Households of 1 person HH, second of 2 person HH, etc.
    int number_retirement_hh, number_disabled_hh, number_refugee_hh, number_other_hh; // These are preferably dividable by 2.
    
    // Number of families in each Household size. Full means two parents plus kids, half means one parent and rest are kids. Rest means all other combinations which are not full or half.
    std::vector<int> full_families_private_hh_each_size, half_families_private_hh_each_size, rest_families_private_hh_each_size;
};

/**
 * Determine the infection state of a person at the beginning of the simulation.
 * The infection states are chosen randomly. They are distributed according to the probabilites set in the example.
 * @return random infection state
 */
epi::InfectionState determine_infection_state(infectionStateParameters *inf_state_par){
    double susceptible = 1 - inf_state_par->exposed - inf_state_par->infected - inf_state_par->carrier - inf_state_par->recovered;
    std::vector<double> weights = {susceptible, inf_state_par->exposed, inf_state_par->carrier, (inf_state_par->infected)/2, (inf_state_par->infected)/2, (inf_state_par->recovered)/2, (inf_state_par->recovered)/2};
    uint32_t state = epi::DiscreteDistribution<size_t>::get_instance()(weights);
    return (epi::InfectionState) state;
}
/*
 * Function that takes a vector of integers and returns a vector which has the same length. This vector has the percentage of the
 * whole vector sum in each row.
 * @param absoluteValueVec Vector with the absolute Values
 * @return  Vector with percent
 */
std::vector<double> get_percentage_of_sum(std::vector<int> absoluteValueVec){
    int sumOfValues;
    std::vector<double> getPercentageOfSum;
    // Sum the values
    sumOfValues = std::accumulate(absoluteValueVec.begin(), absoluteValueVec.end(), 0);
    // Get the percentage of each entry of the total sum.
    getPercentageOfSum.reserve(absoluteValueVec.size());
    for (int i = 0; i<absoluteValueVec.size(); i++) {
        getPercentageOfSum.push_back(((double)absoluteValueVec.at(i))/((double)sumOfValues));
    }
    // Return
    return getPercentageOfSum;
}

/*
 * Function that sets up the demografic parameters
 */
void setup_parameters(demograficParameters *demo_par){
    // Set general parameters
    // Set total population
    demo_par->total_population = 81848;
    
    // Set number of people in private / community HH
    demo_par->people_in_priv_hh = 80615;
    demo_par->people_in_comm_hh = 1233;
    
    // Set number of People in Community HH
    demo_par->people_in_retirement_hh = 744;
    demo_par->people_in_disabled_hh = 194;
    demo_par->people_in_refugee_hh = 74;
    demo_par->people_in_other_hh = 222;
    
    // Set number of People in private Homes
    demo_par->people_in_priv_hh_each_size.resize(5);
    demo_par->people_in_priv_hh_each_size.at(0) = 15387;
    demo_par->people_in_priv_hh_each_size.at(1) = 27562;
    demo_par->people_in_priv_hh_each_size.at(2) = 14856;
    demo_par->people_in_priv_hh_each_size.at(3) = 15132;
    demo_par->people_in_priv_hh_each_size.at(4) = 7170;
    
    // Set age distribution where it is possible
    // General Private HH
    demo_par->age_distr_priv_hh.resize(6);
    demo_par->age_distr_priv_hh.at(0) = 5000;
    demo_par->age_distr_priv_hh.at(1) = 6000;
    demo_par->age_distr_priv_hh.at(2) = 19307;
    demo_par->age_distr_priv_hh.at(3) = 29542;
    demo_par->age_distr_priv_hh.at(4) = 16098;
    demo_par->age_distr_priv_hh.at(5) = 6838;
    
    // general Community HH
    demo_par->age_distr_comm_hh.resize(6);
    demo_par->age_distr_comm_hh.at(0) = 30;
    demo_par->age_distr_comm_hh.at(1) = 50;
    demo_par->age_distr_comm_hh.at(2) = 70;
    demo_par->age_distr_comm_hh.at(3) = 100;
    demo_par->age_distr_comm_hh.at(4) = 450;
    demo_par->age_distr_comm_hh.at(5) = 533;
    
    // Private Households
    // 1-member HH
    demo_par->age_distr_priv_hh_size1.resize(6);
    demo_par->age_distr_priv_hh_size1.at(0) = 0;
    demo_par->age_distr_priv_hh_size1.at(1) = 0;
    demo_par->age_distr_priv_hh_size1.at(2) = 4364;
    demo_par->age_distr_priv_hh_size1.at(3) = 7283;
    demo_par->age_distr_priv_hh_size1.at(4) = 4100;
    demo_par->age_distr_priv_hh_size1.at(5) = 1800;
    
    // Retirement Homes
    demo_par->age_distr_comm_hh_retirement.resize(6);
    demo_par->age_distr_comm_hh_retirement.at(0) = 0;
    demo_par->age_distr_comm_hh_retirement.at(1) = 0;
    demo_par->age_distr_comm_hh_retirement.at(2) = 1;
    demo_par->age_distr_comm_hh_retirement.at(3) = 30;
    demo_par->age_distr_comm_hh_retirement.at(4) = 185;
    demo_par->age_distr_comm_hh_retirement.at(5) = 530;
    
    // Disabled Homes
    demo_par->age_distr_comm_hh_disabled.resize(6);
    demo_par->age_distr_comm_hh_disabled.at(0) = 2;
    demo_par->age_distr_comm_hh_disabled.at(1) = 6;
    demo_par->age_distr_comm_hh_disabled.at(2) = 13;
    demo_par->age_distr_comm_hh_disabled.at(3) = 42;
    demo_par->age_distr_comm_hh_disabled.at(4) = 97;
    demo_par->age_distr_comm_hh_disabled.at(5) = 32;
    
    // Refugee Homes
    demo_par->age_distr_comm_hh_refugee.resize(6);
    demo_par->age_distr_comm_hh_refugee.at(0) = 25;
    demo_par->age_distr_comm_hh_refugee.at(1) = 12;
    demo_par->age_distr_comm_hh_refugee.at(2) = 25;
    demo_par->age_distr_comm_hh_refugee.at(3) = 9;
    demo_par->age_distr_comm_hh_refugee.at(4) = 1;
    demo_par->age_distr_comm_hh_refugee.at(5) = 1;
    
    // Other Homes
    demo_par->age_distr_comm_hh_other.resize(6);
    demo_par->age_distr_comm_hh_other.at(0) = 30;
    demo_par->age_distr_comm_hh_other.at(1) = 40;
    demo_par->age_distr_comm_hh_other.at(2) = 72;
    demo_par->age_distr_comm_hh_other.at(3) = 40;
    demo_par->age_distr_comm_hh_other.at(4) = 30;
    demo_par->age_distr_comm_hh_other.at(5) = 10;
    
    // Number of HH
    demo_par->number_priv_hh_each_size.resize(5);
    demo_par->number_priv_hh_each_size.at(0) = 15387;
    demo_par->number_priv_hh_each_size.at(1) = 13781;
    demo_par->number_priv_hh_each_size.at(2) = 4952;
    demo_par->number_priv_hh_each_size.at(3) = 3783;
    demo_par->number_priv_hh_each_size.at(4) = 1434;
    
    // Number of each community HH, preferably dividable by 2.
    demo_par->number_retirement_hh = 16;
    demo_par->number_disabled_hh = 8;
    demo_par->number_refugee_hh  = 6;
    demo_par->number_other_hh = 20;
    
    demo_par->number_priv_hh = std::accumulate(demo_par->number_priv_hh_each_size.begin(), demo_par->number_priv_hh_each_size.end(), 0);
    demo_par->number_comm_hh = demo_par->number_retirement_hh + demo_par->number_disabled_hh + demo_par->number_refugee_hh+demo_par->number_other_hh;
    
    //Number of full/half/rest families
    demo_par->full_families_private_hh_each_size.resize(5);
    demo_par->full_families_private_hh_each_size.at(0) = 0;
    demo_par->full_families_private_hh_each_size.at(1) = 11850;
    demo_par->full_families_private_hh_each_size.at(2) = 4155;
    demo_par->full_families_private_hh_each_size.at(3) = 3551;
    demo_par->full_families_private_hh_each_size.at(4) = 1245;
    
    demo_par->half_families_private_hh_each_size.resize(5);
    demo_par->half_families_private_hh_each_size.at(0) = 0;
    demo_par->half_families_private_hh_each_size.at(1) = 1765;
    demo_par->half_families_private_hh_each_size.at(2) = 662;
    demo_par->half_families_private_hh_each_size.at(3) = 110;
    demo_par->half_families_private_hh_each_size.at(4) = 80;
    
    demo_par->rest_families_private_hh_each_size.resize(5);
    demo_par->rest_families_private_hh_each_size.at(0) = 0;
    demo_par->rest_families_private_hh_each_size.at(1) = 166;
    demo_par->rest_families_private_hh_each_size.at(2) = 175;
    demo_par->rest_families_private_hh_each_size.at(3) = 122;
    demo_par->rest_families_private_hh_each_size.at(4) = 109;
}

void setup_infection_state_parameters(infectionStateParameters *infStatePar){
    // This function sets up the infection state parameters
    infStatePar->carrier = 0.005;
    infStatePar->exposed = 0.01;
    infStatePar->infected = 0.008;
    infStatePar->recovered = 0.001;
}

void substract_from_pool(poolToChooseFrom *pool, int count, int age_group){
    pool->pool_ages.at(age_group) -= count;
    pool->total_pool_size -= count;
}

// This function picks a age group out of a vector of age groups. It is weighted based on the number of people in each age group. E.g. the age group with the most people gets picked most likely. (Discretely distributed).
int pick_age_group_from_age_group_vector(std::vector<int> age_groups){
    std::vector<double> weights = get_percentage_of_sum(age_groups);
    uint32_t age_group = epi::DiscreteDistribution<size_t>::get_instance()(weights);
    return age_group;
}

// This function picks a person from a pool of people.
epi::AbmAgeGroup pick_age_from_pool(poolToChooseFrom *pool){
    int age_group = pick_age_group_from_age_group_vector(pool->pool_ages);
    substract_from_pool(pool,1,age_group);
    return (epi::AbmAgeGroup) age_group;
}

// Private households functions.
// This function picks the next available children or parent out of the pools and if they aren't available it chooses the next person out of a suitable age group.
epi::AbmAgeGroup pick_next_available_ph(privateHHParameters *ph_par, bool parents, bool children){
    
    int index_sufficient_agegroup=0;
    int age_group;
    
    if(children){
        for (auto i = 0; i<ph_par->children_pool.pool_ages.size(); i++) {
            if( ph_par->parents_pool.pool_ages.at(i) !=0 ){
                age_group = i;
                substract_from_pool(&ph_par->parents_pool,1,age_group);
                return (epi::AbmAgeGroup) age_group;
            }
            else if (ph_par->rest_pool.pool_ages.at(i) !=0 ) {
                age_group = i;
                substract_from_pool(&ph_par->rest_pool,1,age_group);
                return (epi::AbmAgeGroup) age_group;
            }
        }
        printf(" Too few persons in pool, please put more in it or code a function that generates some more.\n");
    }
    
    if(parents){
        for (auto i = 2; i<ph_par->parents_pool.pool_ages.size(); i++) {
            if( ph_par->children_pool.pool_ages.at(i) !=0 ){
                age_group = i;
                substract_from_pool(&ph_par->children_pool,1,age_group);
                return (epi::AbmAgeGroup) age_group;
            }
            else if (ph_par->rest_pool.pool_ages.at(i) !=0 ) {
                age_group = i;
                substract_from_pool(&ph_par->rest_pool,1,age_group);
                return (epi::AbmAgeGroup) age_group;
            }
        }
        printf(" Too few persons in pool, please put more in it or code a function that generates some more.\n");
    }
}

// This function picks a child / parent / or rest for private Households.
epi::AbmAgeGroup pick_age_group_for_ph(privateHHParameters *ph_par, bool child, bool parent){
    int age_group;
    
    if(child){
        if(ph_par->children_pool.total_pool_size > 0){
            return (epi::AbmAgeGroup) pick_age_from_pool(&ph_par->children_pool);
        }else{
            return pick_next_available_ph(ph_par, parent, child);
        }
    }else if (parent){
        if(ph_par->parents_pool.total_pool_size > 0){
            return (epi::AbmAgeGroup) pick_age_from_pool(&ph_par->parents_pool);
        }else{
            return pick_next_available_ph(ph_par, parent, child);
        }
    }else{
        if(ph_par->rest_pool.total_pool_size > 0){
            return (epi::AbmAgeGroup)  pick_age_from_pool(&ph_par->rest_pool);
        }else{
            printf("This shouldn't happen. Can't pick person.");
        }
    }
}

// This function puts all of the children / parents into rest for the "Rest" distribution after families have been distributed in each household.
void update_pools_after_family_distribution_ph(privateHHParameters *ph_par){
    for (auto i = 0; i < ph_par->rest_pool.pool_ages.size(); i++) {
        if(ph_par->children_pool.pool_ages.at(i) > 0){
            ph_par->rest_pool.pool_ages.at(i) += ph_par->children_pool.pool_ages.at(i);
            ph_par->children_pool.total_pool_size -=  ph_par->children_pool.pool_ages.at(i);
            ph_par->children_pool.pool_ages.at(i) = 0;
        }
        if(ph_par->parents_pool.pool_ages.at(i) > 0){
            ph_par->rest_pool.pool_ages.at(i) += ph_par->parents_pool.pool_ages.at(i);
            ph_par->parents_pool.total_pool_size -=  ph_par->parents_pool.pool_ages.at(i);
            ph_par->parents_pool.pool_ages.at(i) = 0;
        }
    }
}

void create_one_person_household_ph(epi::World& world, demograficParameters *demo_par, privateHHParameters *ph_par, infectionStateParameters *inf_state_par){
    // Create all the 1-person Household because the data is clear on this.
    int num_people_oneP_HH = demo_par->people_in_priv_hh_each_size.at(0);
    int number_people_this_agegroup = 0;
    int age_groups = 6;
    
    epi::AbmAgeGroup age_group;
    
    for(auto i = 0; i < age_groups; i++){
        number_people_this_agegroup = demo_par->age_distr_priv_hh_size1.at(i);
        if(number_people_this_agegroup > 0){
            age_group = (epi::AbmAgeGroup) i;
            for (auto j = 0; j < number_people_this_agegroup; j++) {
                auto home = world.add_location(epi::LocationType::Home);
                auto& p = world.add_person(home, determine_infection_state(inf_state_par), age_group);
                substract_from_pool(&ph_par->rest_pool, 1, i);
            }
        }
    }
}

// This function calcs how many extra people are added to rest of 5+ person HH.
int how_many_extras_ph(privateHHParameters *ph_par, demograficParameters *demo_par){
    int extras;
    int number_of_5plus_hh_to_compute = demo_par->rest_families_private_hh_each_size.at(4) - ph_par->index_for_extras_in_5plus_households - 1
    ;
    ++ph_par->index_for_extras_in_5plus_households;
    int remaining_extra_people = ph_par->rest_pool.total_pool_size - number_of_5plus_hh_to_compute * 5;

    if(remaining_extra_people>0){
        extras = remaining_extra_people/number_of_5plus_hh_to_compute;
        if(extras + 1 > remaining_extra_people){
            return remaining_extra_people;
        }else{
            return extras + 1;
        }
    }else{
        extras = 0;
        return extras;
    }
}

void create_full_families_of_ph(epi::World& world, demograficParameters *demo_par, privateHHParameters *ph_par, infectionStateParameters *inf_state_par, int household_size){
    // Full families mean two parents and children, e.g. household_size 4 equals 2 parents and 2 children,  household_size 5 means 2 parents and 3 children.
    if(household_size <2 || household_size > 5 )
        printf("This is not an acceptable housheold size");
    // First we create a household and two parents.
    // Household
    auto home = world.add_location(epi::LocationType::Home);
    // Two parents
    world.add_person(home, determine_infection_state(inf_state_par), pick_age_group_for_ph(ph_par, false, true));
    world.add_person(home, determine_infection_state(inf_state_par), pick_age_group_for_ph(ph_par, false, true));
    // The rest are children.
    int number_of_children = household_size-2;
    for (auto i = 0 ; i < number_of_children; i++) {
        auto child = pick_age_group_for_ph(ph_par, true, false);
        world.add_person(home, determine_infection_state(inf_state_par), child);
    }
}

void create_half_families_of_ph(epi::World& world, demograficParameters *demo_par, privateHHParameters *ph_par, infectionStateParameters *inf_state_par, int household_size){
    // Half families mean one parent and children, e.g. household_size 4 equals 1 parents and 3 children,  household_size 5 means 1 parents and 4 children
    if(household_size <2 || household_size > 5 )
        printf("This is not an acceptable househeold size");
    // First we create a household and two parents.
    // Household
    auto home = world.add_location(epi::LocationType::Home);
    // One parent
    auto parent1 = pick_age_group_for_ph(ph_par, true, false);
    auto& p = world.add_person(home, determine_infection_state(inf_state_par), parent1);
    // The rest are children.
    int number_of_children = household_size-1;
    for (auto i = 0 ; i < number_of_children; i++) {
        auto child = pick_age_group_for_ph(ph_par, false, true);
        world.add_person(home, determine_infection_state(inf_state_par), child);
    }
}

void create_rest_of_ph(epi::World& world, demograficParameters *demo_par, privateHHParameters *ph_par, infectionStateParameters *inf_state_par, int household_size){
    if(household_size <2 || household_size > 5 )
        printf("This is not an acceptable househeold size");
    auto home = world.add_location(epi::LocationType::Home);
    for (auto i = 0 ; i < household_size; i++) {
        auto person = pick_age_group_for_ph(ph_par, false, false);
        world.add_person(home, determine_infection_state(inf_state_par), person);
    }
    // If the household size is 5 there are extra persons in the "Rest" group so we distribute them into this group.
    if (household_size  == 5) {
        int extras = how_many_extras_ph(ph_par, demo_par);
        for (auto i = 0 ; i < extras; i++) {
            auto person = pick_age_group_for_ph(ph_par, false, false);
            world.add_person(home, determine_infection_state(inf_state_par), person);
        }
    }
}

void create_more_person_households_ph(epi::World& world, demograficParameters *demo_par, privateHHParameters *ph_par, infectionStateParameters *inf_state_par){
    /* Here we create the more than one person households.
     * We go through this like this: First we create "families", e.g. half and full families in each of the households of x size.
     * Then we create the rest with the rest of the pool.
     */
    
    // Fill the full families.
    for (auto household_size = 2; household_size < 6; household_size++){
        for (int i = 0; i<demo_par->full_families_private_hh_each_size.at(household_size-1); i++) {
            create_full_families_of_ph(world, demo_par,ph_par, inf_state_par, household_size);
        }
    }
    
    // Fill the half families.
    for (auto household_size = 2; household_size < 6; household_size++){
        for (int i = 0; i<demo_par->half_families_private_hh_each_size.at(household_size-1); i++) {
            create_half_families_of_ph(world, demo_par, ph_par, inf_state_par, household_size);
        }
    }
    
    update_pools_after_family_distribution_ph(ph_par);
    
    // Fill the rest of the households
    for (auto household_size = 2; household_size < 6; household_size++){
        for (int i = 0; i < demo_par->rest_families_private_hh_each_size.at(household_size-1); i++) {
            create_rest_of_ph(world, demo_par, ph_par, inf_state_par, household_size);
        }
    }
}

// This function sets up the private household pools (parent, child and rest) from the statistical data.
void setup_ph(demograficParameters *demo_par, privateHHParameters *ph_par){
    // This function takes into account how many are needed for the one person household whose age distribution is available. After this consideration Age group 15 to 34 is distributed 60/40 into parents/children. This is a guessed distribution.
    int rest_age_group_15to34 = demo_par->age_distr_priv_hh_size1.at(2);
    int rest_age_group_35to59 = demo_par->age_distr_priv_hh_size1.at(3);
    
    // Children Pool
    ph_par->children_pool.pool_ages.resize(6,0);
    ph_par->children_pool.pool_ages.at(0)=demo_par->age_distr_priv_hh.at(0);
    ph_par->children_pool.pool_ages.at(1)=demo_par->age_distr_priv_hh.at(1);
    ph_par->children_pool.pool_ages.at(2)=0.4*(demo_par->age_distr_priv_hh.at(2)-rest_age_group_15to34);
    ph_par->children_pool.total_pool_size = std::accumulate( ph_par->children_pool.pool_ages.begin(),  ph_par->children_pool.pool_ages.end(), 0);
    
    // Parents Pool
    ph_par->parents_pool.pool_ages.resize(6,0);
    ph_par->parents_pool.pool_ages.at(2)=0.6*(demo_par->age_distr_priv_hh.at(2)-rest_age_group_15to34);
    ph_par->parents_pool.pool_ages.at(3)=demo_par->age_distr_priv_hh.at(3)-rest_age_group_35to59;
    ph_par->parents_pool.total_pool_size = std::accumulate( ph_par->parents_pool.pool_ages.begin(),  ph_par->parents_pool.pool_ages.end(), 0);
    
    // Rest Pool
    ph_par->rest_pool.pool_ages.resize(6,0);
    for (auto i = 0; i <  ph_par->rest_pool.pool_ages.size(); i++) {
        ph_par->rest_pool.pool_ages.at(i)=demo_par->age_distr_priv_hh.at(i) - ph_par->children_pool.pool_ages.at(i)- ph_par->parents_pool.pool_ages.at(i);
    }
    ph_par->rest_pool.total_pool_size = std::accumulate( ph_par->rest_pool.pool_ages.begin(),  ph_par->rest_pool.pool_ages.end(), 0);
    
    //Index for extras
    ph_par->index_for_extras_in_5plus_households = 0;
}

void create_private_households(epi::World& world, demograficParameters *demo_par, infectionStateParameters *inf_state_par){
    privateHHParameters ph_par;
    setup_ph(demo_par, &ph_par);
    create_one_person_household_ph(world, demo_par, &ph_par, inf_state_par);
    create_more_person_households_ph(world, demo_par, &ph_par, inf_state_par);
}

// Community households functions.
int how_many_extras_ch(communityHhParameters *ch_par){
    int extras;
    int number_of_communites_to_compute = (ch_par->number_of_households/2)-ch_par->index_extra_persons; //Divided by two because we always add a pair and extras only are in one.
    
    ++ch_par->index_extra_persons;
    if(ch_par->number_of_extra_persons_left > 0){
        extras = ch_par->number_of_extra_persons_left/number_of_communites_to_compute;
        if(extras + 1 > ch_par->number_of_extra_persons_left){
            extras = ch_par->number_of_extra_persons_left;
            ch_par->number_of_extra_persons_left = 0;
            return extras;
        }else{
            extras = extras + 1;
            ch_par->number_of_extra_persons_left -= extras;
            return extras;
        }
    }else{
        extras = 0;
        return extras;
    }
}

void add_pair_ch(epi::World& world, communityHhParameters *ch_par, infectionStateParameters *inf_state_par){
    // This functions adds the number of households to the worlds. In the following way:
    // We need to distribute x people to y households. The avarage rounded down is in ch_par->mean_rounded_down.
    // We add 2 Households at the same time. First we calculate a random number z which is distributed with a standart distribution deviation of ch_par->deviation_from_mean and add mean_rounded_down-z and mean_rounded_down+z+extra.
    // Extra are the remaining people which need to be added after the mean is rounded down and is calculated in the function how many extras.
    std::normal_distribution<double> distribution(ch_par->mean_rounded_down,ch_par->deviation_from_mean);
    std::default_random_engine generator;
    int normal_distr_random_int = (int) distribution(generator);
    if(normal_distr_random_int > ch_par->mean_rounded_down){
        normal_distr_random_int = ch_par->mean_rounded_down - 1;
    }
    
    int hh_size1 = normal_distr_random_int;
    int hh_size2 =  ch_par->mean_rounded_down + (ch_par->mean_rounded_down - normal_distr_random_int);
    if(hh_size1 > hh_size2){ // Add extras to bigger hh size.
        hh_size1 += how_many_extras_ch(ch_par);
    }else{
        hh_size2 += how_many_extras_ch(ch_par);
    }
    
    auto home1 = world.add_location(epi::LocationType::Home);
    for ( auto i = 0; i < hh_size1; i++) {
        world.add_person(home1, determine_infection_state(inf_state_par), pick_age_from_pool(&ch_par->age_pool) );
    }
    
    auto home2 = world.add_location(epi::LocationType::Home);
    for ( auto i = 0; i < hh_size2; i++) {
        world.add_person(home2, determine_infection_state(inf_state_par), pick_age_from_pool(&ch_par->age_pool));
    }
}

void add_ch(epi::World& world, communityHhParameters *ch_par, infectionStateParameters *inf_state_par){
    int number_of_pairs_to_calculate = ch_par->number_of_households / 2 ; // This should always be dividable by 2.
    for (auto i = 0; i < number_of_pairs_to_calculate; i++) {
        add_pair_ch(world, ch_par, inf_state_par);
    }
}

void setup_specific_ch(communityHhParameters *ch_par, std::vector<int> age_pool_of_people, int number_of_households, double relative_deviation_from_mean_in_percent){
    ch_par->age_pool.pool_ages = age_pool_of_people;
    ch_par->number_of_people = std::accumulate(age_pool_of_people.begin(), age_pool_of_people.end(),0) ;
    ch_par->age_pool.total_pool_size =   ch_par->number_of_people;
    ch_par->number_of_households = number_of_households;
    
    ch_par->mean_rounded_down = ch_par->number_of_people/number_of_households;
    ch_par->number_of_extra_persons_left = ch_par->number_of_people - (ch_par->mean_rounded_down*number_of_households);
    ch_par->index_extra_persons = 0;
    
    ch_par->deviation_from_mean = (int) ((double) ch_par->mean_rounded_down * relative_deviation_from_mean_in_percent) ;
}

void create_community_households(epi::World& world, demograficParameters *demo_par, infectionStateParameters *inf_state_par){
    // Create the 4 different community households.
    // Retirement
    communityHhParameters ch_par_retirement;
    double retirement_deviation = 0.1;
    setup_specific_ch(&ch_par_retirement,demo_par->age_distr_comm_hh_retirement, demo_par->number_retirement_hh, retirement_deviation);
    add_ch(world, &ch_par_retirement, inf_state_par);
    // Disabled
    communityHhParameters ch_par_disabled;
    double disabled_deviation = 0.1;
    setup_specific_ch(&ch_par_disabled, demo_par->age_distr_comm_hh_disabled, demo_par->number_disabled_hh, disabled_deviation);
    add_ch(world, &ch_par_disabled, inf_state_par);
    // Refugee
    communityHhParameters ch_par_refugee;
    double refugee_deviation = 0.1;
    setup_specific_ch(&ch_par_refugee, demo_par->age_distr_comm_hh_refugee, demo_par->number_refugee_hh, refugee_deviation);
    add_ch(world, &ch_par_refugee, inf_state_par);
    // Others
    communityHhParameters ch_par_others;
    double others_deviation = 0.5;
    setup_specific_ch(&ch_par_others, demo_par->age_distr_comm_hh_other, demo_par->number_other_hh, others_deviation);
    add_ch(world, &ch_par_others, inf_state_par);
}

// This function creates the world from statistical data saved in demo_par
void create_world_from_statistical_data(epi::World& world, infectionStateParameters *inf_state_par){
    demograficParameters demo_par;
    setup_parameters(&demo_par);
    create_private_households(world, &demo_par, inf_state_par);
    create_community_households(world, &demo_par, inf_state_par);
}

/**
 * Add locations to the world and assign locations to the people.
 */
void create_assign_locations (epi::World& world){
    // Add one social event with 100 effective contacts.
    // Effective contacs limit the number of people that a person can infect while being at this location.
    auto event  = world.add_location(epi::LocationType::SocialEvent);
    world.get_individualized_location(event).get_infection_parameters().set<epi::EffectiveContacts>(100);
    
    // Add hospital and ICU with 5 effective contacs.
    auto hospital = world.add_location(epi::LocationType::Hospital);
    world.get_individualized_location(hospital).get_infection_parameters().set<epi::EffectiveContacts>(5);
    auto icu = world.add_location(epi::LocationType::ICU);
    world.get_individualized_location(icu).get_infection_parameters().set<epi::EffectiveContacts>(5);
    
    // Add schools, workplaces and shops.
    // At every school there are 600 students. The effective contacs are 40.
    // At every workplace work 100 people (needs to be varified), effective contacts are 40.
    // Add one supermarked per 15.000 people, effective constacts are assumed to be 20.
    auto shop  = world.add_location(epi::LocationType::BasicsShop);
    world.get_individualized_location(shop).get_infection_parameters().set<epi::EffectiveContacts>(20);
    
    auto school  = world.add_location(epi::LocationType::School);
    world.get_individualized_location(school).get_infection_parameters().set<epi::EffectiveContacts>(40);
    
    auto work  = world.add_location(epi::LocationType::Work);
    world.get_individualized_location(work).get_infection_parameters().set<epi::EffectiveContacts>(40);
    int counter_school = 0;
    int counter_work = 0;
    int counter_shop = 0;
    //Assign locations to the people
    auto persons = world.get_persons();
    for (auto& person: persons){
        //assign shop and event
        person.set_assigned_location(event);
        person.set_assigned_location(shop);
        counter_shop++;
        //assign hospital and ICU
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
        //assign work/school to people depending on their age
        if (person.get_age() == epi::AbmAgeGroup::Age5to14){
            person.set_assigned_location(school);
            counter_school++;
        }
        if (person.get_age() == epi::AbmAgeGroup::Age15to34 || person.get_age() == epi::AbmAgeGroup::Age35to59){
            person.set_assigned_location(work);
            counter_work++;
        }
        //add new school/work/shop if needed
        if (counter_school == 600){
            counter_school = 0;
            school  = world.add_location(epi::LocationType::School);
            world.get_individualized_location(school).get_infection_parameters().set<epi::EffectiveContacts>(40);
        }
        if (counter_work == 100){
            counter_work = 0;
            work  = world.add_location(epi::LocationType::Work);
            world.get_individualized_location(work).get_infection_parameters().set<epi::EffectiveContacts>(40);
        }
        if (counter_shop == 15000){
            counter_shop = 0;
            shop  = world.add_location(epi::LocationType::BasicsShop);
            world.get_individualized_location(shop).get_infection_parameters().set<epi::EffectiveContacts>(20);
        }
    }
}

int main(){
    // Setup infection state parameters
    infectionStateParameters inf_state_par;
    setup_infection_state_parameters(&inf_state_par);
    
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
    
    // Setup demografic parameters from statistical data and create the world from this data.
    auto world = epi::World(abm_params);
    create_world_from_statistical_data(world, &inf_state_par);
    
    // Add locations and assign locations to the people.
    create_assign_locations(world);
    
    auto t0 = epi::TimePoint(0);
    auto t_lockdown = epi::TimePoint(0)+epi::days(20);
    auto tmax = epi::TimePoint(0) + epi::days(60);
    
    // During the lockdown, 60% of people work from home and schools are closed for 90% of students.
    // Social events are very rare.
    epi::set_home_office(t_lockdown, 0.6, world.get_migration_parameters());
    epi::set_school_closure(t_lockdown, 0.9, world.get_migration_parameters());
    epi::close_social_events(t_lockdown, 0.9, world.get_migration_parameters());
    auto sim  = epi::AbmSimulation(t0, std::move(world));
    
    sim.advance(tmax);
    
    // The results are saved in a table with 9 rows.
    // The first row is t = time, the others correspond to the number of people with a certain infection state at this time:
    // S = Susceptible, E = Exposed, C= Carrier, I_d = Infected_Detected, I_u = Infected_Undetected, R_C = Recovered_Carrier, R_I = Recovered_Infected, D = Dead
    // E.g. the following gnuplot skrips plots detected infections and deaths.
    // plot "abm.txt" using 1:5 with lines title "infected (detected)", "abm.txt" using 1:9 with lines title "dead"
    // set xlabel "days"
    // set ylabel "number of people"
    // set title "ABM Example"
    // set output "abm.png"
    // set terminal png
    // replot
    auto f_abm = fopen("abm.txt", "w");
    fprintf(f_abm, "# t S E C I_d I_u R_C R_I D\n");
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









