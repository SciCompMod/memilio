//
//  world_test.cpp
//  
//
//  Created by Sascha Korf on 04.08.21.
//

#include "epidemiology/abm/abm.h"



/**
 * Determine the infection state of a person at the beginning of the simulation.
 * The infection states are chosen randomly. They are distributed according to the probabilites set in the example.
 * @return random infection state
 */
epi::InfectionState determine_infection_state(double exposed, double infected, double carrier, double recovered){
    double susceptible = 1 - exposed - infected - carrier - recovered;
    std::vector<double> weights = {susceptible, exposed, carrier, infected/2, infected/2, recovered/2, recovered/2};
    uint32_t state = epi::DiscreteDistribution<size_t>::get_instance()(weights);
    return (epi::InfectionState) state;
}

struct poolToChooseFrom{
    
    // This is a pool to chose from
    
    std::vector<int> pool_ages;
    
    int total_pool_size;
};

struct privateHhPools{
    poolToChooseFrom children_pool, parents_pool, rest_pool;
    int index_for_extras_in_5plus_households;
};
/*
 * Struct for defining the demografic parameters important for the building of world with these parameters.
 */
struct demograficParameters
{
    /* We will divide the parameters in safe and unsafe parameters
     * Safe: parameters we can safely get from reliable sources e.g. total population or the absolute number of people between 20-35
     * Unsafe: these parameters aren't as reliable. e.g. the age distribution of people in a 3 people Households
     * Other: we have additional parameters as avarageHH size of a 5 plus HH
     */
    
    /* Furthermore we will have the following dividing:
     * We will divide the general population in: 1. private Households 2. Community Homes
     * 1. We will divide the private Households (HH) into 1A 1 Person HH, 1B 2 person HH, 1C 3 Person HH, 1D 4 Person HH
     * 2. We will divide the community HH into 2A (Retirement Homes) 2B (Disabled Homes), 2C (Refugee Homes), 2D (Other Homes)
     */
    
    // -- Safe Parameters -- //
    int total_population; // Total population
    
    // people in different hh
    int people_in_priv_hh, people_in_comm_hh;
    int people_in_retirement_hh, people_in_disabled_hh, people_in_refugee_hh, people_in_other_hh;
    std::vector<int> people_in_priv_hh_each_size;
    
    
    // Age distributions
    //TODO: maybe use eigenvec or eigen::vector or eigen::matrix, and customIndexArray
    std::vector<int> age_distr_priv_hh, age_distr_comm_hh; // General distribution in first dividing above
    std::vector<int> age_distr_priv_hh_size1; // Only for Household of size 1 the ageDistr is safe
    std::vector<int> age_distr_comm_hh_retirement, age_distr_comm_hh_disabled, age_distr_comm_hh_refugee, age_distr_comm_hh_other; // These are kind of safe
    
    // Number of Household's of each Type
    int number_priv_hh, number_comm_hh;
    std::vector<int> number_priv_hh_each_size;
    int number_retirement_hh, number_disabled_hh, number_refugee_hh, number_other_hh;
    
    //Other parameters to be set
    std::vector<int> full_families_hh_size, half_families_hh_size, rest_families_hh_size;
};

struct communityHhParameters{
  
    poolToChooseFrom age_pool;
    int number_of_people, number_of_households;
    int deviation_from_mean;
    int mean_rounded_down;
    int number_of_extra_persons, index_extra_persons;
    
    
};
/*
 * Function that takes a vector of integers and returns a vector which has the same length. This vector has the percentage of the
 * whole vector sum in each row.
 * @param absoluteValueVec Vector with the absolute Values
 * @return  Vector with percent
 */
std::vector<double> get_percentage_of_sum(std::vector<int> absoluteValueVec){
    //TODO: lookup in utility
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

void setup_parameters(demograficParameters *demo_par){
    
    
    /*
     * Function that sets up the safe parameters
     */
    
    
    
    // Set general parameters
    // Set total population
    demo_par->total_population = 81848;
    
    
    // Set number of people in private / community HH
    

    demo_par->people_in_priv_hh = 80615;
    demo_par->people_in_comm_hh = 1233;
    
    
    
    // Set quotient of People in Community HH
    
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
    
    
    demo_par->number_retirement_hh=16;
    demo_par->number_disabled_hh=8;
    demo_par->number_refugee_hh=2;
    demo_par->number_other_hh=2;
    
    demo_par->number_priv_hh = std::accumulate(demo_par->number_priv_hh_each_size.begin(), demo_par->number_priv_hh_each_size.end(), 0);
    demo_par->number_comm_hh = demo_par->number_retirement_hh+demo_par->number_disabled_hh+demo_par->number_refugee_hh+demo_par->number_other_hh;
    
    
    //Number of full/half/rest families
    
    
    demo_par->full_families_hh_size.resize(5);
    demo_par->full_families_hh_size.at(0) = 0;
    demo_par->full_families_hh_size.at(1) = 11850;
    demo_par->full_families_hh_size.at(2) = 4155;
    demo_par->full_families_hh_size.at(3) = 3551;
    demo_par->full_families_hh_size.at(4) = 1245;
    
    demo_par->half_families_hh_size.resize(5);
    demo_par->half_families_hh_size.at(0) = 0;
    demo_par->half_families_hh_size.at(1) = 1765;
    demo_par->half_families_hh_size.at(2) = 662;
    demo_par->half_families_hh_size.at(3) = 110;
    demo_par->half_families_hh_size.at(4) = 80;
    
    demo_par->rest_families_hh_size.resize(5);
    demo_par->rest_families_hh_size.at(0) = 0;
    demo_par->rest_families_hh_size.at(1) = 166;
    demo_par->rest_families_hh_size.at(2) = 175;
    demo_par->rest_families_hh_size.at(3) = 122;
    demo_par->rest_families_hh_size.at(4) = 109;
    
    
    
}

void substract_from_pool(poolToChooseFrom *pool, int count, int age_group){
        
    pool->pool_ages.at(age_group) -= count;
    pool->total_pool_size -= count;
    
}

epi::AbmAgeGroup pick_next_available(privateHhPools *hh_pool, bool parents, bool children){
    
    int index_sufficient_agegroup=0;
    int age_group;
    
    if(children){
        for (auto i = 0; i<hh_pool->children_pool.pool_ages.size(); i++) {
            if( hh_pool->parents_pool.pool_ages.at(i) !=0 ){
                age_group = i;
                substract_from_pool(&hh_pool->parents_pool,1,age_group);
                return (epi::AbmAgeGroup) age_group;
            }
            
            if (hh_pool->rest_pool.pool_ages.at(i) !=0 ) {
                age_group = i;
                substract_from_pool(&hh_pool->rest_pool,1,age_group);
                return (epi::AbmAgeGroup) age_group;
            }
        }
        
    }
    
    if(parents){
        for (auto i = 2; i<hh_pool->parents_pool.pool_ages.size(); i++) {
            if( hh_pool->children_pool.pool_ages.at(i) !=0 ){
                age_group = i;
                substract_from_pool(&hh_pool->children_pool,1,age_group);
                return (epi::AbmAgeGroup) age_group;
            }
            
            if (hh_pool->rest_pool.pool_ages.at(i) !=0 ) {
                age_group = i;
                substract_from_pool(&hh_pool->rest_pool,1,age_group);
                return (epi::AbmAgeGroup) age_group;
            }
        }
        
    }
    
    
    
}

void setup_pool_private_households(demograficParameters *demo_par, privateHhPools *hh_pool){
    
    int num_househ; // Number of persons is equal to number of HH
    
    int rest_age_group_15to34 = demo_par->age_distr_priv_hh_size1.at(2);
    int rest_age_group_35to59 = demo_par->age_distr_priv_hh_size1.at(3);
    
    // Setup the pool from which the families and rest is taken from
    hh_pool->children_pool.pool_ages.resize(6,0);
    hh_pool->children_pool.pool_ages.at(0)=demo_par->age_distr_priv_hh.at(0);
    hh_pool->children_pool.pool_ages.at(1)=demo_par->age_distr_priv_hh.at(1);
    hh_pool->children_pool.pool_ages.at(2)=0.4*(demo_par->age_distr_priv_hh.at(2)-rest_age_group_15to34);
    
    hh_pool->children_pool.total_pool_size = std::accumulate( hh_pool->children_pool.pool_ages.begin(),  hh_pool->children_pool.pool_ages.end(), 0);
    
    // Pool in each age
    hh_pool->parents_pool.pool_ages.resize(6,0);
    hh_pool->parents_pool.pool_ages.at(2)=0.6*(demo_par->age_distr_priv_hh.at(2)-rest_age_group_15to34);
    hh_pool->parents_pool.pool_ages.at(3)=demo_par->age_distr_priv_hh.at(3)-rest_age_group_35to59;
    
    hh_pool->parents_pool.total_pool_size = std::accumulate( hh_pool->parents_pool.pool_ages.begin(),  hh_pool->parents_pool.pool_ages.end(), 0);
    
    hh_pool->rest_pool.pool_ages.resize(6,0);
    for (auto i = 0; i <  hh_pool->rest_pool.pool_ages.size(); i++) {
        hh_pool->rest_pool.pool_ages.at(i)=demo_par->age_distr_priv_hh.at(i)- hh_pool->children_pool.pool_ages.at(i)- hh_pool->parents_pool.pool_ages.at(i);
    }
    hh_pool->rest_pool.total_pool_size = std::accumulate( hh_pool->rest_pool.pool_ages.begin(),  hh_pool->rest_pool.pool_ages.end(), 0);
    
    hh_pool->index_for_extras_in_5plus_households = 0;
}

int determine_age_group_number_of_agegroups(std::vector<int> age_groups){
    
    std::vector<double> weights = get_percentage_of_sum(age_groups);
    uint32_t age_group = epi::DiscreteDistribution<size_t>::get_instance()(weights);
    return age_group;
}

epi::AbmAgeGroup pick_age_from_pool(poolToChooseFrom *pool){
    int age_group = determine_age_group_number_of_agegroups(pool->pool_ages);
    substract_from_pool(pool,1,age_group);
    return (epi::AbmAgeGroup) age_group;
}

epi::AbmAgeGroup pick_age_group_from_ph_pool(privateHhPools *hh_pool,bool child, bool parent){
    
    int age_group;
    
    if(child){
        if(hh_pool->children_pool.total_pool_size>0){
            age_group = determine_age_group_number_of_agegroups(hh_pool->children_pool.pool_ages );
            substract_from_pool(&hh_pool->children_pool,1,age_group);
            return (epi::AbmAgeGroup) age_group;
        }else{
            return pick_next_available(hh_pool, parent, child);
        }
    }else if (parent){
        if(hh_pool->parents_pool.total_pool_size>0){
            age_group = determine_age_group_number_of_agegroups(hh_pool->parents_pool.pool_ages );
            substract_from_pool(&hh_pool->parents_pool,1,age_group);
            return (epi::AbmAgeGroup) age_group;
        }else{
            return pick_next_available(hh_pool, parent, child);
        }
    }else{
        if(hh_pool->rest_pool.total_pool_size>0){
            age_group = determine_age_group_number_of_agegroups(hh_pool->rest_pool.pool_ages );
            substract_from_pool(&hh_pool->rest_pool,1,age_group);
            return (epi::AbmAgeGroup) age_group;
        }else{
            printf("This shouldn't happen.");
        }
    }
}

void update_pools_after_family_distribution(privateHhPools *hh_pool){
    
    for (auto i = 0; i < hh_pool->rest_pool.pool_ages.size(); i++) {
        if(hh_pool->children_pool.pool_ages.at(i) > 0){
            hh_pool->rest_pool.pool_ages.at(i) += hh_pool->children_pool.pool_ages.at(i);
            hh_pool->children_pool.total_pool_size -=  hh_pool->children_pool.pool_ages.at(i);
            hh_pool->children_pool.pool_ages.at(i) = 0;
            
        }
        if(hh_pool->parents_pool.pool_ages.at(i) > 0){
            hh_pool->rest_pool.pool_ages.at(i) += hh_pool->parents_pool.pool_ages.at(i);
            hh_pool->parents_pool.total_pool_size -=  hh_pool->parents_pool.pool_ages.at(i);
            hh_pool->parents_pool.pool_ages.at(i) = 0;
        }
        
    }
}

int how_many_extras(privateHhPools *hh_pool, demograficParameters *demo_par){
    
    int number_of_5plus_hh_to_compute = demo_par->rest_families_hh_size.at(4) - hh_pool->index_for_extras_in_5plus_households-1;
    hh_pool->index_for_extras_in_5plus_households++;
    
    int remaining_extra_people = hh_pool->rest_pool.total_pool_size - number_of_5plus_hh_to_compute * 5;
    int extras;
    
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

int how_many_extras_community( int hh_to_compute, int extras_left){
    
    
    
    
}



void add_community_pair_household_to_world(epi::World& world, communityHhParameters *chh_pool, double exposed_pct, double infected_pct, double carrier_pct, double recovered_pct){
    // This functions adds the number of households to the worlds. In the following way:
    // We need to distrinute x people to y households. The avarage rounded down is in chh_pool->mean_rounded_down.
    // We add 2 Households at the same time. First we calculate a random number z which is distributed with a standart distribution deviation of chh_pool->deviation_from_mean and add mean_rounded_down-z and mean_rounded_down+z+extra.
    // Extra are the remaining people which need to be added after the mean is rounded down and is calculated in the function how many extras.
    std::normal_distribution<double> distribution(chh_pool->mean_rounded_down,chh_pool->deviation_from_mean);
    std::default_random_engine generator;
    int normal_distr_random_int = (int) distribution(generator);
    if(normal_distr_random_int > chh_pool->mean_rounded_down){
        normal_distr_random_int = chh_pool->mean_rounded_down - 1;
    }
    
    int hh_size1 = normal_distr_random_int;
    int hh_size2 =  chh_pool->mean_rounded_down + (chh_pool->mean_rounded_down - normal_distr_random_int);
    //hh_size1 += to be implemented TODO
    //TODO: see above
    
    auto home1 = world.add_location(epi::LocationType::Home);
    for ( auto i = 0; i < hh_size1; i++) {
        world.add_person(home1, determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct), pick_age_from_pool(&chh_pool->age_pool) );
    }
    
    auto home2 = world.add_location(epi::LocationType::Home);
    for ( auto i = 0; i < hh_size1; i++) {
        world.add_person(home2, determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct), pick_age_from_pool(&chh_pool->age_pool));
    }
}

void add_community_households(epi::World& world, communityHhParameters *chh_pool, double exposed_pct, double infected_pct, double carrier_pct, double recovered_pct){
    int number_of_pairs_to_calculate = chh_pool->number_of_households / 2 ; // This should always be dividable by 2.
    for (auto i = 0; i < number_of_pairs_to_calculate; i++) {
        add_community_pair_household_to_world(world, chh_pool, exposed_pct, infected_pct, carrier_pct, recovered_pct);
    }
}

void create_specific_community_household(communityHhParameters *chh_pool, std::vector<int> age_pool_of_people, int number_of_households, double relative_deviation_from_mean_in_percent){
    chh_pool->age_pool.pool_ages = age_pool_of_people;
    chh_pool->number_of_people = std::accumulate(age_pool_of_people.begin(), age_pool_of_people.end(),0) ;
    chh_pool->age_pool.total_pool_size =   chh_pool->number_of_people;
    chh_pool->number_of_households = number_of_households;
    
    chh_pool->mean_rounded_down = chh_pool->number_of_people/number_of_households;
    chh_pool->number_of_extra_persons = chh_pool->number_of_people - (chh_pool->mean_rounded_down*number_of_households);
    chh_pool->index_extra_persons = 0;
    
    chh_pool->deviation_from_mean = (int) ((double) chh_pool->mean_rounded_down * relative_deviation_from_mean_in_percent) ;
}

void create_community_households(epi::World& world, demograficParameters *demo_par, double exposed_pct, double infected_pct, double carrier_pct, double recovered_pct){
    //Create the 4 different community households
    //Retirement
    communityHhParameters chh_pool_retirement;
    create_specific_community_household(&chh_pool_retirement,demo_par->age_distr_comm_hh_retirement, demo_par->number_retirement_hh, 0.1);
    add_community_households(world, &chh_pool_retirement, exposed_pct, infected_pct, carrier_pct, recovered_pct);
    //Disabled
    communityHhParameters chh_pool_disabled;
    create_specific_community_household(&chh_pool_disabled, demo_par->age_distr_comm_hh_disabled, demo_par->number_disabled_hh, 0.1);
    add_community_households(world, &chh_pool_disabled, exposed_pct, infected_pct, carrier_pct, recovered_pct);
    //Refugee
    communityHhParameters chh_pool_refugee;
    create_specific_community_household(&chh_pool_refugee, demo_par->age_distr_comm_hh_refugee, demo_par->number_refugee_hh, 0.1);
    add_community_households(world, &chh_pool_refugee, exposed_pct, infected_pct, carrier_pct, recovered_pct);
    //Others
    communityHhParameters chh_pool_others;
    create_specific_community_household(&chh_pool_others, demo_par->age_distr_comm_hh_other, demo_par->number_other_hh, 0.5);
    add_community_households(world, &chh_pool_others, exposed_pct, infected_pct, carrier_pct, recovered_pct);
}

void create_full_families_of_private_households(epi::World& world, demograficParameters *demo_par, privateHhPools *hh_pool, double exposed_pct, double infected_pct, double carrier_pct, double recovered_pct, int household_size){
    
    
    //Full families mean two parents and children, e.g. household_size 4 equals 2 parents and 2 children,  household_size 5 means 2 parents and 3 children
    
    if(household_size <2 || household_size > 5 )
        printf("This is not an acceptable housheold size");
    
    //first we create a household and two parents
    
    //Household
    auto home = world.add_location(epi::LocationType::Home);
    
    //Two parents
    world.add_person(home, determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct), pick_age_group_from_ph_pool(hh_pool, false, true));
    world.add_person(home, determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct), pick_age_group_from_ph_pool(hh_pool, false, true));
    
    // The rest are children
    
    int number_of_children = household_size-2;
    
    for (auto i = 0 ; i < number_of_children; i++) {
        auto child = pick_age_group_from_ph_pool(hh_pool, true, false);
        world.add_person(home, determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct), child);
    }
    
    
    
    
}

void create_half_families_of_private_households(epi::World& world, demograficParameters *demo_par, privateHhPools *hh_pool, double exposed_pct, double infected_pct, double carrier_pct, double recovered_pct, int household_size){
    
    
    // Half families mean one parent and children, e.g. household_size 4 equals 1 parents and 3 children,  household_size 5 means 1 parents and 4 children
    
    if(household_size <2 || household_size > 5 )
        printf("This is not an acceptable househeold size");
    
    //first we create a household and two parents
    
    //Household
    auto home = world.add_location(epi::LocationType::Home);
    
    //One  parents
    auto parent1 = pick_age_group_from_ph_pool(hh_pool, true, false);
    auto& p = world.add_person(home, determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct), parent1);
    
    // The rest are children
    
    int number_of_children = household_size-1;
    
    for (auto i = 0 ; i < number_of_children; i++) {
        auto child = pick_age_group_from_ph_pool(hh_pool, false, true);
        world.add_person(home, determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct), child);
    }
    
    
    
    
}

void create_rest_of_private_households(epi::World& world, demograficParameters *demo_par, privateHhPools *hh_pool, double exposed_pct, double infected_pct, double carrier_pct, double recovered_pct, int household_size){
    
    if(household_size <2 || household_size > 5 )
        printf("This is not an acceptable househeold size");
    
    //first we create a household and two parents
    
    //Household
    
   
    auto home = world.add_location(epi::LocationType::Home);
    
    for (auto i = 0 ; i < household_size; i++) {
        auto person = pick_age_group_from_ph_pool(hh_pool, false, false);
        world.add_person(home, determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct), person);
    }
    
    if (household_size  == 5) {
        int extras = how_many_extras(hh_pool, demo_par);
        for (auto i = 0 ; i < extras; i++) {
            auto person = pick_age_group_from_ph_pool(hh_pool, false, false);
            world.add_person(home, determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct), person);
        }
    }



    
    
    
    
}

void create_more_person_households(epi::World& world, demograficParameters *demo_par, privateHhPools *hh_pool, double exposed_pct, double infected_pct, double carrier_pct, double recovered_pct){
    
    

    // Here we create the more than one person households
    
    // We go through this like this: First we create "families" in each of the households of x size
    // Then we create the rest with the rest of the pool
    
    // For each Household we cycle through the families and then fill the rest
    
    
    // Fill the families
    for (auto household_size = 2; household_size < 6; household_size++){
        for (int i = 0; i<demo_par->full_families_hh_size.at(household_size-1); i++) {
            create_full_families_of_private_households(world, demo_par, hh_pool, exposed_pct, infected_pct, carrier_pct,recovered_pct, household_size);
        }
      
      
    }
    
    // Fill the half families
    for (auto household_size = 2; household_size < 6; household_size++){
        for (int i = 0; i<demo_par->half_families_hh_size.at(household_size-1); i++) {
            create_half_families_of_private_households(world, demo_par, hh_pool, exposed_pct, infected_pct, carrier_pct,recovered_pct, household_size);
        }
    }
    
    
    
    update_pools_after_family_distribution(hh_pool);
    // Fill the rest of the households
    
    for (auto household_size = 2; household_size < 6; household_size++){
        for (int i = 0; i<demo_par->rest_families_hh_size.at(household_size-1); i++) {
            create_rest_of_private_households(world, demo_par, hh_pool, exposed_pct, infected_pct, carrier_pct,recovered_pct, household_size);
        }
    }
    
    
    
    
    
}

void create_one_person_household(epi::World& world, demograficParameters *demo_par, privateHhPools *hh_pool, double exposed_pct, double infected_pct, double carrier_pct, double recovered_pct){
    
    //Create all the 1-person Household because the data is clear on this
    
    int num_people_oneP_HH = demo_par->number_priv_hh_each_size.at(0); // Number of persons is equal to number of HH
    int number_people_this_agegroup = 0;
    int age_groups = 6;
    
    epi::AbmAgeGroup age_group;
    
    
    for(auto i = 0; i < age_groups; i++){
        
        number_people_this_agegroup = demo_par->age_distr_priv_hh_size1.at(i);
        
        if(number_people_this_agegroup > 0){
            age_group = (epi::AbmAgeGroup) i;
            
            
            
        
            for (auto j = 0; j < number_people_this_agegroup; j++) {
                auto home = world.add_location(epi::LocationType::Home);
                auto& p = world.add_person(home, determine_infection_state(exposed_pct, infected_pct, carrier_pct, recovered_pct), age_group);
                substract_from_pool(&hh_pool->rest_pool, 1, i);
            }
            
        }
        
        
       
    }
    
    
}

void create_world(epi::World& world, demograficParameters *demo_par, double exposed_pct, double infected_pct, double carrier_pct, double recovered_pct){
    privateHhPools hh_pool;
    setup_pool_private_households(demo_par, &hh_pool);
    create_one_person_household(world, demo_par, &hh_pool, exposed_pct, infected_pct, carrier_pct, recovered_pct);
    create_more_person_households(world, demo_par, &hh_pool, exposed_pct, infected_pct, carrier_pct, recovered_pct);
    create_community_households(world,demo_par, exposed_pct, infected_pct, carrier_pct, recovered_pct);
}

void setup_world(epi::World& world, demograficParameters *demo_par, double exposed_pct, double infected_pct, double carrier_pct, double recovered_pct){
    // Setup parameters
    setup_parameters(demo_par);
    // Create the world object
    create_world(world, demo_par, exposed_pct, infected_pct, carrier_pct, recovered_pct);
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
    
    
    double exposed_pct = 0.01, infected_pct = 0.008, carrier_pct = 0.005, recovered_pct = 0.001;
    
    
    //Set global infection parameters (similar to infection parameters in SECIR model) and initialize the world
    epi::GlobalInfectionParameters abm_params;
    abm_params.set<epi::IncubationPeriod>({{epi::AbmAgeGroup::Count}, 3.});
    abm_params.set<epi::SusceptibleToExposedByCarrier>({{epi::AbmAgeGroup::Count}, 0.02});
    abm_params.set<epi::SusceptibleToExposedByInfected>({{epi::AbmAgeGroup::Count}, 0.02});
    abm_params.set<epi::CarrierToInfected>({{epi::AbmAgeGroup::Count}, 0.15});
    abm_params.set<epi::CarrierToRecovered>({{epi::AbmAgeGroup::Count}, 0.15});
    abm_params.set<epi::InfectedToRecovered>({{epi::AbmAgeGroup::Count}, 0.2});
    abm_params.set<epi::InfectedToDead>({{epi::AbmAgeGroup::Count}, 0.02});
    abm_params.set<epi::RecoveredToSusceptible>({{epi::AbmAgeGroup::Count}, 0.});
    
    auto world    = epi::World(abm_params);
    
    demograficParameters demo_par;
    
    setup_world(world, &demo_par, exposed_pct, infected_pct, carrier_pct, recovered_pct);
    
    //Add locations and assign locations to the people.
    //create_assign_locations(world);

    auto t0   = epi::TimePoint(0);
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









