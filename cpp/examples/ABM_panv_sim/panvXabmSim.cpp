#include <fstream>
#include <vector>
#include <iostream>
#include <chrono>
#include <string>
#include <map>
#include "abm/abm.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/uncertain_value.h"
#include "boost/filesystem.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"
#include "abm/vaccine.h"
#include "abm/common_abm_loggers.h"
#include "memilio/utils/miompi.h"
#include "memilio/io/binary_serializer.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "abm/household.h"
#include "memilio/utils/random_number_generator.h"
#include "abm/common_abm_loggers.h"
#include "boost/filesystem.hpp"

// Number of age groups and their definitions
size_t num_age_groups        = 6;
const auto age_group_0_to_4   = mio::AgeGroup(0);
const auto age_group_5_to_14  = mio::AgeGroup(1);
const auto age_group_15_to_34 = mio::AgeGroup(2);
const auto age_group_35_to_59 = mio::AgeGroup(3);
const auto age_group_60_to_79 = mio::AgeGroup(4);
const auto age_group_80_plus  = mio::AgeGroup(5);

struct LogInfectionStatePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    /** 
     * @brief Log the TimeSeries of the number of Person%s in an #InfectionState.
     * @param[in] sim The simulation of the abm.
     * @return A pair of the TimePoint and the TimeSeries of the number of Person%s in an #InfectionState.
     */
    static Type log(const mio::abm::Simulation& sim)
    {

        Eigen::VectorXd sum = Eigen::VectorXd::Zero(
            Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups()));
        auto curr_time     = sim.get_time();
        const auto persons = sim.get_world().get_persons();

        for (auto i = size_t(0); i < persons.size(); ++i) {
            auto& p = persons[i];
            if (p.get_should_be_logged()) {
                auto index = (((size_t)(mio::abm::InfectionState::Count)) * ((uint32_t)p.get_age().get())) +
                             ((uint32_t)p.get_infection_state(curr_time));
                // PRAGMA_OMP(atomic)
                sum[index] += 1;
            }
        }
        return std::make_pair(curr_time, sum);
    }
};
struct LogInfectionPerLocationTypePerAgeGroup : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    /** 
     * @brief Log the TimeSeries of the number of Person%s in an #InfectionState.
     * @param[in] sim The simulation of the abm.
     * @return A pair of the TimePoint and the TimeSeries of the number of Person%s in an #InfectionState.
     */
    static Type log(const mio::abm::Simulation& sim)
    {

        Eigen::VectorXd sum = Eigen::VectorXd::Zero(
            Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups()));
        auto curr_time     = sim.get_time();
        auto prev_time     = sim.get_prev_time();
        const auto persons = sim.get_world().get_persons();

        // PRAGMA_OMP(parallel for)
        for (auto i = size_t(0); i < persons.size(); ++i) {
            auto& p = persons[i];
            if (p.get_should_be_logged()) {
                // PRAGMA_OMP(atomic)
                if ((p.get_infection_state(prev_time) != mio::abm::InfectionState::Exposed) &&
                    (p.get_infection_state(curr_time) == mio::abm::InfectionState::Exposed)) {
                    auto index = (((size_t)(mio::abm::LocationType::Count)) * ((uint32_t)p.get_age().get())) +
                                 ((uint32_t)p.get_location().get_type());
                    sum[index] += 1;
                }
            }
        }
        return std::make_pair(curr_time, sum);
    }
};


/**
 * Helper function to convert log-normal parameters from mean and std to mu and sigma
 */
std::pair<double, double> get_my_and_sigma(std::pair<double, double> mean_and_std)
{
    auto mean    = mean_and_std.first;
    auto stddev  = mean_and_std.second;
    double my    = log(mean * mean / sqrt(mean * mean + stddev * stddev));
    double sigma = sqrt(log(1 + stddev * stddev / (mean * mean)));
    return {my, sigma};
}

/**
 * Sets the parameters for the disease model
 */
void set_parameters(mio::abm::Parameters& params)
{
    auto incubation_period_my_sigma          = get_my_and_sigma({4.5, 1.5});
    params.get<mio::abm::IncubationPeriod>() = {incubation_period_my_sigma.first, incubation_period_my_sigma.second};

    auto InfectedNoSymptoms_to_symptoms_my_sigma             = get_my_and_sigma({1.1, 0.9});
    params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>() = {InfectedNoSymptoms_to_symptoms_my_sigma.first,
                                                                InfectedNoSymptoms_to_symptoms_my_sigma.second};

    auto TimeInfectedNoSymptomsToRecovered_my_sigma           = get_my_and_sigma({8.0, 2.0});
    params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>() = {TimeInfectedNoSymptomsToRecovered_my_sigma.first,
                                                                 TimeInfectedNoSymptomsToRecovered_my_sigma.second};

    auto TimeInfectedSymptomsToSevere_my_sigma           = get_my_and_sigma({6.6, 4.9});
    params.get<mio::abm::TimeInfectedSymptomsToSevere>() = {TimeInfectedSymptomsToSevere_my_sigma.first,
                                                            TimeInfectedSymptomsToSevere_my_sigma.second};

    auto TimeInfectedSymptomsToRecovered_my_sigma           = get_my_and_sigma({8.0, 2.0});
    params.get<mio::abm::TimeInfectedSymptomsToRecovered>() = {TimeInfectedSymptomsToRecovered_my_sigma.first,
                                                               TimeInfectedSymptomsToRecovered_my_sigma.second};

    auto TimeInfectedSevereToCritical_my_sigma           = get_my_and_sigma({1.5, 2.0});
    params.get<mio::abm::TimeInfectedSevereToCritical>() = {TimeInfectedSevereToCritical_my_sigma.first,
                                                            TimeInfectedSevereToCritical_my_sigma.second};

    auto TimeInfectedSevereToRecovered_my_sigma           = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedSevereToRecovered>() = {TimeInfectedSevereToRecovered_my_sigma.first,
                                                             TimeInfectedSevereToRecovered_my_sigma.second};

    auto TimeInfectedCriticalToDead_my_sigma           = get_my_and_sigma({10.7, 4.8});
    params.get<mio::abm::TimeInfectedCriticalToDead>() = {TimeInfectedCriticalToDead_my_sigma.first,
                                                          TimeInfectedCriticalToDead_my_sigma.second};

    auto TimeInfectedCriticalToRecovered_my_sigma           = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedCriticalToRecovered>() = {TimeInfectedCriticalToRecovered_my_sigma.first,
                                                               TimeInfectedCriticalToRecovered_my_sigma.second};

    // Set percentage parameters
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]  = 0.75;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = 0.75;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.8;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.8;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.8;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]  = 0.8;

    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.01;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]  = 0.01;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.02;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.07;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.3;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]  = 0.6;

    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.1;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]  = 0.1;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.1;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.17;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.6;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]  = 0.8;

    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.055;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]  = 0.055;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.14;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.28;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.55;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]  = 0.75;

    // Set infection parameters
    params.get<mio::abm::InfectionRateFromViralShed>()[{mio::abm::VirusVariant::Wildtype}] = 1;
    params.get<mio::abm::AerosolTransmissionRates>() = 0.0;
    
    // Set contact parameter (people of same age group meet more often)
    // params.get<mio::abm::AgeGroupGotoSocialEvent>() = true;
}

void set_local_parameters(mio::abm::World& world)
{
    const int n_age_groups = (int)world.parameters.get_num_groups();

    // setting this up in matrix-form would be much nicer,
    // but we somehow can't construct Eigen object with initializer lists
    /* baseline_home
        0.4413 0.4504 1.2383 0.8033 0.0494 0.0017
        0.0485 0.7616 0.6532 1.1614 0.0256 0.0013
        0.1800 0.1795 0.8806 0.6413 0.0429 0.0032
        0.0495 0.2639 0.5189 0.8277 0.0679 0.0014
        0.0087 0.0394 0.1417 0.3834 0.7064 0.0447
        0.0292 0.0648 0.1248 0.4179 0.3497 0.1544
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_home(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_home[{age_group_0_to_4, age_group_0_to_4}]     = 0.4413;
    contacts_home[{age_group_0_to_4, age_group_5_to_14}]    = 0.0504;
    contacts_home[{age_group_0_to_4, age_group_15_to_34}]   = 1.2383;
    contacts_home[{age_group_0_to_4, age_group_35_to_59}]   = 0.8033;
    contacts_home[{age_group_0_to_4, age_group_60_to_79}]   = 0.0494;
    contacts_home[{age_group_0_to_4, age_group_80_plus}]    = 0.0017;
    contacts_home[{age_group_5_to_14, age_group_0_to_4}]    = 0.0485;
    contacts_home[{age_group_5_to_14, age_group_5_to_14}]   = 0.7616;
    contacts_home[{age_group_5_to_14, age_group_15_to_34}]  = 0.6532;
    contacts_home[{age_group_5_to_14, age_group_35_to_59}]  = 1.1614;
    contacts_home[{age_group_5_to_14, age_group_60_to_79}]  = 0.0256;
    contacts_home[{age_group_5_to_14, age_group_80_plus}]   = 0.0013;
    contacts_home[{age_group_15_to_34, age_group_0_to_4}]   = 0.1800;
    contacts_home[{age_group_15_to_34, age_group_5_to_14}]  = 0.1795;
    contacts_home[{age_group_15_to_34, age_group_15_to_34}] = 0.8806;
    contacts_home[{age_group_15_to_34, age_group_35_to_59}] = 0.6413;
    contacts_home[{age_group_15_to_34, age_group_60_to_79}] = 0.0429;
    contacts_home[{age_group_15_to_34, age_group_80_plus}]  = 0.0032;
    contacts_home[{age_group_35_to_59, age_group_0_to_4}]   = 0.0495;
    contacts_home[{age_group_35_to_59, age_group_5_to_14}]  = 0.2639;
    contacts_home[{age_group_35_to_59, age_group_15_to_34}] = 0.5189;
    contacts_home[{age_group_35_to_59, age_group_35_to_59}] = 0.8277;
    contacts_home[{age_group_35_to_59, age_group_60_to_79}] = 0.0679;
    contacts_home[{age_group_35_to_59, age_group_80_plus}]  = 0.0014;
    contacts_home[{age_group_60_to_79, age_group_0_to_4}]   = 0.0087;
    contacts_home[{age_group_60_to_79, age_group_5_to_14}]  = 0.0394;
    contacts_home[{age_group_60_to_79, age_group_15_to_34}] = 0.1417;
    contacts_home[{age_group_60_to_79, age_group_35_to_59}] = 0.3834;
    contacts_home[{age_group_60_to_79, age_group_60_to_79}] = 0.7064;
    contacts_home[{age_group_60_to_79, age_group_80_plus}]  = 0.0447;
    contacts_home[{age_group_80_plus, age_group_0_to_4}]    = 0.0292;
    contacts_home[{age_group_80_plus, age_group_5_to_14}]   = 0.0648;
    contacts_home[{age_group_80_plus, age_group_15_to_34}]  = 0.1248;
    contacts_home[{age_group_80_plus, age_group_35_to_59}]  = 0.4179;
    contacts_home[{age_group_80_plus, age_group_60_to_79}]  = 0.3497;
    contacts_home[{age_group_80_plus, age_group_80_plus}]   = 0.1544;

    /* baseline_school
        1.1165 0.2741 0.2235 0.1028 0.0007 0.0000
        0.1627 1.9412 0.2431 0.1780 0.0130 0.0000
        0.0148 0.1646 1.1266 0.0923 0.0074 0.0000
        0.0367 0.1843 0.3265 0.0502 0.0021 0.0005
        0.0004 0.0370 0.0115 0.0014 0.0039 0.0000
        0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_school(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_school[{age_group_0_to_4, age_group_0_to_4}]     = 1.1165;
    contacts_school[{age_group_0_to_4, age_group_5_to_14}]    = 0.2741;
    contacts_school[{age_group_0_to_4, age_group_15_to_34}]   = 0.2235;
    contacts_school[{age_group_0_to_4, age_group_35_to_59}]   = 0.1028;
    contacts_school[{age_group_0_to_4, age_group_60_to_79}]   = 0.0007;
    contacts_school[{age_group_0_to_4, age_group_80_plus}]    = 0.0000;
    contacts_school[{age_group_5_to_14, age_group_0_to_4}]    = 0.1627;
    contacts_school[{age_group_5_to_14, age_group_5_to_14}]   = 1.9412;
    contacts_school[{age_group_5_to_14, age_group_15_to_34}]  = 0.2431;
    contacts_school[{age_group_5_to_14, age_group_35_to_59}]  = 0.1780;
    contacts_school[{age_group_5_to_14, age_group_60_to_79}]  = 0.0130;
    contacts_school[{age_group_5_to_14, age_group_80_plus}]   = 0.0000;
    contacts_school[{age_group_15_to_34, age_group_0_to_4}]   = 0.0148;
    contacts_school[{age_group_15_to_34, age_group_5_to_14}]  = 0.1646;
    contacts_school[{age_group_15_to_34, age_group_15_to_34}] = 1.1266;
    contacts_school[{age_group_15_to_34, age_group_35_to_59}] = 0.0923;
    contacts_school[{age_group_15_to_34, age_group_60_to_79}] = 0.0074;
    contacts_school[{age_group_15_to_34, age_group_80_plus}]  = 0.0000;
    contacts_school[{age_group_35_to_59, age_group_0_to_4}]   = 0.0367;
    contacts_school[{age_group_35_to_59, age_group_5_to_14}]  = 0.1843;
    contacts_school[{age_group_35_to_59, age_group_15_to_34}] = 0.3265;
    contacts_school[{age_group_35_to_59, age_group_35_to_59}] = 0.0502;
    contacts_school[{age_group_35_to_59, age_group_60_to_79}] = 0.0021;
    contacts_school[{age_group_35_to_59, age_group_80_plus}]  = 0.0005;
    contacts_school[{age_group_60_to_79, age_group_0_to_4}]   = 0.0004;
    contacts_school[{age_group_60_to_79, age_group_5_to_14}]  = 0.0370;
    contacts_school[{age_group_60_to_79, age_group_15_to_34}] = 0.0115;
    contacts_school[{age_group_60_to_79, age_group_35_to_59}] = 0.0014;
    contacts_school[{age_group_60_to_79, age_group_60_to_79}] = 0.0039;
    contacts_school[{age_group_60_to_79, age_group_80_plus}]  = 0.0000;
    contacts_school[{age_group_80_plus, age_group_0_to_4}]    = 0.0000;
    contacts_school[{age_group_80_plus, age_group_5_to_14}]   = 0.0000;
    contacts_school[{age_group_80_plus, age_group_15_to_34}]  = 0.0000;
    contacts_school[{age_group_80_plus, age_group_35_to_59}]  = 0.0000;
    contacts_school[{age_group_80_plus, age_group_60_to_79}]  = 0.0000;
    contacts_school[{age_group_80_plus, age_group_80_plus}]   = 0.0000;

    /* baseline_work
        0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
        0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
        0.0000 0.0127 1.7570 1.6050 0.0133 0.0000
        0.0000 0.0020 1.0311 2.3166 0.0098 0.0000
        0.0000 0.0002 0.0194 0.0325 0.0003 0.0000
        0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_work(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_work[{age_group_0_to_4, age_group_0_to_4}]     = 0.0000;
    contacts_work[{age_group_0_to_4, age_group_5_to_14}]    = 0.0000;
    contacts_work[{age_group_0_to_4, age_group_15_to_34}]   = 0.0000;
    contacts_work[{age_group_0_to_4, age_group_35_to_59}]   = 0.0000;
    contacts_work[{age_group_0_to_4, age_group_60_to_79}]   = 0.0000;
    contacts_work[{age_group_0_to_4, age_group_80_plus}]    = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_0_to_4}]    = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_5_to_14}]   = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_15_to_34}]  = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_35_to_59}]  = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_60_to_79}]  = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_80_plus}]   = 0.0000;
    contacts_work[{age_group_15_to_34, age_group_0_to_4}]   = 0.0000;
    contacts_work[{age_group_15_to_34, age_group_5_to_14}]  = 0.0127;
    contacts_work[{age_group_15_to_34, age_group_15_to_34}] = 1.7570;
    contacts_work[{age_group_15_to_34, age_group_35_to_59}] = 1.6050;
    contacts_work[{age_group_15_to_34, age_group_60_to_79}] = 0.0133;
    contacts_work[{age_group_15_to_34, age_group_80_plus}]  = 0.0000;
    contacts_work[{age_group_35_to_59, age_group_0_to_4}]   = 0.0000;
    contacts_work[{age_group_35_to_59, age_group_5_to_14}]  = 0.0020;
    contacts_work[{age_group_35_to_59, age_group_15_to_34}] = 1.0311;
    contacts_work[{age_group_35_to_59, age_group_35_to_59}] = 2.3166;
    contacts_work[{age_group_35_to_59, age_group_60_to_79}] = 0.0098;
    contacts_work[{age_group_35_to_59, age_group_80_plus}]  = 0.0000;
    contacts_work[{age_group_60_to_79, age_group_0_to_4}]   = 0.0000;
    contacts_work[{age_group_60_to_79, age_group_5_to_14}]  = 0.0002;
    contacts_work[{age_group_60_to_79, age_group_15_to_34}] = 0.0194;
    contacts_work[{age_group_60_to_79, age_group_35_to_59}] = 0.0325;
    contacts_work[{age_group_60_to_79, age_group_60_to_79}] = 0.0003;
    contacts_work[{age_group_60_to_79, age_group_80_plus}]  = 0.0000;
    contacts_work[{age_group_80_plus, age_group_0_to_4}]    = 0.0000;
    contacts_work[{age_group_80_plus, age_group_5_to_14}]   = 0.0000;
    contacts_work[{age_group_80_plus, age_group_15_to_34}]  = 0.0000;
    contacts_work[{age_group_80_plus, age_group_35_to_59}]  = 0.0000;
    contacts_work[{age_group_80_plus, age_group_60_to_79}]  = 0.0000;
    contacts_work[{age_group_80_plus, age_group_80_plus}]   = 0.0000;

    /* baseline_other
        0.5170 0.3997 0.7957 0.9958 0.3239 0.0428
        0.0632 0.9121 0.3254 0.4731 0.2355 0.0148
        0.0336 0.1604 1.7529 0.8622 0.1440 0.0077
        0.0204 0.1444 0.5738 1.2127 0.3433 0.0178
        0.0371 0.0393 0.4171 0.9666 0.7495 0.0257
        0.0791 0.0800 0.3480 0.5588 0.2769 0.0180
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_other(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_other[{age_group_0_to_4, age_group_0_to_4}]     = 0.5170;
    contacts_other[{age_group_0_to_4, age_group_5_to_14}]    = 0.3997;
    contacts_other[{age_group_0_to_4, age_group_15_to_34}]   = 0.7957;
    contacts_other[{age_group_0_to_4, age_group_35_to_59}]   = 0.9958;
    contacts_other[{age_group_0_to_4, age_group_60_to_79}]   = 0.3239;
    contacts_other[{age_group_0_to_4, age_group_80_plus}]    = 0.0428;
    contacts_other[{age_group_5_to_14, age_group_0_to_4}]    = 0.0632;
    contacts_other[{age_group_5_to_14, age_group_5_to_14}]   = 0.9121;
    contacts_other[{age_group_5_to_14, age_group_15_to_34}]  = 0.3254;
    contacts_other[{age_group_5_to_14, age_group_35_to_59}]  = 0.4731;
    contacts_other[{age_group_5_to_14, age_group_60_to_79}]  = 0.2355;
    contacts_other[{age_group_5_to_14, age_group_80_plus}]   = 0.0148;
    contacts_other[{age_group_15_to_34, age_group_0_to_4}]   = 0.0336;
    contacts_other[{age_group_15_to_34, age_group_5_to_14}]  = 0.1604;
    contacts_other[{age_group_15_to_34, age_group_15_to_34}] = 1.7529;
    contacts_other[{age_group_15_to_34, age_group_35_to_59}] = 0.8622;
    contacts_other[{age_group_15_to_34, age_group_60_to_79}] = 0.1440;
    contacts_other[{age_group_15_to_34, age_group_80_plus}]  = 0.0077;
    contacts_other[{age_group_35_to_59, age_group_0_to_4}]   = 0.0204;
    contacts_other[{age_group_35_to_59, age_group_5_to_14}]  = 0.1444;
    contacts_other[{age_group_35_to_59, age_group_15_to_34}] = 0.5738;
    contacts_other[{age_group_35_to_59, age_group_35_to_59}] = 1.2127;
    contacts_other[{age_group_35_to_59, age_group_60_to_79}] = 0.3433;
    contacts_other[{age_group_35_to_59, age_group_80_plus}]  = 0.0178;
    contacts_other[{age_group_60_to_79, age_group_0_to_4}]   = 0.0371;
    contacts_other[{age_group_60_to_79, age_group_5_to_14}]  = 0.0393;
    contacts_other[{age_group_60_to_79, age_group_15_to_34}] = 0.4171;
    contacts_other[{age_group_60_to_79, age_group_35_to_59}] = 0.9666;
    contacts_other[{age_group_60_to_79, age_group_60_to_79}] = 0.7495;
    contacts_other[{age_group_60_to_79, age_group_80_plus}]  = 0.0257;
    contacts_other[{age_group_80_plus, age_group_0_to_4}]    = 0.0791;
    contacts_other[{age_group_80_plus, age_group_5_to_14}]   = 0.0800;
    contacts_other[{age_group_80_plus, age_group_15_to_34}]  = 0.3480;
    contacts_other[{age_group_80_plus, age_group_35_to_59}]  = 0.5588;
    contacts_other[{age_group_80_plus, age_group_60_to_79}]  = 0.2769;
    contacts_other[{age_group_80_plus, age_group_80_plus}]   = 0.0180;

    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_random(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 1.0);

    for (auto& loc : world.get_locations()) {
        switch (loc.get_type()) {
        case mio::abm::LocationType::Home:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_home;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1.6; //15 hours
            break;
        case mio::abm::LocationType::School:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_school;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 12.0; //2 hours
            break;
        case mio::abm::LocationType::Work:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_work;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 8.0; // 3 hours
            break;
        case mio::abm::LocationType::SocialEvent:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1.2;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 8.0; // 3 hours
            break;
        case mio::abm::LocationType::BasicsShop:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1.8;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 12.0; // 2 hours
            break;
        default:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_random;
            break;
        }
    }
}



/**
 * Read infection data from file and return a map of infected person IDs
 */
std::map<uint32_t, bool> read_infection_data(const std::string& filename)
{
    std::map<uint32_t, bool> infected_status;
    std::map<std::string, std::vector<uint32_t>> tables;
    
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return infected_status;
    }
    
    std::string line;
    // Skip header line
    std::getline(file, line);
    
    // Read each line
    while (std::getline(file, line)) {
        // Parse tab-separated values
        std::istringstream iss(line);
        uint32_t id;
        std::string table;
        int infected;
        
        if (iss >> id >> table >> infected) {
            infected_status[id] = (infected == 1);
            tables[table].push_back(id);
        }
    }
    
    file.close();
    return infected_status;
}


/**
 * Creates a restaurant simulation world from infection data file
 * and adds additional persons up to number_of_persons if specified
 */
mio::abm::World create_world_from_file(const std::string& infection_data_file, int number_of_persons = 0)
{
    // Create world with 6 age groups
    auto world = mio::abm::World(num_age_groups);
    
    // Set infection parameters
    set_parameters(world.parameters);
    
    // Read infection data
    auto infected_status = read_infection_data(infection_data_file);
    
    // Create a hospital and ICU for severe cases
    auto hospital = world.add_location(mio::abm::LocationType::Hospital);
    auto icu = world.add_location(mio::abm::LocationType::ICU);
    
    // Create common locations
    auto event = world.add_location(mio::abm::LocationType::SocialEvent);
    auto shop = world.add_location(mio::abm::LocationType::BasicsShop);
    auto school = world.add_location(mio::abm::LocationType::School);
    auto work = world.add_location(mio::abm::LocationType::Work);
    
    // Create homes (one per table) for the initial persons from the file
    std::map<std::string, mio::abm::LocationId> table_homes;
    for (const auto& [id, infected] : infected_status) {
        // Create home location for each person
        auto home = world.add_location(mio::abm::LocationType::Home);
        table_homes[std::to_string(id)] = home;
    }
    
    // Create people from the infection data and assign them to locations
    std::vector<double> weights_age = {0.05, 0.25, 0.2, 0.4, 0.07, 0.03};
    std::vector<double> weight_home_size = {0.2, 0.3, 0.2, 0.2, 0.1};
    int initial_person_count = 0;
    
    for (const auto& [id, infected] : infected_status) {
        mio::AgeGroup age = mio::AgeGroup(
            mio::DiscreteDistribution<size_t>::get_instance()(
                world.get_rng(), weights_age
            )
        );
        
        // Create person and assign them to their home
        auto& person = world.add_person(table_homes[std::to_string(id)], age);
        person.set_assigned_location(hospital);
        person.set_assigned_location(icu);
        person.set_assigned_location(event);
        person.set_assigned_location(world.get_individualized_location(table_homes[std::to_string(id)]));
        person.set_assigned_location(shop);
        
        // Assign school or work based on age
        if (age == age_group_5_to_14) {
            person.set_assigned_location(school);
        }
        if (age == age_group_15_to_34 || age == age_group_35_to_59) {
            person.set_assigned_location(work);
        }
        
        // If infected according to data, set them as exposed
        if (infected) {
            auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);
            person.add_new_infection(mio::abm::Infection(
                rng, 
                mio::abm::VirusVariant::Wildtype, 
                person.get_age(),
                world.parameters, 
                mio::abm::TimePoint(0), 
                mio::abm::InfectionState::Exposed
            ));
        }
        
        initial_person_count++;
    }
    
    // If a higher number of persons is requested, create additional people
    if (number_of_persons > initial_person_count) {
        // Create additional homes with 1-5 persons per home
        std::vector<mio::abm::LocationId> additional_homes;
        std::vector<int> persons_per_home;
        
        // First create homes and determine how many people in each
        int remaining_persons = number_of_persons - initial_person_count;
        while (remaining_persons > 0) {
            // Randomly decide how many people in this home (1-5)
            int persons_in_home = std::min(remaining_persons, 
                                          1 + static_cast<int>(
                                                mio::DiscreteDistribution<size_t>::get_instance()(
                                                    world.get_rng(), weight_home_size)));
            
            // Create a new home
            auto home = world.add_location(mio::abm::LocationType::Home);
            additional_homes.push_back(home);
            persons_per_home.push_back(persons_in_home);
            
            remaining_persons -= persons_in_home;
        }
        
        // Now create the people and assign them to homes
        for (size_t i = 0; i < additional_homes.size(); ++i) {
            for (int j = 0; j < persons_per_home[i]; ++j) {
                // Randomly determine age based on weights
                mio::AgeGroup age = mio::AgeGroup(
                    mio::DiscreteDistribution<size_t>::get_instance()(
                        world.get_rng(), weights_age
                    )
                );
                
                // Create person and assign them to their home
                auto& person = world.add_person(additional_homes[i], age);
                
                // Assign to common locations
                person.set_assigned_location(world.get_individualized_location(additional_homes[i]));
                person.set_assigned_location(hospital);
                person.set_assigned_location(icu);
                person.set_assigned_location(event);
                person.set_assigned_location(shop);
                
                // Assign school or work based on age
                if (age == age_group_5_to_14) {
                    person.set_assigned_location(school);
                }
                if (age == age_group_15_to_34 || age == age_group_35_to_59) {
                    person.set_assigned_location(work);
                }
            }
        }
    }
    
    set_local_parameters(world);
    std::cout << "Created world with " << world.get_persons().size() << " people" << std::endl;
    std::cout << "Number of homes: " << 
        std::count_if(world.get_locations().begin(), world.get_locations().end(), 
                      [](const mio::abm::Location& loc) { 
                         return loc.get_type() == mio::abm::LocationType::Home; 
                      }) << std::endl;
    
    return world;
}


/**
 * Template function for reading a world from a text file (to be implemented later)
 */
template<typename DataFormat>
mio::abm::World read_world_from_file(const std::string& filename)
{
    mio::unused(filename);
    auto world = mio::abm::World(num_age_groups);
    set_parameters(world.parameters);
    return world;
}


// From: https://en.cppreference.com/w/cpp/chrono/c/strftime
const std::string currentDateTime()
{
    // Example of the very popular RFC 3339 format UTC time
    std::time_t time = std::time({});
    char timeString[std::size("yyyy-mm-dddddddd")];
    std::strftime(std::data(timeString), std::size(timeString), "%F%H%M%S", std::gmtime(&time));
    return timeString;
}

mio::IOResult<bool> create_result_folders(std::string const& result_dir, int n_params = 1)
{
    std::string inf_p_loc_t_p_ag      = result_dir + "/infection_per_location_type_per_age_group/";
    std::string inf_state_p_ag        = result_dir + "/infection_state_per_age_group/";

    BOOST_OUTCOME_TRY(mio::create_directory(result_dir));
    BOOST_OUTCOME_TRY(mio::create_directory(inf_p_loc_t_p_ag));
    BOOST_OUTCOME_TRY(mio::create_directory(inf_state_p_ag));
    if (n_params > 0) {
        // we create n_param folders in each of the subfolders
        for (int i = 0; i < n_params; i++) {
            BOOST_OUTCOME_TRY(mio::create_directory(inf_p_loc_t_p_ag + "/" + std::to_string(i)));
            BOOST_OUTCOME_TRY(mio::create_directory(inf_state_p_ag + "/" + std::to_string(i)));
        }
    }
    return mio::success();
}

mio::IOResult<bool> copy_result_folder(std::string const& from_dir, std::string const& to_dir)
{
    BOOST_OUTCOME_TRY(mio::create_directory(to_dir));
    fs::copy(from_dir, to_dir, fs::copy_options::overwrite_existing | fs::copy_options::recursive);
    return mio::success();
}

mio::IOResult<void> main_flow()
{
    // Default infection data file path
    std::string infection_data_file = "/Users/saschakorf/Nosynch/Arbeit/memilio/memilio/cpp/examples/PanVadere/restaurant/simulation_runs/lu-2020_airflow_inlet_right_outlet_left/infections.txt";
    std::string input_dir = "/Users/saschakorf/Nosynch/Arbeit/memilio/memilio/cpp/examples/results";

    std::string precomputed_dir = input_dir + "/results";
    std::string result_dir      = input_dir + "/results_" + currentDateTime();
    auto created = create_result_folders(result_dir);

    
    std::cout << "Creating restaurant simulation from: " << infection_data_file << std::endl;
    
    int n_persons= 1000;

    // Create the simulation world
    auto world = create_world_from_file(infection_data_file, n_persons);
    
    // Set up simulation timeframe
    auto t0 = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(20);
    
    auto ensemble_infection_per_loc_type =
        std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of infection per location type results
    ensemble_infection_per_loc_type.reserve(size_t(1));

    auto ensemble_infection_state_per_age_group =
        std::vector<std::vector<mio::TimeSeries<ScalarType>>>{}; // Vector of infection state per age group results
    ensemble_infection_state_per_age_group.reserve(size_t(1));

    auto ensemble_params = std::vector<std::vector<mio::abm::World>>{}; // Vector of all worlds
    ensemble_params.reserve(size_t(1));


    // Initialize simulation
    auto sim = mio::abm::Simulation(t0, std::move(world));

    mio::History<mio::abm::TimeSeriesWriter, LogInfectionPerLocationTypePerAgeGroup>
            historyInfectionPerLocationType{
                Eigen::Index((size_t)mio::abm::LocationType::Count * sim.get_world().parameters.get_num_groups())};
    mio::History<mio::abm::TimeSeriesWriter, LogInfectionStatePerAgeGroup> historyInfectionStatePerAgeGroup{
            Eigen::Index((size_t)mio::abm::InfectionState::Count * sim.get_world().parameters.get_num_groups())};
    mio::History<mio::abm::TimeSeriesWriter, mio::abm::LogInfectionState> historyTimeSeries{
        Eigen::Index(mio::abm::InfectionState::Count)
    };
    
    // Run simulation and collect data
    sim.advance(tmax, historyInfectionPerLocationType, historyInfectionStatePerAgeGroup, historyTimeSeries);


    auto temp_sim_infection_per_loc_tpye =
            std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionPerLocationType.get_log())};
    auto temp_sim_infection_state_per_age_group =
        std::vector<mio::TimeSeries<ScalarType>>{std::get<0>(historyInfectionStatePerAgeGroup.get_log())};
        // Push result of the simulation back to the result vector
    ensemble_infection_per_loc_type.emplace_back(temp_sim_infection_per_loc_tpye);
    ensemble_infection_state_per_age_group.emplace_back(temp_sim_infection_state_per_age_group);
    ensemble_params.emplace_back(std::vector<mio::abm::World>{sim.get_world()});
    
    BOOST_OUTCOME_TRY(save_results(ensemble_infection_state_per_age_group, ensemble_params, {0},
            (fs::path)result_dir / "infection_state_per_age_group" / "0", true));
    BOOST_OUTCOME_TRY(save_results(ensemble_infection_per_loc_type, ensemble_params, {0},
            (fs::path)result_dir / "infection_per_location_type_per_age_group" / "0", true));
    

    // copy results into a fixed name folder to have easier access
    std::string last_run_dir = input_dir + "/results_last_run";
    auto copied              = copy_result_folder(result_dir, last_run_dir);

    
    // Write results to file
    std::ofstream outfile("panvXabm_results.txt");
    std::get<0>(historyTimeSeries.get_log())
        .print_table({"S", "E", "I_NS", "I_Sy", "I_Sev", "I_Crit", "R", "D"}, 7, 4, outfile);
    
    std::cout << "Results written to panvXabm_results.txt" << std::endl;
    

    
    return mio::success();
}

int main()
{
    auto result = main_flow();
    if (result.has_error()) {
        std::cerr << "Error: " << result.error().message() << std::endl;
        return 1;
    }
    return 0;
}