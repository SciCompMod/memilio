/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding, Sascha Korf
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
#include "abm/result_simulation.h"
#include "abm/household.h"
#include "abm/lockdown_rules.h"
#include "abm/model.h"
#include "abm/time.h"
#include "abm/protection_event.h"

#include "boost/filesystem.hpp"
#include "boost/algorithm/string/split.hpp"
#include "boost/algorithm/string/classification.hpp"

#include "memilio/compartments/parameter_studies.h"
#include "memilio/data/analyze_result.h"
#include "memilio/io/io.h"
#include "memilio/io/directories.h"
#include "memilio/io/result_io.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/stl_util.h"

#include <string>

namespace fs                    = std::filesystem;
constexpr size_t num_age_groups = 6;

const auto age_group_0_to_4   = mio::AgeGroup(0);
const auto age_group_5_to_14  = mio::AgeGroup(1);
const auto age_group_15_to_34 = mio::AgeGroup(2);
const auto age_group_35_to_59 = mio::AgeGroup(3);
const auto age_group_60_to_79 = mio::AgeGroup(4);
const auto age_group_80_plus  = mio::AgeGroup(5);

std::pair<double, double> get_my_and_sigma(std::pair<double, double> mean_and_std)
{
    auto mean    = mean_and_std.first;
    auto stddev  = mean_and_std.second;
    double my    = log(mean * mean / sqrt(mean * mean + stddev * stddev));
    double sigma = sqrt(log(1 + stddev * stddev / (mean * mean)));
    return {my, sigma};
}

void set_parameters(mio::abm::Parameters& params)
{
    // Set the Time parameters for the infection same for every age group for now

    auto my_and_sigma_exposed = get_my_and_sigma({4.5, 1.5});
    params.get<mio::abm::TimeExposedToNoSymptoms>() =
        mio::ParameterDistributionLogNormal(my_and_sigma_exposed.first, my_and_sigma_exposed.second);

    auto my_and_sigma_no_symptoms_to_symptoms                = get_my_and_sigma({1.1, 0.9});
    params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>() = mio::ParameterDistributionLogNormal(
        my_and_sigma_no_symptoms_to_symptoms.first, my_and_sigma_no_symptoms_to_symptoms.second);

    auto my_and_sigma_no_symptoms_to_recovered                = get_my_and_sigma({8.0, 2.0});
    params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>() = mio::ParameterDistributionLogNormal(
        my_and_sigma_no_symptoms_to_recovered.first, my_and_sigma_no_symptoms_to_recovered.second);

    auto my_and_sigma_symptoms_to_severe                 = get_my_and_sigma({6.6, 4.9});
    params.get<mio::abm::TimeInfectedSymptomsToSevere>() = mio::ParameterDistributionLogNormal(
        my_and_sigma_symptoms_to_severe.first, my_and_sigma_symptoms_to_severe.second);

    auto my_and_sigma_symptoms_to_recovered                 = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedSymptomsToRecovered>() = mio::ParameterDistributionLogNormal(
        my_and_sigma_symptoms_to_recovered.first, my_and_sigma_symptoms_to_recovered.second);

    auto my_and_sigma_severe_to_critical                 = get_my_and_sigma({1.5, 2.0});
    params.get<mio::abm::TimeInfectedSevereToCritical>() = mio::ParameterDistributionLogNormal(
        my_and_sigma_severe_to_critical.first, my_and_sigma_severe_to_critical.second);

    auto my_and_sigma_severe_to_recovered                 = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedSevereToRecovered>() = mio::ParameterDistributionLogNormal(
        my_and_sigma_severe_to_recovered.first, my_and_sigma_severe_to_recovered.second);

    auto my_and_sigma_severe_to_dead = get_my_and_sigma({10.7, 4.8});
    params.get<mio::abm::TimeInfectedSevereToDead>() =
        mio::ParameterDistributionLogNormal(my_and_sigma_severe_to_dead.first, my_and_sigma_severe_to_dead.second);

    auto my_and_sigma_critical_to_dead = get_my_and_sigma({10.7, 4.8});
    params.get<mio::abm::TimeInfectedCriticalToDead>() =
        mio::ParameterDistributionLogNormal(my_and_sigma_critical_to_dead.first, my_and_sigma_critical_to_dead.second);

    auto my_and_sigma_critical_to_recovered                 = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedCriticalToRecovered>() = mio::ParameterDistributionLogNormal(
        my_and_sigma_critical_to_recovered.first, my_and_sigma_critical_to_recovered.second);

    //Set testing parameters
    auto pcr_test_values = mio::abm::TestParameters{0.9, 0.995, mio::abm::hours(24), mio::abm::TestType::PCR};
    auto antigen_test_values =
        mio::abm::TestParameters{0.9, 0.9999, mio::abm::minutes(15), mio::abm::TestType::Antigen};
    auto generic_test_values = mio::abm::TestParameters{0.7, 0.95, mio::abm::minutes(30), mio::abm::TestType::Generic};

    params.get<mio::abm::TestData>()[mio::abm::TestType::PCR]     = pcr_test_values;
    params.get<mio::abm::TestData>()[mio::abm::TestType::Antigen] = antigen_test_values;
    params.get<mio::abm::TestData>()[mio::abm::TestType::Generic] = generic_test_values;

    // Set percentage parameters
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]  = 0.50;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}] = 0.55;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] =
        0.60;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] =
        0.70;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] =
        0.83;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}] = 0.90;

    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.02;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]  = 0.03;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.04;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.07;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.17;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]  = 0.24;

    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.1;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]  = 0.11;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.12;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.14;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.33;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]  = 0.62;

    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_0_to_4}]   = 0.12;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_5_to_14}]  = 0.13;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34}] = 0.15;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59}] = 0.26;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79}] = 0.40;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus}]  = 0.48;

    // Set infection parameters
    // Set protection level against an severe infection.
    params.get<mio::abm::SeverityProtectionFactor>() =
        mio::TimeSeriesFunctor<ScalarType>{mio::TimeSeriesFunctorType::LinearInterpolation, {{0, 0.8}, {150, 0.8}}};
    params.get<mio::abm::HighViralLoadProtectionFactor>() =
        mio::TimeSeriesFunctor<ScalarType>{mio::TimeSeriesFunctorType::LinearInterpolation, {{0, 0.0}, {150, 0.0}}};

    // Set PAIS parameters
    // PAISProbability is 0.0 by default. Set only values that are different from 0.0 here.
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34,
                                             mio::abm::Sex::Female, mio::abm::VaccinationClass::Zero}] = 0.108;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34, mio::abm::Sex::Male,
                                             mio::abm::VaccinationClass::Zero}]                        = 0.051;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59,
                                             mio::abm::Sex::Female, mio::abm::VaccinationClass::Zero}] = 0.136;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59, mio::abm::Sex::Male,
                                             mio::abm::VaccinationClass::Zero}]                        = 0.105;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79,
                                             mio::abm::Sex::Female, mio::abm::VaccinationClass::Zero}] = 0.188;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79, mio::abm::Sex::Male,
                                             mio::abm::VaccinationClass::Zero}]                        = 0.108;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus, mio::abm::Sex::Female,
                                             mio::abm::VaccinationClass::Zero}]                        = 0.192;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus, mio::abm::Sex::Male,
                                             mio::abm::VaccinationClass::Zero}]                        = 0.175;

    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34,
                                             mio::abm::Sex::Female, mio::abm::VaccinationClass::OneOrTwo}] = 0.167;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34, mio::abm::Sex::Male,
                                             mio::abm::VaccinationClass::OneOrTwo}]                        = 0.098;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59,
                                             mio::abm::Sex::Female, mio::abm::VaccinationClass::OneOrTwo}] = 0.189;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59, mio::abm::Sex::Male,
                                             mio::abm::VaccinationClass::OneOrTwo}]                        = 0.133;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79,
                                             mio::abm::Sex::Female, mio::abm::VaccinationClass::OneOrTwo}] = 0.166;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79, mio::abm::Sex::Male,
                                             mio::abm::VaccinationClass::OneOrTwo}]                        = 0.109;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus, mio::abm::Sex::Female,
                                             mio::abm::VaccinationClass::OneOrTwo}]                        = 0.137;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus, mio::abm::Sex::Male,
                                             mio::abm::VaccinationClass::OneOrTwo}]                        = 0.078;

    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34,
                                             mio::abm::Sex::Female, mio::abm::VaccinationClass::ThreeOrMore}] = 0.156;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_15_to_34, mio::abm::Sex::Male,
                                             mio::abm::VaccinationClass::ThreeOrMore}]                        = 0.107;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59,
                                             mio::abm::Sex::Female, mio::abm::VaccinationClass::ThreeOrMore}] = 0.193;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_35_to_59, mio::abm::Sex::Male,
                                             mio::abm::VaccinationClass::ThreeOrMore}]                        = 0.122;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79,
                                             mio::abm::Sex::Female, mio::abm::VaccinationClass::ThreeOrMore}] = 0.163;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_60_to_79, mio::abm::Sex::Male,
                                             mio::abm::VaccinationClass::ThreeOrMore}]                        = 0.106;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus, mio::abm::Sex::Female,
                                             mio::abm::VaccinationClass::ThreeOrMore}]                        = 0.154;
    params.get<mio::abm::PAISProbability>()[{mio::abm::VirusVariant::Wildtype, age_group_80_plus, mio::abm::Sex::Male,
                                             mio::abm::VaccinationClass::ThreeOrMore}]                        = 0.101;

    params.get<mio::abm::PAISProbabilitySeverityFactor>()[{mio::abm::VirusVariant::Wildtype,
                                                           mio::abm::VaccinationClass::Zero}]          = 0.267;
    params.get<mio::abm::PAISProbabilitySeverityFactor>()[{mio::abm::VirusVariant::Wildtype,
                                                           mio::abm::VaccinationClass::OneOrTwo}]      = 0.247;
    params.get<mio::abm::PAISProbabilitySeverityFactor>()[{mio::abm::VirusVariant::Wildtype,
                                                           mio::abm::VaccinationClass::ThreeOrMore}]   = 0.217;
    params.get<mio::abm::PAISProtectionAtSecondInfection>()[{mio::abm::VirusVariant::Wildtype,
                                                             mio::abm::VaccinationClass::Zero}]        = 0.038;
    params.get<mio::abm::PAISProtectionAtSecondInfection>()[{mio::abm::VirusVariant::Wildtype,
                                                             mio::abm::VaccinationClass::OneOrTwo}]    = 0.116;
    params.get<mio::abm::PAISProtectionAtSecondInfection>()[{mio::abm::VirusVariant::Wildtype,
                                                             mio::abm::VaccinationClass::ThreeOrMore}] = 0.118;

    Eigen::MatrixXd pais_transition_matrix = Eigen::MatrixXd::Zero(
        static_cast<Eigen::Index>(mio::abm::PAISState::Count), static_cast<Eigen::Index>(mio::abm::PAISState::Count));
    pais_transition_matrix(static_cast<Eigen::Index>(mio::abm::PAISState::Healthy),
                           static_cast<Eigen::Index>(mio::abm::PAISState::Healthy)) = 0.9925839491519293;
    pais_transition_matrix(static_cast<Eigen::Index>(mio::abm::PAISState::Healthy),
                           static_cast<Eigen::Index>(mio::abm::PAISState::Medium))  = 0.00026973522588106033;
    pais_transition_matrix(static_cast<Eigen::Index>(mio::abm::PAISState::Healthy),
                           static_cast<Eigen::Index>(mio::abm::PAISState::Severe))  = 0.007146315622189754;
    pais_transition_matrix(static_cast<Eigen::Index>(mio::abm::PAISState::Medium),
                           static_cast<Eigen::Index>(mio::abm::PAISState::Healthy)) = 1.879727103567679e-5;
    pais_transition_matrix(static_cast<Eigen::Index>(mio::abm::PAISState::Medium),
                           static_cast<Eigen::Index>(mio::abm::PAISState::Medium))  = 1.1042796363187963e-9;
    pais_transition_matrix(static_cast<Eigen::Index>(mio::abm::PAISState::Medium),
                           static_cast<Eigen::Index>(mio::abm::PAISState::Severe))  = 0.9999812016246846;
    pais_transition_matrix(static_cast<Eigen::Index>(mio::abm::PAISState::Severe),
                           static_cast<Eigen::Index>(mio::abm::PAISState::Healthy)) = 0.23940186497365926;
    pais_transition_matrix(static_cast<Eigen::Index>(mio::abm::PAISState::Severe),
                           static_cast<Eigen::Index>(mio::abm::PAISState::Medium))  = 0.6049572234530622;
    pais_transition_matrix(static_cast<Eigen::Index>(mio::abm::PAISState::Severe),
                           static_cast<Eigen::Index>(mio::abm::PAISState::Severe))  = 0.15564091157327845;
    params.get<mio::abm::PAISTransitionMatrix>()                                    = pais_transition_matrix;

    //Set other parameters
    params.get<mio::abm::AerosolTransmissionRates>()                                       = 0.0;
    params.get<mio::abm::InfectionRateFromViralShed>()[{mio::abm::VirusVariant::Wildtype}] = 1.6;
    //params.get<mio::abm::MaskProtection>()             = 0.5;
    //params.get<mio::abm::MaskProtection>({mio::abm::MaskType::None}) = 0.0;
    params.get<mio::abm::QuarantineDuration>()      = mio::abm::days(10);
    params.get<mio::abm::QuarantineEffectiveness>() = 0.5;
}

// set location specific parameters
void set_local_parameters(mio::abm::Model& model)
{
    const int n_age_groups = (int)model.parameters.get_num_groups();

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
    mio::ContactMatrix<ScalarType> contacts_home(static_cast<Eigen::Index>(n_age_groups));
    contacts_home.get_baseline()(age_group_0_to_4.get(), age_group_0_to_4.get())     = 0.4413;
    contacts_home.get_baseline()(age_group_0_to_4.get(), age_group_5_to_14.get())    = 0.0504;
    contacts_home.get_baseline()(age_group_0_to_4.get(), age_group_15_to_34.get())   = 1.2383;
    contacts_home.get_baseline()(age_group_0_to_4.get(), age_group_35_to_59.get())   = 0.8033;
    contacts_home.get_baseline()(age_group_0_to_4.get(), age_group_60_to_79.get())   = 0.0494;
    contacts_home.get_baseline()(age_group_0_to_4.get(), age_group_80_plus.get())    = 0.0017;
    contacts_home.get_baseline()(age_group_5_to_14.get(), age_group_0_to_4.get())    = 0.0485;
    contacts_home.get_baseline()(age_group_5_to_14.get(), age_group_5_to_14.get())   = 0.7616;
    contacts_home.get_baseline()(age_group_5_to_14.get(), age_group_15_to_34.get())  = 0.6532;
    contacts_home.get_baseline()(age_group_5_to_14.get(), age_group_35_to_59.get())  = 1.1614;
    contacts_home.get_baseline()(age_group_5_to_14.get(), age_group_60_to_79.get())  = 0.0256;
    contacts_home.get_baseline()(age_group_5_to_14.get(), age_group_80_plus.get())   = 0.0013;
    contacts_home.get_baseline()(age_group_15_to_34.get(), age_group_0_to_4.get())   = 0.1800;
    contacts_home.get_baseline()(age_group_15_to_34.get(), age_group_5_to_14.get())  = 0.1795;
    contacts_home.get_baseline()(age_group_15_to_34.get(), age_group_15_to_34.get()) = 0.8806;
    contacts_home.get_baseline()(age_group_15_to_34.get(), age_group_35_to_59.get()) = 0.6413;
    contacts_home.get_baseline()(age_group_15_to_34.get(), age_group_60_to_79.get()) = 0.0429;
    contacts_home.get_baseline()(age_group_15_to_34.get(), age_group_80_plus.get())  = 0.0032;
    contacts_home.get_baseline()(age_group_35_to_59.get(), age_group_0_to_4.get())   = 0.0495;
    contacts_home.get_baseline()(age_group_35_to_59.get(), age_group_5_to_14.get())  = 0.2639;
    contacts_home.get_baseline()(age_group_35_to_59.get(), age_group_15_to_34.get()) = 0.5189;
    contacts_home.get_baseline()(age_group_35_to_59.get(), age_group_35_to_59.get()) = 0.8277;
    contacts_home.get_baseline()(age_group_35_to_59.get(), age_group_60_to_79.get()) = 0.0679;
    contacts_home.get_baseline()(age_group_35_to_59.get(), age_group_80_plus.get())  = 0.0014;
    contacts_home.get_baseline()(age_group_60_to_79.get(), age_group_0_to_4.get())   = 0.0087;
    contacts_home.get_baseline()(age_group_60_to_79.get(), age_group_5_to_14.get())  = 0.0394;
    contacts_home.get_baseline()(age_group_60_to_79.get(), age_group_15_to_34.get()) = 0.1417;
    contacts_home.get_baseline()(age_group_60_to_79.get(), age_group_35_to_59.get()) = 0.3834;
    contacts_home.get_baseline()(age_group_60_to_79.get(), age_group_60_to_79.get()) = 0.7064;
    contacts_home.get_baseline()(age_group_60_to_79.get(), age_group_80_plus.get())  = 0.0447;
    contacts_home.get_baseline()(age_group_80_plus.get(), age_group_0_to_4.get())    = 0.0292;
    contacts_home.get_baseline()(age_group_80_plus.get(), age_group_5_to_14.get())   = 0.0648;
    contacts_home.get_baseline()(age_group_80_plus.get(), age_group_15_to_34.get())  = 0.1248;
    contacts_home.get_baseline()(age_group_80_plus.get(), age_group_35_to_59.get())  = 0.4179;
    contacts_home.get_baseline()(age_group_80_plus.get(), age_group_60_to_79.get())  = 0.3497;
    contacts_home.get_baseline()(age_group_80_plus.get(), age_group_80_plus.get())   = 0.1544;

    /* baseline_school
        1.1165 0.2741 0.2235 0.1028 0.0007 0.0000
        0.1627 1.9412 0.2431 0.1780 0.0130 0.0000
        0.0148 0.1646 1.1266 0.0923 0.0074 0.0000
        0.0367 0.1843 0.3265 0.0502 0.0021 0.0005
        0.0004 0.0370 0.0115 0.0014 0.0039 0.0000
        0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    */
    mio::ContactMatrix<ScalarType> contacts_school(static_cast<Eigen::Index>(n_age_groups));
    contacts_school.get_baseline()(age_group_0_to_4.get(), age_group_0_to_4.get())     = 1.1165;
    contacts_school.get_baseline()(age_group_0_to_4.get(), age_group_5_to_14.get())    = 0.2741;
    contacts_school.get_baseline()(age_group_0_to_4.get(), age_group_15_to_34.get())   = 0.2235;
    contacts_school.get_baseline()(age_group_0_to_4.get(), age_group_35_to_59.get())   = 0.1028;
    contacts_school.get_baseline()(age_group_0_to_4.get(), age_group_60_to_79.get())   = 0.0007;
    contacts_school.get_baseline()(age_group_0_to_4.get(), age_group_80_plus.get())    = 0.0000;
    contacts_school.get_baseline()(age_group_5_to_14.get(), age_group_0_to_4.get())    = 0.1627;
    contacts_school.get_baseline()(age_group_5_to_14.get(), age_group_5_to_14.get())   = 1.9412;
    contacts_school.get_baseline()(age_group_5_to_14.get(), age_group_15_to_34.get())  = 0.2431;
    contacts_school.get_baseline()(age_group_5_to_14.get(), age_group_35_to_59.get())  = 0.1780;
    contacts_school.get_baseline()(age_group_5_to_14.get(), age_group_60_to_79.get())  = 0.0130;
    contacts_school.get_baseline()(age_group_5_to_14.get(), age_group_80_plus.get())   = 0.0000;
    contacts_school.get_baseline()(age_group_15_to_34.get(), age_group_0_to_4.get())   = 0.0148;
    contacts_school.get_baseline()(age_group_15_to_34.get(), age_group_5_to_14.get())  = 0.1646;
    contacts_school.get_baseline()(age_group_15_to_34.get(), age_group_15_to_34.get()) = 1.1266;
    contacts_school.get_baseline()(age_group_15_to_34.get(), age_group_35_to_59.get()) = 0.0923;
    contacts_school.get_baseline()(age_group_15_to_34.get(), age_group_60_to_79.get()) = 0.0074;
    contacts_school.get_baseline()(age_group_15_to_34.get(), age_group_80_plus.get())  = 0.0000;
    contacts_school.get_baseline()(age_group_35_to_59.get(), age_group_0_to_4.get())   = 0.0367;
    contacts_school.get_baseline()(age_group_35_to_59.get(), age_group_5_to_14.get())  = 0.1843;
    contacts_school.get_baseline()(age_group_35_to_59.get(), age_group_15_to_34.get()) = 0.3265;
    contacts_school.get_baseline()(age_group_35_to_59.get(), age_group_35_to_59.get()) = 0.0502;
    contacts_school.get_baseline()(age_group_35_to_59.get(), age_group_60_to_79.get()) = 0.0021;
    contacts_school.get_baseline()(age_group_35_to_59.get(), age_group_80_plus.get())  = 0.0005;
    contacts_school.get_baseline()(age_group_60_to_79.get(), age_group_0_to_4.get())   = 0.0004;
    contacts_school.get_baseline()(age_group_60_to_79.get(), age_group_5_to_14.get())  = 0.0370;
    contacts_school.get_baseline()(age_group_60_to_79.get(), age_group_15_to_34.get()) = 0.0115;
    contacts_school.get_baseline()(age_group_60_to_79.get(), age_group_35_to_59.get()) = 0.0014;
    contacts_school.get_baseline()(age_group_60_to_79.get(), age_group_60_to_79.get()) = 0.0039;
    contacts_school.get_baseline()(age_group_60_to_79.get(), age_group_80_plus.get())  = 0.0000;
    contacts_school.get_baseline()(age_group_80_plus.get(), age_group_0_to_4.get())    = 0.0000;
    contacts_school.get_baseline()(age_group_80_plus.get(), age_group_5_to_14.get())   = 0.0000;
    contacts_school.get_baseline()(age_group_80_plus.get(), age_group_15_to_34.get())  = 0.0000;
    contacts_school.get_baseline()(age_group_80_plus.get(), age_group_35_to_59.get())  = 0.0000;
    contacts_school.get_baseline()(age_group_80_plus.get(), age_group_60_to_79.get())  = 0.0000;
    contacts_school.get_baseline()(age_group_80_plus.get(), age_group_80_plus.get())   = 0.0000;

    /* baseline_work
        0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
        0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
        0.0000 0.0127 1.7570 1.6050 0.0133 0.0000
        0.0000 0.0020 1.0311 2.3166 0.0098 0.0000
        0.0000 0.0002 0.0194 0.0325 0.0003 0.0000
        0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
    */
    mio::ContactMatrix<ScalarType> contacts_work(static_cast<Eigen::Index>(n_age_groups));
    contacts_work.get_baseline()(age_group_0_to_4.get(), age_group_0_to_4.get())     = 0.0000;
    contacts_work.get_baseline()(age_group_0_to_4.get(), age_group_5_to_14.get())    = 0.0000;
    contacts_work.get_baseline()(age_group_0_to_4.get(), age_group_15_to_34.get())   = 0.0000;
    contacts_work.get_baseline()(age_group_0_to_4.get(), age_group_35_to_59.get())   = 0.0000;
    contacts_work.get_baseline()(age_group_0_to_4.get(), age_group_60_to_79.get())   = 0.0000;
    contacts_work.get_baseline()(age_group_0_to_4.get(), age_group_80_plus.get())    = 0.0000;
    contacts_work.get_baseline()(age_group_5_to_14.get(), age_group_0_to_4.get())    = 0.0000;
    contacts_work.get_baseline()(age_group_5_to_14.get(), age_group_5_to_14.get())   = 0.0000;
    contacts_work.get_baseline()(age_group_5_to_14.get(), age_group_15_to_34.get())  = 0.0000;
    contacts_work.get_baseline()(age_group_5_to_14.get(), age_group_35_to_59.get())  = 0.0000;
    contacts_work.get_baseline()(age_group_5_to_14.get(), age_group_60_to_79.get())  = 0.0000;
    contacts_work.get_baseline()(age_group_5_to_14.get(), age_group_80_plus.get())   = 0.0000;
    contacts_work.get_baseline()(age_group_15_to_34.get(), age_group_0_to_4.get())   = 0.0000;
    contacts_work.get_baseline()(age_group_15_to_34.get(), age_group_5_to_14.get())  = 0.0127;
    contacts_work.get_baseline()(age_group_15_to_34.get(), age_group_15_to_34.get()) = 1.7570;
    contacts_work.get_baseline()(age_group_15_to_34.get(), age_group_35_to_59.get()) = 1.6050;
    contacts_work.get_baseline()(age_group_15_to_34.get(), age_group_60_to_79.get()) = 0.0133;
    contacts_work.get_baseline()(age_group_15_to_34.get(), age_group_80_plus.get())  = 0.0000;
    contacts_work.get_baseline()(age_group_35_to_59.get(), age_group_0_to_4.get())   = 0.0000;
    contacts_work.get_baseline()(age_group_35_to_59.get(), age_group_5_to_14.get())  = 0.0020;
    contacts_work.get_baseline()(age_group_35_to_59.get(), age_group_15_to_34.get()) = 1.0311;
    contacts_work.get_baseline()(age_group_35_to_59.get(), age_group_35_to_59.get()) = 2.3166;
    contacts_work.get_baseline()(age_group_35_to_59.get(), age_group_60_to_79.get()) = 0.0098;
    contacts_work.get_baseline()(age_group_35_to_59.get(), age_group_80_plus.get())  = 0.0000;
    contacts_work.get_baseline()(age_group_60_to_79.get(), age_group_0_to_4.get())   = 0.0000;
    contacts_work.get_baseline()(age_group_60_to_79.get(), age_group_5_to_14.get())  = 0.0002;
    contacts_work.get_baseline()(age_group_60_to_79.get(), age_group_15_to_34.get()) = 0.0194;
    contacts_work.get_baseline()(age_group_60_to_79.get(), age_group_35_to_59.get()) = 0.0325;
    contacts_work.get_baseline()(age_group_60_to_79.get(), age_group_60_to_79.get()) = 0.0003;
    contacts_work.get_baseline()(age_group_60_to_79.get(), age_group_80_plus.get())  = 0.0000;
    contacts_work.get_baseline()(age_group_80_plus.get(), age_group_0_to_4.get())    = 0.0000;
    contacts_work.get_baseline()(age_group_80_plus.get(), age_group_5_to_14.get())   = 0.0000;
    contacts_work.get_baseline()(age_group_80_plus.get(), age_group_15_to_34.get())  = 0.0000;
    contacts_work.get_baseline()(age_group_80_plus.get(), age_group_35_to_59.get())  = 0.0000;
    contacts_work.get_baseline()(age_group_80_plus.get(), age_group_60_to_79.get())  = 0.0000;
    contacts_work.get_baseline()(age_group_80_plus.get(), age_group_80_plus.get())   = 0.0000;

    /* baseline_other
        0.5170 0.3997 0.7957 0.9958 0.3239 0.0428
        0.0632 0.9121 0.3254 0.4731 0.2355 0.0148
        0.0336 0.1604 1.7529 0.8622 0.1440 0.0077
        0.0204 0.1444 0.5738 1.2127 0.3433 0.0178
        0.0371 0.0393 0.4171 0.9666 0.7495 0.0257
        0.0791 0.0800 0.3480 0.5588 0.2769 0.0180
    */
    mio::ContactMatrix<ScalarType> contacts_other(static_cast<Eigen::Index>(n_age_groups));
    contacts_other.get_baseline()(age_group_0_to_4.get(), age_group_0_to_4.get())     = 0.5170;
    contacts_other.get_baseline()(age_group_0_to_4.get(), age_group_5_to_14.get())    = 0.3997;
    contacts_other.get_baseline()(age_group_0_to_4.get(), age_group_15_to_34.get())   = 0.7957;
    contacts_other.get_baseline()(age_group_0_to_4.get(), age_group_35_to_59.get())   = 0.9958;
    contacts_other.get_baseline()(age_group_0_to_4.get(), age_group_60_to_79.get())   = 0.3239;
    contacts_other.get_baseline()(age_group_0_to_4.get(), age_group_80_plus.get())    = 0.0428;
    contacts_other.get_baseline()(age_group_5_to_14.get(), age_group_0_to_4.get())    = 0.0632;
    contacts_other.get_baseline()(age_group_5_to_14.get(), age_group_5_to_14.get())   = 0.9121;
    contacts_other.get_baseline()(age_group_5_to_14.get(), age_group_15_to_34.get())  = 0.3254;
    contacts_other.get_baseline()(age_group_5_to_14.get(), age_group_35_to_59.get())  = 0.4731;
    contacts_other.get_baseline()(age_group_5_to_14.get(), age_group_60_to_79.get())  = 0.2355;
    contacts_other.get_baseline()(age_group_5_to_14.get(), age_group_80_plus.get())   = 0.0148;
    contacts_other.get_baseline()(age_group_15_to_34.get(), age_group_0_to_4.get())   = 0.0336;
    contacts_other.get_baseline()(age_group_15_to_34.get(), age_group_5_to_14.get())  = 0.1604;
    contacts_other.get_baseline()(age_group_15_to_34.get(), age_group_15_to_34.get()) = 1.7529;
    contacts_other.get_baseline()(age_group_15_to_34.get(), age_group_35_to_59.get()) = 0.8622;
    contacts_other.get_baseline()(age_group_15_to_34.get(), age_group_60_to_79.get()) = 0.1440;
    contacts_other.get_baseline()(age_group_15_to_34.get(), age_group_80_plus.get())  = 0.0077;
    contacts_other.get_baseline()(age_group_35_to_59.get(), age_group_0_to_4.get())   = 0.0204;
    contacts_other.get_baseline()(age_group_35_to_59.get(), age_group_5_to_14.get())  = 0.1444;
    contacts_other.get_baseline()(age_group_35_to_59.get(), age_group_15_to_34.get()) = 0.5738;
    contacts_other.get_baseline()(age_group_35_to_59.get(), age_group_35_to_59.get()) = 1.2127;
    contacts_other.get_baseline()(age_group_35_to_59.get(), age_group_60_to_79.get()) = 0.3433;
    contacts_other.get_baseline()(age_group_35_to_59.get(), age_group_80_plus.get())  = 0.0178;
    contacts_other.get_baseline()(age_group_60_to_79.get(), age_group_0_to_4.get())   = 0.0371;
    contacts_other.get_baseline()(age_group_60_to_79.get(), age_group_5_to_14.get())  = 0.0393;
    contacts_other.get_baseline()(age_group_60_to_79.get(), age_group_15_to_34.get()) = 0.4171;
    contacts_other.get_baseline()(age_group_60_to_79.get(), age_group_35_to_59.get()) = 0.9666;
    contacts_other.get_baseline()(age_group_60_to_79.get(), age_group_60_to_79.get()) = 0.7495;
    contacts_other.get_baseline()(age_group_60_to_79.get(), age_group_80_plus.get())  = 0.0257;
    contacts_other.get_baseline()(age_group_80_plus.get(), age_group_0_to_4.get())    = 0.0791;
    contacts_other.get_baseline()(age_group_80_plus.get(), age_group_5_to_14.get())   = 0.0800;
    contacts_other.get_baseline()(age_group_80_plus.get(), age_group_15_to_34.get())  = 0.3480;
    contacts_other.get_baseline()(age_group_80_plus.get(), age_group_35_to_59.get())  = 0.5588;
    contacts_other.get_baseline()(age_group_80_plus.get(), age_group_60_to_79.get())  = 0.2769;
    contacts_other.get_baseline()(age_group_80_plus.get(), age_group_80_plus.get())   = 0.0180;

    mio::ContactMatrix<ScalarType> contacts_random(static_cast<Eigen::Index>(n_age_groups));

    for (auto& loc : model.get_locations()) {
        switch (loc.get_type()) {
        case mio::abm::LocationType::Home:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_home;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 1.6; //15 hours
            break;
        case mio::abm::LocationType::School:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_school;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 12.0; //2 hours
            break;
        case mio::abm::LocationType::Work:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_work;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 8.0; // 3 hours
            break;
        case mio::abm::LocationType::SocialEvent:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 1.2;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 8.0; // 3 hours
            break;
        case mio::abm::LocationType::BasicsShop:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 0.8;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 12.0; // 2 hours
            break;
        default:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_random;
            break;
        }
    }
}

// map age groups from file to model age groups
mio::AgeGroup get_age_group_from_string(const std::string& age_group_string)
{
    if (age_group_string == "18-20") {
        return age_group_15_to_34;
    }
    if (age_group_string == "20-29") {
        return age_group_15_to_34;
    }
    else if (age_group_string == "30-39") {
        // randomly 15-34 or 35-59
        auto rng = mio::RandomNumberGenerator();
        if (mio::UniformDistribution<ScalarType>::get_instance()(rng, 0.0, 1.0) < 0.5) {
            return age_group_15_to_34;
        }
        return age_group_35_to_59;
    }
    else if (age_group_string == "40-49") {
        return age_group_35_to_59;
    }
    else if (age_group_string == "50-59") {
        return age_group_35_to_59;
    }
    else if (age_group_string == "60-69") {
        return age_group_60_to_79;
    }
    else if (age_group_string == "70+") {
        // randomly 60-79 or 80+
        auto rng = mio::RandomNumberGenerator();
        if (mio::UniformDistribution<ScalarType>::get_instance()(rng, 0.0, 1.0) < 0.5) {
            return age_group_60_to_79;
        }
        return age_group_80_plus;
    }
    else {
        mio::log_error("Invalid age group string.");
        return age_group_0_to_4; // default value
    }
}

// convert date string to TimePoint
mio::abm::TimePoint get_time_point_from_string(const std::string& date_string)
{
    // split the string into year, month and day and convert to time_t
    std::vector<std::string> values;
    boost::split(values, date_string, boost::is_any_of("-"));
    struct tm datetime{};
    datetime.tm_year  = std::stoi(values[0]) - 1900;
    datetime.tm_mon   = std::stoi(values[1]) - 1;
    datetime.tm_mday  = std::stoi(values[2]);
    datetime.tm_isdst = -1;
    time_t timestamp  = mktime(&datetime);
    return mio::abm::TimePoint((int)timestamp);
}

// add infection and vaccination dates to vectors
void add_infection_dates_and_vaccinations(std::vector<std::vector<mio::abm::TimePoint>>& infections,
                                          std::vector<std::vector<mio::abm::TimePoint>>& vaccinations,
                                          const std::vector<std::string>& values,
                                          const std::map<std::string, uint32_t>& index, mio::abm::TimePoint start_date)
{
    std::vector<mio::abm::TimePoint> infection_dates;
    if (values[index.at("ih_infection_1_date")] != "NA") {
        mio::abm::TimePoint infection_date = get_time_point_from_string(values[index.at("ih_infection_1_date")]);
        if (infection_date < start_date) {
            infection_dates.push_back(infection_date);
        }
    }
    if (values[index.at("ih_infection_2_date")] != "NA") {
        mio::abm::TimePoint infection_date = get_time_point_from_string(values[index.at("ih_infection_2_date")]);
        if (infection_date < start_date) {
            infection_dates.push_back(infection_date);
        }
    }
    if (values[index.at("ih_infection_3_date")] != "NA") {
        mio::abm::TimePoint infection_date = get_time_point_from_string(values[index.at("ih_infection_3_date")]);
        if (infection_date < start_date) {
            infection_dates.push_back(infection_date);
        }
    }
    if (values[index.at("ih_infection_4_date")] != "NA") {
        mio::abm::TimePoint infection_date = get_time_point_from_string(values[index.at("ih_infection_4_date")]);
        if (infection_date < start_date) {
            infection_dates.push_back(infection_date);
        }
    }
    infections.push_back(infection_dates);

    std::vector<mio::abm::TimePoint> vaccination_dates;
    if (values[index.at("ih_vaccine_1_date")] != "NA") {
        mio::abm::TimePoint vaccination_date = get_time_point_from_string(values[index.at("ih_vaccine_1_date")]);
        if (vaccination_date < start_date) {
            vaccination_dates.push_back(vaccination_date);
        }
    }
    if (values[index.at("ih_vaccine_2_date")] != "NA") {
        mio::abm::TimePoint vaccination_date = get_time_point_from_string(values[index.at("ih_vaccine_2_date")]);
        if (vaccination_date < start_date) {
            vaccination_dates.push_back(vaccination_date);
        }
    }
    if (values[index.at("ih_vaccine_3_date")] != "NA") {
        mio::abm::TimePoint vaccination_date = get_time_point_from_string(values[index.at("ih_vaccine_3_date")]);
        if (vaccination_date < start_date) {
            vaccination_dates.push_back(vaccination_date);
        }
    }
    vaccinations.push_back(vaccination_dates);
}

void add_infection(mio::abm::Model& model, mio::abm::PersonId pid, const mio::abm::TimePoint infection_date)
{
    mio::abm::PersonalRandomNumberGenerator prng(model.get_rng(), model.get_person(pid));
    model.get_person(pid).add_new_infection(mio::abm::Infection(prng, mio::abm::VirusVariant::Wildtype,
                                                                model.get_person(pid).get_age(), model.parameters,
                                                                infection_date),
                                            prng, infection_date, model.parameters);
}

// add vaccination by date
void add_vaccination(mio::abm::Model& model, mio::abm::PersonId pid, const mio::abm::TimePoint vaccination_date)
{
    model.get_person(pid).add_new_vaccination(mio::abm::ProtectionType::GenericVaccine, vaccination_date);
}

// add infection by date
void init_infections_and_vaccinations(mio::abm::Model& model, mio::abm::PersonId pid,
                                      const std::vector<mio::abm::TimePoint>& infection_dates,
                                      const std::vector<mio::abm::TimePoint>& vaccination_dates)
{
    for (const auto& infection_date : infection_dates) {
        add_infection(model, pid, infection_date);
    }
    for (const auto& vaccination_date : vaccination_dates) {
        add_vaccination(model, pid, vaccination_date);
    }
}

/** @brief Initializes infections and vaccinations in the model based on the data file. Infections and vaccinations are sampled for each person in the model from the given data, depending on their age group.
 * @param[in] model The ABM model in which infections and vaccinations should be initialized.
 * @param[in] start_date The start date of the simulation. Infections and vaccinations that happened after this date are not initialized in the model.
 */
void initialize_infections_and_vaccinations(mio::abm::Model& model, const std::string& filename,
                                            const mio::abm::TimePoint start_date)
{
    if (!fs::exists(fs::path(filename))) {
        mio::log_error("Cannot read in data. File does not exist.");
    }
    // File pointer
    std::fstream fin;

    // Open an existing file
    fin.open(filename, std::ios::in);
    std::vector<int32_t> row;
    std::vector<std::string> row_string;
    std::string line;

    // Read the Titles from the Data file
    std::getline(fin, line);
    line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
    std::vector<std::string> titles;
    boost::split(titles, line, boost::is_any_of(","));
    uint32_t count_of_titles              = 0;
    std::map<std::string, uint32_t> index = {};
    for (auto const& title : titles) {
        index.insert({title, count_of_titles});
        row_string.push_back(title);
        count_of_titles++;
    }

    mio::CustomIndexArray<std::vector<std::vector<mio::abm::TimePoint>>, mio::AgeGroup> infections(
        mio::AgeGroup(num_age_groups), std::vector<std::vector<mio::abm::TimePoint>>{});
    mio::CustomIndexArray<std::vector<std::vector<mio::abm::TimePoint>>, mio::AgeGroup> vaccinations(
        mio::AgeGroup(num_age_groups), std::vector<std::vector<mio::abm::TimePoint>>{});

    // Read the data from the file and save sampled infections and vaccinations
    while (std::getline(fin, line)) {
        line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());
        std::vector<std::string> values;
        boost::split(values, line, boost::is_any_of(","));
        mio::AgeGroup age_group = get_age_group_from_string(values[index["age_group"]]);

        add_infection_dates_and_vaccinations(infections[age_group], vaccinations[age_group], values, index, start_date);
    }

    // Distribute persons into age groups
    mio::CustomIndexArray<std::vector<mio::abm::PersonId>, mio::AgeGroup> persons_by_age_group(
        mio::AgeGroup(num_age_groups), std::vector<mio::abm::PersonId>{});

    for (auto& person : model.get_persons()) {
        persons_by_age_group[person.get_age()].push_back(person.get_id());
    }

    // Randomly sample from the infections and vaccinations for each person in the model
    for (auto age_group : mio::AgeGroup(num_age_groups)) {
        if (persons_by_age_group[age_group].size() < infections[age_group].size()) {
            mio::log_warning("Fewer persons in age group {} than given in the sample file.", age_group.get());
        }
        if (infections[age_group].size() == 0) {
            continue; // no infections/vaccinations in this age group from data
        }
        for (mio::abm::PersonId pid : persons_by_age_group[age_group]) {
            // randomly select a sample
            auto rng                 = mio::RandomNumberGenerator();
            auto random_sample_index = mio::DiscreteDistribution<size_t>::get_instance()(
                rng, std::vector<double>(infections[age_group].size(), 1.0 / infections[age_group].size()));
            init_infections_and_vaccinations(model, pid, infections[age_group][random_sample_index],
                                             vaccinations[age_group][random_sample_index]);
        }
    }
}

/// An ABM setup taken from abm_minimal.cpp.
mio::abm::Model make_model(const std::string& infections_vaccinations_file, const mio::abm::TimePoint start_date,
                           const mio::RandomNumberGenerator& rng)
{
    // Create the model with 6 age groups.
    const int model_id = 15502; // Halle county id
    auto model         = mio::abm::Model(num_age_groups, model_id);
    model.get_rng()    = rng;

    // Set the age groups that can go to school; here this is AgeGroup(1) (i.e. 5-14)
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()                    = false;
    model.parameters.get<mio::abm::AgeGroupGotoSchool>()[age_group_5_to_14] = true;
    // Set the age groups that can go to work; here these are AgeGroup(2) and AgeGroup(3) (i.e. 15-34 and 35-59)
    model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple({age_group_15_to_34, age_group_35_to_59}, true);

    set_parameters(model.parameters);
    // Check if the parameters satisfy their constraints.
    model.parameters.check_constraints();

    //std::vector<size_t> age_distribution{1000, 1000, 1000, 1000, 1000, 1000}; // artificial age distribution
    std::vector<size_t> age_distribution_male{4788, 10890, 29272, 35381, 21757, 7503}; // Halle age distribution male
    std::vector<size_t> age_distribution_female{4439,  9571,  28188,
                                                34020, 28597, 12361}; // Halle age distribution female
    const size_t total_population =
        std::accumulate(age_distribution_male.begin(), age_distribution_male.end(), (size_t)0) +
        std::accumulate(age_distribution_female.begin(), age_distribution_female.end(), (size_t)0);

    auto random_household_member = mio::abm::HouseholdMember(num_age_groups);
    random_household_member.set_age_weight(age_group_0_to_4, (int)age_distribution_male[age_group_0_to_4.get()] +
                                                                 (int)age_distribution_female[age_group_0_to_4.get()]);
    random_household_member.set_age_weight(age_group_5_to_14,
                                           (int)age_distribution_male[age_group_5_to_14.get()] +
                                               (int)age_distribution_female[age_group_5_to_14.get()]);
    random_household_member.set_age_weight(age_group_15_to_34,
                                           (int)age_distribution_male[age_group_15_to_34.get()] +
                                               (int)age_distribution_female[age_group_15_to_34.get()]);
    random_household_member.set_age_weight(age_group_35_to_59,
                                           (int)age_distribution_male[age_group_35_to_59.get()] +
                                               (int)age_distribution_female[age_group_35_to_59.get()]);
    random_household_member.set_age_weight(age_group_60_to_79,
                                           (int)age_distribution_male[age_group_60_to_79.get()] +
                                               (int)age_distribution_female[age_group_60_to_79.get()]);
    random_household_member.set_age_weight(age_group_80_plus,
                                           (int)age_distribution_male[age_group_80_plus.get()] +
                                               (int)age_distribution_female[age_group_80_plus.get()]);

    // Add households with people drawn randomly from the age distribution.
    size_t persons_added = 0;

    while (total_population > persons_added) {
        auto household = mio::abm::Household();
        // random household size between 1 and 6 (weighted)
        auto hh_size = mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(),
                                                                         std::vector<ScalarType>{40, 20, 20, 10, 5, 5});
        household.add_members(random_household_member, (int)hh_size);
        add_household_to_model(model, household);
        persons_added += hh_size;
    }

    // Add one social event with 5 maximum contacts.
    // Maximum contacts limit the number of people that a person can infect while being at this location.
    auto event = model.add_location(mio::abm::LocationType::SocialEvent);
    model.get_location(event).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // Add hospital and ICU with 5 maximum contacs.
    auto hospital = model.add_location(mio::abm::LocationType::Hospital);
    model.get_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    auto icu = model.add_location(mio::abm::LocationType::ICU);
    model.get_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
    // Add one supermarket, maximum constacts are assumed to be 20.
    auto shop = model.add_location(mio::abm::LocationType::BasicsShop);
    model.get_location(shop).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // At every school, the maximum contacts are 20.
    auto school = model.add_location(mio::abm::LocationType::School);
    model.get_location(school).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
    // At every workplace, maximum contacts are 20.
    auto work = model.add_location(mio::abm::LocationType::Work);
    model.get_location(work).get_infection_parameters().set<mio::abm::MaximumContacts>(20);

    // People can get tested at work (and do this with 0.01 probability) from day 0 to day 20.
    auto validity_period       = mio::abm::days(2);
    auto probability           = 0.01;
    auto start_date_test       = start_date;
    auto end_date_test         = start_date + mio::abm::days(20);
    auto test_type             = mio::abm::TestType::Antigen;
    auto test_parameters       = model.parameters.get<mio::abm::TestData>()[test_type];
    auto testing_criteria_work = mio::abm::TestingCriteria();
    auto testing_scheme_work   = mio::abm::TestingScheme(testing_criteria_work, validity_period, start_date_test,
                                                         end_date_test, test_parameters, probability);
    model.get_testing_strategy().add_scheme(mio::abm::LocationType::Work, testing_scheme_work);

    // People test at home with 0.01 probability.
    probability              = 0.01;
    end_date_test            = start_date + mio::abm::days(100);
    auto testing_scheme_home = mio::abm::TestingScheme(testing_criteria_work, validity_period, start_date_test,
                                                       end_date_test, test_parameters, probability);
    model.get_testing_strategy().add_scheme(mio::abm::LocationType::Home, testing_scheme_home);

    // Add infections and vaccinations.
    initialize_infections_and_vaccinations(model, infections_vaccinations_file, start_date);

    // OLD:
    // Assign infection state to each person.
    // The infection states are chosen randomly with the following discrete distribution
    /*
    std::vector<ScalarType> infection_distribution{0.99, 0.002, 0.002, 0.002, 0.002, 0.002, 0.0, 0.0};
    for (auto& person : model.get_persons()) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
        auto person_rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(person_rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model.parameters, start_date, infection_state),
                                     person_rng, start_date, model.parameters);
        }
    }

    // Assign vaccinations to each person. The number of vaccinations is chosen randomly with the following discrete distribution
    std::vector<ScalarType> vaccination_distribution{0.5, 0.4, 0.1};
    // The date of the last vaccination is randomly chosen within a year before the start date of the simulation and six weeks back in case of two vaccinations.
    for (auto& person : model.get_persons()) {
        auto person_rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        size_t vaccination_count =
            mio::DiscreteDistribution<size_t>::get_instance()(person_rng, vaccination_distribution);
        if (vaccination_count > 0) {
            auto last_vaccination_date =
                start_date - mio::abm::days(365) +
                mio::abm::days(mio::UniformDistribution<ScalarType>::get_instance()(person_rng) * 365);
            if (vaccination_count > 1) {
                auto second_last_vaccination_date =
                    last_vaccination_date - mio::abm::days(42) +
                    mio::abm::days(mio::UniformDistribution<ScalarType>::get_instance()(person_rng) * 42);
                person.add_new_vaccination(mio::abm::ProtectionType::GenericVaccine, second_last_vaccination_date);
            }
            person.add_new_vaccination(mio::abm::ProtectionType::GenericVaccine, last_vaccination_date);
        }
    }
    */

    // Assign locations to the people
    for (auto& person : model.get_persons()) {
        const auto id = person.get_id();
        //assign shop and event
        model.assign_location(id, event);
        model.assign_location(id, shop);
        //assign hospital and ICU
        model.assign_location(id, hospital);
        model.assign_location(id, icu);
        //assign work/school to people depending on their age
        if (person.get_age() == age_group_5_to_14) {
            model.assign_location(id, school);
        }
        if (person.get_age() == age_group_15_to_34 || person.get_age() == age_group_35_to_59) {
            model.assign_location(id, work);
        }
    }

    // During the lockdown, social events are closed for 90% of people.
    auto t_lockdown = start_date + mio::abm::days(10);
    mio::abm::close_social_events(t_lockdown, 0.9, model.parameters);

    set_local_parameters(model);

    return model;
}

/**
 * @brief Logger to log the recovered persons for each age group.
 */
struct LogRecovered : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;

    static Type log(const mio::abm::Simulation<>& sim)
    {
        Type recovered{};
        Eigen::VectorXd age_group_counts =
            Eigen::VectorXd::Zero(Eigen::Index(sim.get_model().parameters.get_num_groups()));
        for (auto& person : sim.get_model().get_persons()) {
            if (person.get_infection_state(sim.get_time()) == mio::abm::InfectionState::Recovered) {
                age_group_counts[person.get_age().get()] += 1;
            }
        }
        recovered.first  = sim.get_time();
        recovered.second = age_group_counts;
        return recovered;
    }
};

/**
 * @brief Logger to log the number of persons for each age group that had PAIS.
 */
struct LogActivePAIS : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;

    static Type log(const mio::abm::Simulation<>& sim)
    {
        Type active_pais{};
        Eigen::VectorXd age_group_counts =
            Eigen::VectorXd::Zero(Eigen::Index(sim.get_model().parameters.get_num_groups()));
        for (auto& person : sim.get_model().get_persons()) {
            if (person.has_active_pais(sim.get_time())) {
                age_group_counts[person.get_age().get()] += 1;
            }
        }
        active_pais.first  = sim.get_time();
        active_pais.second = age_group_counts;
        return active_pais;
    }
};

/**
 * @brief Logger to log the recovered persons for each age group that were detected.
 */
struct LogRecoveredDetected : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;

    static Type log(const mio::abm::Simulation<>& sim)
    {
        Type recovered_detected{};
        Eigen::VectorXd age_group_counts =
            Eigen::VectorXd::Zero(Eigen::Index(sim.get_model().parameters.get_num_groups()));
        for (auto& person : sim.get_model().get_persons()) {
            if (person.get_infection_state(sim.get_time()) == mio::abm::InfectionState::Recovered &&
                person.get_infection().is_detected()) {
                age_group_counts[person.get_age().get()] += 1;
            }
        }
        recovered_detected.first  = sim.get_time();
        recovered_detected.second = age_group_counts;
        return recovered_detected;
    }
};

/**
 * @brief Logger to log vaccination counts for each age group for all agents in the simulation.
 */
struct LogVaccinationCounts : mio::LogAlways {
    using Type = std::pair<mio::abm::TimePoint, Eigen::VectorXd>;
    /** 
     * @brief Log the vaccination counts for each age group for all agents in the simulation.
     * @param[in] sim The simulation of the ABM.
     * @return A pair of the TimePoint and the TimeSeries of the vaccination counts for each age group for all agents in the simulation.
     * The TimeSeries is a CustomIndexArray with Index mio::AgeGroup and a value of a tuple, where each tuple contains the following information:
     * -# The count of people that have 0 vaccinations.
     * -# The count of people that have 1 vaccination.
     * -# The count of people that have 2 or more vaccinations.
    */
    static Type log(const mio::abm::Simulation<>& sim)
    {
        Type vaccination_counts{};
        Eigen::VectorXd age_group_counts = Eigen::VectorXd::Zero(
            Eigen::Index(sim.get_model().parameters.get_num_groups() *
                         static_cast<size_t>(mio::abm::VaccinationClass::Count))); // 3 values per age group

        for (auto& person : sim.get_model().get_persons()) {
            auto age_group = person.get_age().get();
            auto n_vaccinations =
                static_cast<size_t>(mio::abm::get_vaccination_class(person.get_vaccinations().size()));
            auto index = age_group + n_vaccinations * sim.get_model().parameters.get_num_groups();
            age_group_counts[index] += 1;
        }
        vaccination_counts.first  = sim.get_time();
        vaccination_counts.second = age_group_counts;
        return vaccination_counts;
    }
};

/**
 * @brief Main function to run the ABM simulation.
 */
int main()
{
    // mio::mpi::init();

    // mio::set_log_level(mio::LogLevel::warn);

    // Set start and end time for the simulation.
    auto t0 = get_time_point_from_string("2022-01-01"); // Start date is 2022-01-01 00:00:00.000

    auto tmax = t0 + mio::abm::days(20);

    const std::string infections_vaccinations_file =
        mio::path_join(mio::base_dir(), "AIMS_Halle_Data/260304_infections_vaccines_halle.csv");
    auto sim =
        mio::abm::Simulation(t0, std::move(make_model(infections_vaccinations_file, t0, mio::thread_local_rng())));

    // Create a history object to store the time series of the number of recovered (detected) persons and the number of vaccinations for each age group.
    mio::History<mio::abm::TimeSeriesWriter, LogRecoveredDetected> historyRecoveredDetected{
        Eigen::Index(sim.get_model().parameters.get_num_groups())};
    mio::History<mio::abm::TimeSeriesWriter, LogRecovered> historyRecovered{
        Eigen::Index(sim.get_model().parameters.get_num_groups())};
    mio::History<mio::abm::TimeSeriesWriter, LogVaccinationCounts> historyVaccinations{
        Eigen::Index(sim.get_model().parameters.get_num_groups() * 3)};
    mio::History<mio::abm::TimeSeriesWriter, LogActivePAIS> historyActivePAIS{
        Eigen::Index(sim.get_model().parameters.get_num_groups())};

    // Run the simulation until tmax with the history object.
    sim.advance(tmax, historyRecovered, historyVaccinations, historyActivePAIS);

    const std::string result_dir = mio::path_join(mio::base_dir(), "AIMS_Halle_Example_Results");
    if (!mio::create_directory(result_dir)) {
        mio::log_error("Could not create result directory \"{}\".", result_dir);
        return 1;
    }

    // The amount of recovered detected persons are written into the file "recovered.txt" as a table with 7 columns.
    // The first column is Time. The other columns correspond to the amount of people within each AgeGroup at this Time.
    std::ofstream outfile(mio::path_join(result_dir, "recovered.txt"));
    std::get<0>(historyRecovered.get_log()).print_table(outfile, {}, 16, 1, ',');
    std::cout << "Results written to recovered.txt" << std::endl;

    // The amount of vaccinated persons are written into the file "vaccinations.txt" as a table with 19 columns.
    // The first column is Time. The other columns correspond to the number of people within each AgeGroup and with the number of Vaccinations (0/1/2) at this Time.
    std::ofstream outfile2(mio::path_join(result_dir, "vaccinations.txt"));
    std::get<0>(historyVaccinations.get_log()).print_table(outfile2, {}, 16, 1, ',');
    std::cout << "Results written to vaccinations.txt" << std::endl;

    // The amount of active PAIS persons are written into the file "active_pais.txt" as a table with 7 columns.
    std::ofstream outfile3(mio::path_join(result_dir, "active_pais.txt"));
    std::get<0>(historyActivePAIS.get_log()).print_table(outfile3, {}, 16, 1, ',');
    std::cout << "Results written to active_pais.txt" << std::endl;

    // mio::mpi::finalize();

    return 0;
}
