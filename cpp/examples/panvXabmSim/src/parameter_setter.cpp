#include "../include/parameter_setter.h"
#include "../include/constants.h"
#include <cmath>

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
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_0_to_4}]   = 0.50;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_5_to_14}]  = 0.55;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_15_to_34}] = 0.60;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_35_to_59}] = 0.70;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_60_to_79}] = 0.83;
    params.get<mio::abm::SymptomsPerInfectedNoSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_80_plus}]  = 0.90;

    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_0_to_4}]   = 0.02;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_5_to_14}]  = 0.03;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_15_to_34}] = 0.04;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_35_to_59}] = 0.07;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_60_to_79}] = 0.17;
    params.get<mio::abm::SeverePerInfectedSymptoms>()[{mio::abm::VirusVariant::Alpha, age_group_80_plus}]  = 0.24;

    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_0_to_4}]   = 0.1;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_5_to_14}]  = 0.11;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_15_to_34}] = 0.12;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_35_to_59}] = 0.14;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_60_to_79}] = 0.33;
    params.get<mio::abm::CriticalPerInfectedSevere>()[{mio::abm::VirusVariant::Alpha, age_group_80_plus}]  = 0.62;

    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_0_to_4}]   = 0.12;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_5_to_14}]  = 0.13;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_15_to_34}] = 0.15;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_35_to_59}] = 0.26;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_60_to_79}] = 0.40;
    params.get<mio::abm::DeathsPerInfectedCritical>()[{mio::abm::VirusVariant::Alpha, age_group_80_plus}]  = 0.48;

    // Set infection parameters
    params.get<mio::abm::InfectionRateFromViralShed>()[{mio::abm::VirusVariant::Alpha}] = 1;
    params.get<mio::abm::AerosolTransmissionRates>()                                    = 0.0;
}

void set_local_parameters_ger(mio::abm::World& world)
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
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1.4; //17 hours //intensity
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 15.0; // Intensity
            break;
        case mio::abm::LocationType::School:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_school;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 4.8; //5h
            break;
        case mio::abm::LocationType::Work:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_work;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 3.0 * 0.5; // 7h
            break;
        case mio::abm::LocationType::SocialEvent:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1.2; //aufteilung
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 2.0; // intensity
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 6.0; // 4 hours
            break;
        case mio::abm::LocationType::BasicsShop:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.8; //aufteilung
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.33; // intensity
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 12.0; // 2 hours
            break;
        default:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_random;
            break;
        }
    }
}

void set_local_parameters_fra(mio::abm::World& world)
{
    const int n_age_groups = (int)world.parameters.get_num_groups();

    // setting this up in matrix-form would be much nicer,
    // but we somehow can't construct Eigen object with initializer lists
    /* baseline_home
        0,6881	0,6771	1,2965	0,9261	0,0337	0,0034
0,2257	1,6804	0,7570	1,4088	0,0235	0,0022
0,2563	0,3517	1,4941	0,7716	0,0381	0,0020
0,2096	0,6996	0,8293	1,2112	0,0640	0,0075
0,1710	0,4608	0,4983	0,6181	0,8598	0,0665
0,1990	0,6331	0,5450	1,1572	0,4189	0,2800
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_home(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_home[{age_group_0_to_4, age_group_0_to_4}]     = 0.6881;
    contacts_home[{age_group_0_to_4, age_group_5_to_14}]    = 0.6771;
    contacts_home[{age_group_0_to_4, age_group_15_to_34}]   = 1.2965;
    contacts_home[{age_group_0_to_4, age_group_35_to_59}]   = 0.9261;
    contacts_home[{age_group_0_to_4, age_group_60_to_79}]   = 0.0337;
    contacts_home[{age_group_0_to_4, age_group_80_plus}]    = 0.0034;
    contacts_home[{age_group_5_to_14, age_group_0_to_4}]    = 0.2257;
    contacts_home[{age_group_5_to_14, age_group_5_to_14}]   = 1.6804;
    contacts_home[{age_group_5_to_14, age_group_15_to_34}]  = 0.7570;
    contacts_home[{age_group_5_to_14, age_group_35_to_59}]  = 1.4088;
    contacts_home[{age_group_5_to_14, age_group_60_to_79}]  = 0.0235;
    contacts_home[{age_group_5_to_14, age_group_80_plus}]   = 0.0022;
    contacts_home[{age_group_15_to_34, age_group_0_to_4}]   = 0.2563;
    contacts_home[{age_group_15_to_34, age_group_5_to_14}]  = 0.3517;
    contacts_home[{age_group_15_to_34, age_group_15_to_34}] = 1.4941;
    contacts_home[{age_group_15_to_34, age_group_35_to_59}] = 0.7716;
    contacts_home[{age_group_15_to_34, age_group_60_to_79}] = 0.0381;
    contacts_home[{age_group_15_to_34, age_group_80_plus}]  = 0.0020;
    contacts_home[{age_group_35_to_59, age_group_0_to_4}]   = 0.2096;
    contacts_home[{age_group_35_to_59, age_group_5_to_14}]  = 0.6996;
    contacts_home[{age_group_35_to_59, age_group_15_to_34}] = 0.8293;
    contacts_home[{age_group_35_to_59, age_group_35_to_59}] = 1.2112;
    contacts_home[{age_group_35_to_59, age_group_60_to_79}] = 0.0640;
    contacts_home[{age_group_35_to_59, age_group_80_plus}]  = 0.0075;
    contacts_home[{age_group_60_to_79, age_group_0_to_4}]   = 0.1710;
    contacts_home[{age_group_60_to_79, age_group_5_to_14}]  = 0.4608;
    contacts_home[{age_group_60_to_79, age_group_15_to_34}] = 0.4983;
    contacts_home[{age_group_60_to_79, age_group_35_to_59}] = 0.6181;
    contacts_home[{age_group_60_to_79, age_group_60_to_79}] = 0.8598;
    contacts_home[{age_group_60_to_79, age_group_80_plus}]  = 0.0665;
    contacts_home[{age_group_80_plus, age_group_0_to_4}]    = 0.1990;
    contacts_home[{age_group_80_plus, age_group_5_to_14}]   = 0.6331;
    contacts_home[{age_group_80_plus, age_group_15_to_34}]  = 0.5450;
    contacts_home[{age_group_80_plus, age_group_35_to_59}]  = 1.1572;
    contacts_home[{age_group_80_plus, age_group_60_to_79}]  = 0.4189;
    contacts_home[{age_group_80_plus, age_group_80_plus}]   = 0.2800;

    /* baseline_school
        2,4233	0,3734	0,3248	0,3342	0,0032	0,0000
0,1972	3,3900	0,1830	0,2888	0,0054	0,0001
0,0317	0,4252	1,3594	0,1850	0,0051	0,0002
0,1477	0,4910	0,4562	0,2218	0,0061	0,0001
0,0207	0,0356	0,0564	0,0603	0,0277	0,0035
0,0000	0,0212	0,0289	0,0000	0,0000	0,0000
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_school(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_school[{age_group_0_to_4, age_group_0_to_4}]     = 2.4233;
    contacts_school[{age_group_0_to_4, age_group_5_to_14}]    = 0.3734;
    contacts_school[{age_group_0_to_4, age_group_15_to_34}]   = 0.3248;
    contacts_school[{age_group_0_to_4, age_group_35_to_59}]   = 0.3342;
    contacts_school[{age_group_0_to_4, age_group_60_to_79}]   = 0.0032;
    contacts_school[{age_group_0_to_4, age_group_80_plus}]    = 0.0000;
    contacts_school[{age_group_5_to_14, age_group_0_to_4}]    = 0.1972;
    contacts_school[{age_group_5_to_14, age_group_5_to_14}]   = 3.3900;
    contacts_school[{age_group_5_to_14, age_group_15_to_34}]  = 0.1830;
    contacts_school[{age_group_5_to_14, age_group_35_to_59}]  = 0.2888;
    contacts_school[{age_group_5_to_14, age_group_60_to_79}]  = 0.0054;
    contacts_school[{age_group_5_to_14, age_group_80_plus}]   = 0.0001;
    contacts_school[{age_group_15_to_34, age_group_0_to_4}]   = 0.0317;
    contacts_school[{age_group_15_to_34, age_group_5_to_14}]  = 0.4252;
    contacts_school[{age_group_15_to_34, age_group_15_to_34}] = 1.3594;
    contacts_school[{age_group_15_to_34, age_group_35_to_59}] = 0.1850;
    contacts_school[{age_group_15_to_34, age_group_60_to_79}] = 0.0051;
    contacts_school[{age_group_15_to_34, age_group_80_plus}]  = 0.0002;
    contacts_school[{age_group_35_to_59, age_group_0_to_4}]   = 0.1477;
    contacts_school[{age_group_35_to_59, age_group_5_to_14}]  = 0.4910;
    contacts_school[{age_group_35_to_59, age_group_15_to_34}] = 0.4562;
    contacts_school[{age_group_35_to_59, age_group_35_to_59}] = 0.2218;
    contacts_school[{age_group_35_to_59, age_group_60_to_79}] = 0.0061;
    contacts_school[{age_group_35_to_59, age_group_80_plus}]  = 0.0001;
    contacts_school[{age_group_60_to_79, age_group_0_to_4}]   = 0.0207;
    contacts_school[{age_group_60_to_79, age_group_5_to_14}]  = 0.0356;
    contacts_school[{age_group_60_to_79, age_group_15_to_34}] = 0.0564;
    contacts_school[{age_group_60_to_79, age_group_35_to_59}] = 0.0603;
    contacts_school[{age_group_60_to_79, age_group_60_to_79}] = 0.0277;
    contacts_school[{age_group_60_to_79, age_group_80_plus}]  = 0.0035;
    contacts_school[{age_group_80_plus, age_group_0_to_4}]    = 0.0000;
    contacts_school[{age_group_80_plus, age_group_5_to_14}]   = 0.0212;
    contacts_school[{age_group_80_plus, age_group_15_to_34}]  = 0.0289;
    contacts_school[{age_group_80_plus, age_group_35_to_59}]  = 0.0000;
    contacts_school[{age_group_80_plus, age_group_60_to_79}]  = 0.0000;
    contacts_school[{age_group_80_plus, age_group_80_plus}]   = 0.0000;

    /* baseline_work
        0,0000	0,0000	0,0000	0,0000	0,0000	0,0000
0,0000	0,0074	0,0222	0,0254	0,0000	0,0000
0,0000	0,0235	2,7163	2,4826	0,0050	0,0000
0,0000	0,0249	2,0297	3,5135	0,0052	0,0000
0,0000	0,0005	0,0133	0,0262	0,0001	0,0000
0,0000	0,0000	0,0001	0,0001	0,0000	0,0000
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
    contacts_work[{age_group_5_to_14, age_group_5_to_14}]   = 0.0074;
    contacts_work[{age_group_5_to_14, age_group_15_to_34}]  = 0.0222;
    contacts_work[{age_group_5_to_14, age_group_35_to_59}]  = 0.0254;
    contacts_work[{age_group_5_to_14, age_group_60_to_79}]  = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_80_plus}]   = 0.0000;
    contacts_work[{age_group_15_to_34, age_group_0_to_4}]   = 0.0000;
    contacts_work[{age_group_15_to_34, age_group_5_to_14}]  = 0.0235;
    contacts_work[{age_group_15_to_34, age_group_15_to_34}] = 2.7163;
    contacts_work[{age_group_15_to_34, age_group_35_to_59}] = 2.4826;
    contacts_work[{age_group_15_to_34, age_group_60_to_79}] = 0.0050;
    contacts_work[{age_group_15_to_34, age_group_80_plus}]  = 0.0000;
    contacts_work[{age_group_35_to_59, age_group_0_to_4}]   = 0.0000;
    contacts_work[{age_group_35_to_59, age_group_5_to_14}]  = 0.0249;
    contacts_work[{age_group_35_to_59, age_group_15_to_34}] = 2.0297;
    contacts_work[{age_group_35_to_59, age_group_35_to_59}] = 3.5135;
    contacts_work[{age_group_35_to_59, age_group_60_to_79}] = 0.0052;
    contacts_work[{age_group_35_to_59, age_group_80_plus}]  = 0.0000;
    contacts_work[{age_group_60_to_79, age_group_0_to_4}]   = 0.0000;
    contacts_work[{age_group_60_to_79, age_group_5_to_14}]  = 0.0005;
    contacts_work[{age_group_60_to_79, age_group_15_to_34}] = 0.0133;
    contacts_work[{age_group_60_to_79, age_group_35_to_59}] = 0.0262;
    contacts_work[{age_group_60_to_79, age_group_60_to_79}] = 0.0001;
    contacts_work[{age_group_60_to_79, age_group_80_plus}]  = 0.0000;
    contacts_work[{age_group_80_plus, age_group_0_to_4}]    = 0.0000;
    contacts_work[{age_group_80_plus, age_group_5_to_14}]   = 0.0000;
    contacts_work[{age_group_80_plus, age_group_15_to_34}]  = 0.0001;
    contacts_work[{age_group_80_plus, age_group_35_to_59}]  = 0.0001;
    contacts_work[{age_group_80_plus, age_group_60_to_79}]  = 0.0000;
    contacts_work[{age_group_80_plus, age_group_80_plus}]   = 0.0000;

    /* baseline_other
        0,6935	0,4678	0,9381	1,1670	0,3841	0,0317
0,2217	2,2229	0,8144	1,0773	0,2950	0,0378
0,1051	0,4410	2,7757	1,2325	0,1793	0,0230
0,0847	0,2358	1,1066	1,8774	0,4578	0,0398
0,0495	0,1284	0,7258	1,5573	1,0805	0,1023
0,0446	0,1150	0,4335	1,0023	0,8198	0,1387
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_other(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_other[{age_group_0_to_4, age_group_0_to_4}]     = 0.6935;
    contacts_other[{age_group_0_to_4, age_group_5_to_14}]    = 0.4678;
    contacts_other[{age_group_0_to_4, age_group_15_to_34}]   = 0.9381;
    contacts_other[{age_group_0_to_4, age_group_35_to_59}]   = 1.1670;
    contacts_other[{age_group_0_to_4, age_group_60_to_79}]   = 0.3841;
    contacts_other[{age_group_0_to_4, age_group_80_plus}]    = 0.0317;
    contacts_other[{age_group_5_to_14, age_group_0_to_4}]    = 0.2217;
    contacts_other[{age_group_5_to_14, age_group_5_to_14}]   = 2.2229;
    contacts_other[{age_group_5_to_14, age_group_15_to_34}]  = 0.8144;
    contacts_other[{age_group_5_to_14, age_group_35_to_59}]  = 1.0773;
    contacts_other[{age_group_5_to_14, age_group_60_to_79}]  = 0.2950;
    contacts_other[{age_group_5_to_14, age_group_80_plus}]   = 0.0378;
    contacts_other[{age_group_15_to_34, age_group_0_to_4}]   = 0.1051;
    contacts_other[{age_group_15_to_34, age_group_5_to_14}]  = 0.4410;
    contacts_other[{age_group_15_to_34, age_group_15_to_34}] = 2.7757;
    contacts_other[{age_group_15_to_34, age_group_35_to_59}] = 1.2325;
    contacts_other[{age_group_15_to_34, age_group_60_to_79}] = 0.1793;
    contacts_other[{age_group_15_to_34, age_group_80_plus}]  = 0.0230;
    contacts_other[{age_group_35_to_59, age_group_0_to_4}]   = 0.0847;
    contacts_other[{age_group_35_to_59, age_group_5_to_14}]  = 0.2358;
    contacts_other[{age_group_35_to_59, age_group_15_to_34}] = 1.1066;
    contacts_other[{age_group_35_to_59, age_group_35_to_59}] = 1.8774;
    contacts_other[{age_group_35_to_59, age_group_60_to_79}] = 0.4578;
    contacts_other[{age_group_35_to_59, age_group_80_plus}]  = 0.0398;
    contacts_other[{age_group_60_to_79, age_group_0_to_4}]   = 0.0495;
    contacts_other[{age_group_60_to_79, age_group_5_to_14}]  = 0.1284;
    contacts_other[{age_group_60_to_79, age_group_15_to_34}] = 0.7258;
    contacts_other[{age_group_60_to_79, age_group_35_to_59}] = 1.5573;
    contacts_other[{age_group_60_to_79, age_group_60_to_79}] = 1.0805;
    contacts_other[{age_group_60_to_79, age_group_80_plus}]  = 0.1023;
    contacts_other[{age_group_80_plus, age_group_0_to_4}]    = 0.0446;
    contacts_other[{age_group_80_plus, age_group_5_to_14}]   = 0.1150;
    contacts_other[{age_group_80_plus, age_group_15_to_34}]  = 0.4335;
    contacts_other[{age_group_80_plus, age_group_35_to_59}]  = 1.0023;
    contacts_other[{age_group_80_plus, age_group_60_to_79}]  = 0.8198;
    contacts_other[{age_group_80_plus, age_group_80_plus}]   = 0.1387;

    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_random(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 1.0);

    for (auto& loc : world.get_locations()) {
        switch (loc.get_type()) {
        case mio::abm::LocationType::Home:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_home;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1.4; //17 hours //intensity
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 15.0; // Intensity
            break;
        case mio::abm::LocationType::School:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_school;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 4.8; //5h
            break;
        case mio::abm::LocationType::Work:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_work;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 3.0 * 0.5; // 7h
            break;
        case mio::abm::LocationType::SocialEvent:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1.2; //aufteilung
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 2.0; // intensity
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 6.0; // 4 hours
            break;
        case mio::abm::LocationType::BasicsShop:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.8; //aufteilung
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.33; // intensity
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 12.0; // 2 hours
            break;
        default:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_random;
            break;
        }
    }
}

void set_local_parameters_usa(mio::abm::World& world)
{
    const int n_age_groups = (int)world.parameters.get_num_groups();

    // setting this up in matrix-form would be much nicer,
    // but we somehow can't construct Eigen object with initializer lists
    /* baseline_home
        0,6197	0,7925	1,1687	0,9485	0,0335	0,0040
        0,2460	1,7677	0,8041	1,3347	0,0308	0,0032
        0,2832	0,4564	1,5297	0,7018	0,0449	0,0028
        0,2613	0,8550	0,9010	1,1937	0,0834	0,0095
        0,2117	0,5553	0,5944	0,7254	0,6861	0,0641
        0,1787	0,5972	0,5194	0,9997	0,2481	0,2714
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_home(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_home[{age_group_0_to_4, age_group_0_to_4}]     = 0.6197;
    contacts_home[{age_group_0_to_4, age_group_5_to_14}]    = 0.7925;
    contacts_home[{age_group_0_to_4, age_group_15_to_34}]   = 1.1687;
    contacts_home[{age_group_0_to_4, age_group_35_to_59}]   = 0.9485;
    contacts_home[{age_group_0_to_4, age_group_60_to_79}]   = 0.0335;
    contacts_home[{age_group_0_to_4, age_group_80_plus}]    = 0.0040;
    contacts_home[{age_group_5_to_14, age_group_0_to_4}]    = 0.2460;
    contacts_home[{age_group_5_to_14, age_group_5_to_14}]   = 1.7677;
    contacts_home[{age_group_5_to_14, age_group_15_to_34}]  = 0.8041;
    contacts_home[{age_group_5_to_14, age_group_35_to_59}]  = 1.3347;
    contacts_home[{age_group_5_to_14, age_group_60_to_79}]  = 0.0308;
    contacts_home[{age_group_5_to_14, age_group_80_plus}]   = 0.0032;
    contacts_home[{age_group_15_to_34, age_group_0_to_4}]   = 0.2832;
    contacts_home[{age_group_15_to_34, age_group_5_to_14}]  = 0.4564;
    contacts_home[{age_group_15_to_34, age_group_15_to_34}] = 1.5297;
    contacts_home[{age_group_15_to_34, age_group_35_to_59}] = 0.7018;
    contacts_home[{age_group_15_to_34, age_group_60_to_79}] = 0.0449;
    contacts_home[{age_group_15_to_34, age_group_80_plus}]  = 0.0028;
    contacts_home[{age_group_35_to_59, age_group_0_to_4}]   = 0.2613;
    contacts_home[{age_group_35_to_59, age_group_5_to_14}]  = 0.8550;
    contacts_home[{age_group_35_to_59, age_group_15_to_34}] = 0.9010;
    contacts_home[{age_group_35_to_59, age_group_35_to_59}] = 1.1937;
    contacts_home[{age_group_35_to_59, age_group_60_to_79}] = 0.0834;
    contacts_home[{age_group_35_to_59, age_group_80_plus}]  = 0.0095;
    contacts_home[{age_group_60_to_79, age_group_0_to_4}]   = 0.2117;
    contacts_home[{age_group_60_to_79, age_group_5_to_14}]  = 0.5553;
    contacts_home[{age_group_60_to_79, age_group_15_to_34}] = 0.5944;
    contacts_home[{age_group_60_to_79, age_group_35_to_59}] = 0.7254;
    contacts_home[{age_group_60_to_79, age_group_60_to_79}] = 0.6861;
    contacts_home[{age_group_60_to_79, age_group_80_plus}]  = 0.0641;
    contacts_home[{age_group_80_plus, age_group_0_to_4}]    = 0.1787;
    contacts_home[{age_group_80_plus, age_group_5_to_14}]   = 0.5972;
    contacts_home[{age_group_80_plus, age_group_15_to_34}]  = 0.5194;
    contacts_home[{age_group_80_plus, age_group_35_to_59}]  = 0.9997;
    contacts_home[{age_group_80_plus, age_group_60_to_79}]  = 0.2481;
    contacts_home[{age_group_80_plus, age_group_80_plus}]   = 0.2714;

    /* baseline_school
        1,1966	0,2696	0,2404	0,2202	0,0036	0,0000
        0,1381	3,9843	0,2120	0,2775	0,0077	0,0001
        0,0240	0,5101	1,8745	0,1804	0,0066	0,0003
        0,0982	0,4768	0,4823	0,1892	0,0071	0,0001
        0,0251	0,0527	0,0906	0,0676	0,0291	0,0044
        0,0000	0,0211	0,0289	0,0000	0,0000	0,0000
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_school(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_school[{age_group_0_to_4, age_group_0_to_4}]     = 1.1966;
    contacts_school[{age_group_0_to_4, age_group_5_to_14}]    = 0.2696;
    contacts_school[{age_group_0_to_4, age_group_15_to_34}]   = 0.2404;
    contacts_school[{age_group_0_to_4, age_group_35_to_59}]   = 0.2202;
    contacts_school[{age_group_0_to_4, age_group_60_to_79}]   = 0.0036;
    contacts_school[{age_group_0_to_4, age_group_80_plus}]    = 0.0000;
    contacts_school[{age_group_5_to_14, age_group_0_to_4}]    = 0.1381;
    contacts_school[{age_group_5_to_14, age_group_5_to_14}]   = 3.9843;
    contacts_school[{age_group_5_to_14, age_group_15_to_34}]  = 0.2120;
    contacts_school[{age_group_5_to_14, age_group_35_to_59}]  = 0.2775;
    contacts_school[{age_group_5_to_14, age_group_60_to_79}]  = 0.0077;
    contacts_school[{age_group_5_to_14, age_group_80_plus}]   = 0.0001;
    contacts_school[{age_group_15_to_34, age_group_0_to_4}]   = 0.0240;
    contacts_school[{age_group_15_to_34, age_group_5_to_14}]  = 0.5101;
    contacts_school[{age_group_15_to_34, age_group_15_to_34}] = 1.8745;
    contacts_school[{age_group_15_to_34, age_group_35_to_59}] = 0.1804;
    contacts_school[{age_group_15_to_34, age_group_60_to_79}] = 0.0066;
    contacts_school[{age_group_15_to_34, age_group_80_plus}]  = 0.0003;
    contacts_school[{age_group_35_to_59, age_group_0_to_4}]   = 0.0982;
    contacts_school[{age_group_35_to_59, age_group_5_to_14}]  = 0.4768;
    contacts_school[{age_group_35_to_59, age_group_15_to_34}] = 0.4823;
    contacts_school[{age_group_35_to_59, age_group_35_to_59}] = 0.1892;
    contacts_school[{age_group_35_to_59, age_group_60_to_79}] = 0.0071;
    contacts_school[{age_group_35_to_59, age_group_80_plus}]  = 0.0001;
    contacts_school[{age_group_60_to_79, age_group_0_to_4}]   = 0.0251;
    contacts_school[{age_group_60_to_79, age_group_5_to_14}]  = 0.0527;
    contacts_school[{age_group_60_to_79, age_group_15_to_34}] = 0.0906;
    contacts_school[{age_group_60_to_79, age_group_35_to_59}] = 0.0676;
    contacts_school[{age_group_60_to_79, age_group_60_to_79}] = 0.0291;
    contacts_school[{age_group_60_to_79, age_group_80_plus}]  = 0.0044;
    contacts_school[{age_group_80_plus, age_group_0_to_4}]    = 0.0000;
    contacts_school[{age_group_80_plus, age_group_5_to_14}]   = 0.0211;
    contacts_school[{age_group_80_plus, age_group_15_to_34}]  = 0.0289;
    contacts_school[{age_group_80_plus, age_group_35_to_59}]  = 0.0000;
    contacts_school[{age_group_80_plus, age_group_60_to_79}]  = 0.0000;
    contacts_school[{age_group_80_plus, age_group_80_plus}]   = 0.0000;

    /* baseline_work
        0,0000	0,0000	0,0000	0,0000	0,0000	0,0000
        0,0000	0,0402	0,0502	0,0548	0,0000	0,0000
        0,0000	0,0526	2,5209	2,3056	0,0387	0,0000
        0,0000	0,0664	1,9258	3,4223	0,0442	0,0000
        0,0000	0,0105	0,1102	0,2350	0,0047	0,0000
        0,0000	0,0000	0,0001	0,0001	0,0000	0,0000
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
    contacts_work[{age_group_5_to_14, age_group_5_to_14}]   = 0.0402;
    contacts_work[{age_group_5_to_14, age_group_15_to_34}]  = 0.0502;
    contacts_work[{age_group_5_to_14, age_group_35_to_59}]  = 0.0548;
    contacts_work[{age_group_5_to_14, age_group_60_to_79}]  = 0.0000;
    contacts_work[{age_group_5_to_14, age_group_80_plus}]   = 0.0000;
    contacts_work[{age_group_15_to_34, age_group_0_to_4}]   = 0.0000;
    contacts_work[{age_group_15_to_34, age_group_5_to_14}]  = 0.0526;
    contacts_work[{age_group_15_to_34, age_group_15_to_34}] = 2.5209;
    contacts_work[{age_group_15_to_34, age_group_35_to_59}] = 2.3056;
    contacts_work[{age_group_15_to_34, age_group_60_to_79}] = 0.0387;
    contacts_work[{age_group_15_to_34, age_group_80_plus}]  = 0.0000;
    contacts_work[{age_group_35_to_59, age_group_0_to_4}]   = 0.0000;
    contacts_work[{age_group_35_to_59, age_group_5_to_14}]  = 0.0664;
    contacts_work[{age_group_35_to_59, age_group_15_to_34}] = 1.9258;
    contacts_work[{age_group_35_to_59, age_group_35_to_59}] = 3.4223;
    contacts_work[{age_group_35_to_59, age_group_60_to_79}] = 0.0442;
    contacts_work[{age_group_35_to_59, age_group_80_plus}]  = 0.0000;
    contacts_work[{age_group_60_to_79, age_group_0_to_4}]   = 0.0000;
    contacts_work[{age_group_60_to_79, age_group_5_to_14}]  = 0.0105;
    contacts_work[{age_group_60_to_79, age_group_15_to_34}] = 0.1102;
    contacts_work[{age_group_60_to_79, age_group_35_to_59}] = 0.2350;
    contacts_work[{age_group_60_to_79, age_group_60_to_79}] = 0.0047;
    contacts_work[{age_group_60_to_79, age_group_80_plus}]  = 0.0000;
    contacts_work[{age_group_80_plus, age_group_0_to_4}]    = 0.0000;
    contacts_work[{age_group_80_plus, age_group_5_to_14}]   = 0.0000;
    contacts_work[{age_group_80_plus, age_group_15_to_34}]  = 0.0001;
    contacts_work[{age_group_80_plus, age_group_35_to_59}]  = 0.0001;
    contacts_work[{age_group_80_plus, age_group_60_to_79}]  = 0.0000;
    contacts_work[{age_group_80_plus, age_group_80_plus}]   = 0.0000;

    /* baseline_other
        0,7819	0,5386	1,0918	1,2430	0,3503	0,0309
        0,2522	2,6569	0,9886	1,1725	0,2698	0,0380
        0,1245	0,5258	3,4633	1,3746	0,1660	0,0232
        0,0913	0,2606	1,2291	1,9120	0,3962	0,0368
        0,0469	0,1202	0,7026	1,3731	0,7760	0,0768
        0,0315	0,0835	0,3139	0,6754	0,4473	0,0847
    */
    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_other(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 0.);
    contacts_other[{age_group_0_to_4, age_group_0_to_4}]     = 0.7819;
    contacts_other[{age_group_0_to_4, age_group_5_to_14}]    = 0.5386;
    contacts_other[{age_group_0_to_4, age_group_15_to_34}]   = 1.0918;
    contacts_other[{age_group_0_to_4, age_group_35_to_59}]   = 1.2430;
    contacts_other[{age_group_0_to_4, age_group_60_to_79}]   = 0.3503;
    contacts_other[{age_group_0_to_4, age_group_80_plus}]    = 0.0309;
    contacts_other[{age_group_5_to_14, age_group_0_to_4}]    = 0.2522;
    contacts_other[{age_group_5_to_14, age_group_5_to_14}]   = 2.6569;
    contacts_other[{age_group_5_to_14, age_group_15_to_34}]  = 0.9886;
    contacts_other[{age_group_5_to_14, age_group_35_to_59}]  = 1.1725;
    contacts_other[{age_group_5_to_14, age_group_60_to_79}]  = 0.2698;
    contacts_other[{age_group_5_to_14, age_group_80_plus}]   = 0.0380;
    contacts_other[{age_group_15_to_34, age_group_0_to_4}]   = 0.1245;
    contacts_other[{age_group_15_to_34, age_group_5_to_14}]  = 0.5258;
    contacts_other[{age_group_15_to_34, age_group_15_to_34}] = 3.4633;
    contacts_other[{age_group_15_to_34, age_group_35_to_59}] = 1.3746;
    contacts_other[{age_group_15_to_34, age_group_60_to_79}] = 0.1660;
    contacts_other[{age_group_15_to_34, age_group_80_plus}]  = 0.0232;
    contacts_other[{age_group_35_to_59, age_group_0_to_4}]   = 0.0913;
    contacts_other[{age_group_35_to_59, age_group_5_to_14}]  = 0.2606;
    contacts_other[{age_group_35_to_59, age_group_15_to_34}] = 1.2291;
    contacts_other[{age_group_35_to_59, age_group_35_to_59}] = 1.9120;
    contacts_other[{age_group_35_to_59, age_group_60_to_79}] = 0.3962;
    contacts_other[{age_group_35_to_59, age_group_80_plus}]  = 0.0368;
    contacts_other[{age_group_60_to_79, age_group_0_to_4}]   = 0.0469;
    contacts_other[{age_group_60_to_79, age_group_5_to_14}]  = 0.1202;
    contacts_other[{age_group_60_to_79, age_group_15_to_34}] = 0.7026;
    contacts_other[{age_group_60_to_79, age_group_35_to_59}] = 1.3731;
    contacts_other[{age_group_60_to_79, age_group_60_to_79}] = 0.7760;
    contacts_other[{age_group_60_to_79, age_group_80_plus}]  = 0.0768;
    contacts_other[{age_group_80_plus, age_group_0_to_4}]    = 0.0315;
    contacts_other[{age_group_80_plus, age_group_5_to_14}]   = 0.0835;
    contacts_other[{age_group_80_plus, age_group_15_to_34}]  = 0.3139;
    contacts_other[{age_group_80_plus, age_group_35_to_59}]  = 0.6754;
    contacts_other[{age_group_80_plus, age_group_60_to_79}]  = 0.4473;
    contacts_other[{age_group_80_plus, age_group_80_plus}]   = 0.0847;

    mio::CustomIndexArray<ScalarType, mio::AgeGroup, mio::AgeGroup> contacts_random(
        {mio::AgeGroup(n_age_groups), mio::AgeGroup(n_age_groups)}, 1.0);

    for (auto& loc : world.get_locations()) {
        switch (loc.get_type()) {
        case mio::abm::LocationType::Home:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_home;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1.4; //17 hours //intensity
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 15.0; // Intensity
            break;
        case mio::abm::LocationType::School:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_school;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 4.8; //5h
            break;
        case mio::abm::LocationType::Work:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_work;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 3.0 * 0.5; // 7h
            break;
        case mio::abm::LocationType::SocialEvent:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 1.2; //aufteilung
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 2.0; // intensity
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 6.0; // 4 hours
            break;
        case mio::abm::LocationType::BasicsShop:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.8; //aufteilung
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 0.33; // intensity
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= 12.0; // 2 hours
            break;
        default:
            loc.get_infection_parameters().get<mio::abm::ContactRates>() = contacts_random;
            break;
        }
    }
}

void set_local_parameters_event(mio::abm::World& world, double contact_rate_multiplier)
{
    set_local_parameters_ger(world);
    for (auto& loc : world.get_locations()) {
        switch (loc.get_type()) {
        case mio::abm::LocationType::Home:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= contact_rate_multiplier; //15 hours
            break;
        case mio::abm::LocationType::School:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= contact_rate_multiplier; //2 hours
            break;
        case mio::abm::LocationType::Work:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= contact_rate_multiplier; // 3 hours
            break;
        case mio::abm::LocationType::SocialEvent:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= contact_rate_multiplier; // 3 hours
            break;
        case mio::abm::LocationType::BasicsShop:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= contact_rate_multiplier; // 2 hours
            break;
        default:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().array() *= contact_rate_multiplier;
            break;
        }
    }
}