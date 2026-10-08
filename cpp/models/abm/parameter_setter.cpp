#include "parameter_setter.h"
#include "memilio/config.h"
#include <cassert>
#include <cmath>
#include <cstdlib>

// Define age group constants
const mio::AgeGroup age_group_0_to_4{0};
const mio::AgeGroup age_group_5_to_14{1};
const mio::AgeGroup age_group_15_to_34{2};
const mio::AgeGroup age_group_35_to_59{3};
const mio::AgeGroup age_group_60_to_79{4};
const mio::AgeGroup age_group_80_plus{5};

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
    auto incubation_period_my_sigma = get_my_and_sigma({4.5, 1.5});
    params.get<mio::abm::TimeExposedToNoSymptoms>() =
        mio::ParameterDistributionLogNormal{incubation_period_my_sigma.first, incubation_period_my_sigma.second};

    auto InfectedNoSymptoms_to_symptoms_my_sigma             = get_my_and_sigma({1.1, 0.9});
    params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>() = mio::ParameterDistributionLogNormal{
        InfectedNoSymptoms_to_symptoms_my_sigma.first, InfectedNoSymptoms_to_symptoms_my_sigma.second};

    auto TimeInfectedNoSymptomsToRecovered_my_sigma           = get_my_and_sigma({8.0, 2.0});
    params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>() = mio::ParameterDistributionLogNormal{
        TimeInfectedNoSymptomsToRecovered_my_sigma.first, TimeInfectedNoSymptomsToRecovered_my_sigma.second};

    auto TimeInfectedSymptomsToSevere_my_sigma           = get_my_and_sigma({6.6, 4.9});
    params.get<mio::abm::TimeInfectedSymptomsToSevere>() = mio::ParameterDistributionLogNormal{
        TimeInfectedSymptomsToSevere_my_sigma.first, TimeInfectedSymptomsToSevere_my_sigma.second};

    auto TimeInfectedSymptomsToRecovered_my_sigma           = get_my_and_sigma({8.0, 2.0});
    params.get<mio::abm::TimeInfectedSymptomsToRecovered>() = mio::ParameterDistributionLogNormal{
        TimeInfectedSymptomsToRecovered_my_sigma.first, TimeInfectedSymptomsToRecovered_my_sigma.second};

    auto TimeInfectedSevereToCritical_my_sigma           = get_my_and_sigma({1.5, 2.0});
    params.get<mio::abm::TimeInfectedSevereToCritical>() = mio::ParameterDistributionLogNormal{
        TimeInfectedSevereToCritical_my_sigma.first, TimeInfectedSevereToCritical_my_sigma.second};

    auto TimeInfectedSevereToRecovered_my_sigma           = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedSevereToRecovered>() = mio::ParameterDistributionLogNormal{
        TimeInfectedSevereToRecovered_my_sigma.first, TimeInfectedSevereToRecovered_my_sigma.second};

    auto TimeInfectedCriticalToDead_my_sigma           = get_my_and_sigma({10.7, 4.8});
    params.get<mio::abm::TimeInfectedCriticalToDead>() = mio::ParameterDistributionLogNormal{
        TimeInfectedCriticalToDead_my_sigma.first, TimeInfectedCriticalToDead_my_sigma.second};

    auto TimeInfectedCriticalToRecovered_my_sigma           = get_my_and_sigma({18.1, 6.3});
    params.get<mio::abm::TimeInfectedCriticalToRecovered>() = mio::ParameterDistributionLogNormal{
        TimeInfectedCriticalToRecovered_my_sigma.first, TimeInfectedCriticalToRecovered_my_sigma.second};

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
    params.get<mio::abm::AerosolTransmissionRates>() = 0.0;
}

void set_local_parameters(mio::abm::Model& world)
{
    const int n_age_groups = (int)world.parameters.get_num_groups();

    // baseline_home
    Eigen::MatrixX<ScalarType> contacts_home{
        {0.4413, 0.0504, 1.2383, 0.8033, 0.0494, 0.0017}, // age_group_0_to_4   -> X
        {0.0485, 0.7616, 0.6532, 1.1614, 0.0256, 0.0013}, // age_group_5_to_14  -> X
        {0.1800, 0.1795, 0.8806, 0.6413, 0.0429, 0.0032}, // age_group_15_to_34 -> X
        {0.0495, 0.2639, 0.5189, 0.8277, 0.0679, 0.0014}, // age_group_35_to_59 -> X
        {0.0087, 0.0394, 0.1417, 0.3834, 0.7064, 0.0447}, // age_group_60_to_79 -> X
        {0.0292, 0.0648, 0.1248, 0.4179, 0.3497, 0.1544} //  age_group_80_plus  -> X
    };

    // baseline_school
    Eigen::MatrixX<ScalarType> contacts_school{
        {1.1165, 0.2741, 0.2235, 0.1028, 0.0007, 0.0000}, // age_group_0_to_4   -> X
        {0.1627, 1.9412, 0.2431, 0.1780, 0.0130, 0.0000}, // age_group_5_to_14  -> X
        {0.0148, 0.1646, 1.1266, 0.0923, 0.0074, 0.0000}, // age_group_15_to_34 -> X
        {0.0367, 0.1843, 0.3265, 0.0502, 0.0021, 0.0005}, // age_group_35_to_59 -> X
        {0.0004, 0.0370, 0.0115, 0.0014, 0.0039, 0.0000}, // age_group_60_to_79 -> X
        {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000} //  age_group_80_plus  -> X
    };

    // baseline_work
    Eigen::MatrixX<ScalarType> contacts_work{
        {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, // age_group_0_to_4   -> X
        {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000}, // age_group_5_to_14  -> X
        {0.0000, 0.0127, 1.7570, 1.6050, 0.0133, 0.0000}, // age_group_15_to_34 -> X
        {0.0000, 0.0020, 1.0311, 2.3166, 0.0098, 0.0000}, // age_group_35_to_59 -> X
        {0.0000, 0.0002, 0.0194, 0.0325, 0.0003, 0.0000}, // age_group_60_to_79 -> X
        {0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000} //  age_group_80_plus  -> X
    };

    // baseline_other
    Eigen::MatrixX<ScalarType> contacts_other{
        {0.5170, 0.3997, 0.7957, 0.9958, 0.3239, 0.0428}, // age_group_0_to_4   -> X
        {0.0632, 0.9121, 0.3254, 0.4731, 0.2355, 0.0148}, // age_group_5_to_14  -> X
        {0.0336, 0.1604, 1.7529, 0.8622, 0.1440, 0.0077}, // age_group_15_to_34 -> X
        {0.0204, 0.1444, 0.5738, 1.2127, 0.3433, 0.0178}, // age_group_35_to_59 -> X
        {0.0371, 0.0393, 0.4171, 0.9666, 0.7495, 0.0257}, // age_group_60_to_79 -> X
        {0.0791, 0.0800, 0.3480, 0.5588, 0.2769, 0.0180} //  age_group_80_plus  -> X
    };

    auto contacts_random = Eigen::MatrixX<ScalarType>::Ones(n_age_groups, n_age_groups);

    for (auto& loc : world.get_locations()) {
        switch (loc.get_type()) {
        case mio::abm::LocationType::Home:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() = contacts_home;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 1.4; //17 hours //intensity
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 15.0; // Intensity
            break;
        case mio::abm::LocationType::School:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() = contacts_school;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 4.8; //5h
            break;
        case mio::abm::LocationType::Work:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() = contacts_work;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 3.0 * 0.5; // 7h
            break;
        case mio::abm::LocationType::SocialEvent:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 1.2; //aufteilung
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 2.0; // intensity
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 6.0; // 4 hours
            break;
        case mio::abm::LocationType::BasicsShop:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() = contacts_other;
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 0.8; //aufteilung
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 0.33; // intensity
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= 12.0; // 2 hours
            break;
        default:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() = contacts_random;
            break;
        }
    }
}

void set_local_parameters_event(mio::abm::Model& world, double contact_rate_multiplier)
{
    set_local_parameters(world);
    for (auto& loc : world.get_locations()) {
        switch (loc.get_type()) {
        case mio::abm::LocationType::Home:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *=
                contact_rate_multiplier; //15 hours
            break;
        case mio::abm::LocationType::School:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *=
                contact_rate_multiplier; //2 hours
            break;
        case mio::abm::LocationType::Work:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *=
                contact_rate_multiplier; // 3 hours
            break;
        case mio::abm::LocationType::SocialEvent:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *=
                contact_rate_multiplier; // 3 hours
            break;
        case mio::abm::LocationType::BasicsShop:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *=
                contact_rate_multiplier; // 2 hours
            break;
        default:
            loc.get_infection_parameters().get<mio::abm::ContactRates>().get_baseline() *= contact_rate_multiplier;
            break;
        }
    }
}
