#pragma once

#include <iostream>
#include <filesystem>
#include <string>
#include <vector>
#include <cstddef>
#include <algorithm>

#include "boost/outcome/try.hpp"
#include "boost/outcome/result.hpp"
#include "boost/optional.hpp"
#include "boost/filesystem.hpp"

#include "memilio/utils/date.h"
#include "memilio/io/mobility_io.h"
#include "models/ode_secirvvs/model.h"

namespace fs = boost::filesystem;

enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count,
};

/**
 * different types of NPI, used as DampingType.
 */
enum class Intervention
{
    Home,
    SchoolClosure,
    HomeOffice,
    GatheringBanFacilitiesClosure,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Count,
};

static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}
                                                                        };

template <typename FP, class Tag>
void set_parameters(mio::osecirvvs::Parameters<FP>& params, const std::vector<FP>& parameter_values)
{
    std::copy(parameter_values.begin(), parameter_values.end(), params.template get<Tag>());
}

template <typename FP, class Tag>
void set_parameters(mio::osecirvvs::Parameters<FP>& params, double parameter_values)
{
    std::fill(params.template get<Tag>().begin(), params.template get<Tag>().end(), FP(parameter_values));
}

class OptimizationModel
{
public:
    OptimizationModel(fs::path data_dir, double t0, double tmax, double num_age_groups);

    const fs::path& data_directory() const;
    double t0() const;
    double tmax() const;
    double num_age_groups() const;

    template <typename FP>
    mio::osecirvvs::Model<FP> create_model();

    // template <typename FP>
    // mio::osecirvvs::Model<FP> create_synthetic_model() const;

private:
    template <typename FP>
    mio::IOResult<void> set_covid_parameters(mio::osecirvvs::Parameters<FP>& params);

    template <typename FP>
    mio::IOResult<void> set_contact_matrices(mio::osecirvvs::Parameters<FP>& params);

    template <typename FP>
    mio::IOResult<void> set_population_data(mio::osecirvvs::Model<FP>& model);

    template <typename FP>
    mio::IOResult<void> set_initial_values(mio::osecirvvs::Model<FP>& model);

    fs::path m_data_directory;
    double m_t0;
    double m_tmax;
    double m_num_age_groups;
};

template <typename FP>
mio::IOResult<void> OptimizationModel::set_covid_parameters(mio::osecirvvs::Parameters<FP>& params)
{
    BOOST_OUTCOME_TRY(auto&& parameter_list, mio::read_json((m_data_directory / "data" / "parameter_list.json").string()));

    std::map<std::string, std::string> id_to_name{};
    for(auto entry:  parameter_list)
    {
        id_to_name[entry["id"].asString()] = entry["name"].asString();
    }

    // only uses group 0 for all parameters
    BOOST_OUTCOME_TRY(auto&& scenario_data_run, mio::read_json((m_data_directory / "data" / "scenario_data_run.json").string()));
    std::map<std::string, double> parameter_values{};
    for(auto parameter:  scenario_data_run["modelParameters"])
    {
        std::string parameter_name = id_to_name[parameter["parameterId"].asString()];
        parameter_values[parameter_name] = 0.5 * (parameter["values"][0]["valueMin"].asDouble() + 
                                                    parameter["values"][0]["valueMax"].asDouble());
    }

    //times
    set_parameters<FP, mio::osecirvvs::TimeExposed<FP>>(params, parameter_values["TimeExposed"]);
    set_parameters<FP, mio::osecirvvs::TimeInfectedNoSymptoms<FP>>(params, parameter_values["TimeInfectedNoSymptoms"]);
    set_parameters<FP, mio::osecirvvs::TimeInfectedSymptoms<FP>>(params, parameter_values["TimeInfectedSymptoms"]);
    set_parameters<FP, mio::osecirvvs::TimeInfectedSevere<FP>>(params, parameter_values["TimeInfectedSevere"]);
    set_parameters<FP, mio::osecirvvs::TimeInfectedCritical<FP>>(params, parameter_values["TimeInfectedCritical"]);

    //probabilities
    set_parameters<FP, mio::osecirvvs::TransmissionProbabilityOnContact<FP>>(params, parameter_values["TransmissionProbabilityOnContact"]);
    set_parameters<FP, mio::osecirvvs::RelativeTransmissionNoSymptoms<FP>>(params, parameter_values["RelativeTransmissionNoSymptoms"]);
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    set_parameters<FP, mio::osecirvvs::RiskOfInfectionFromSymptomatic<FP>>(params, parameter_values["RiskOfInfectionFromSymptomatic"]);
    set_parameters<FP, mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<FP>>(params, parameter_values["MaxRiskOfInfectionFromSymptomatic"]);
    set_parameters<FP, mio::osecirvvs::RecoveredPerInfectedNoSymptoms<FP>>(params, parameter_values["RecoveredPerInfectedNoSymptoms"]);
    set_parameters<FP, mio::osecirvvs::SeverePerInfectedSymptoms<FP>>(params, parameter_values["SeverePerInfectedSymptoms"]);
    set_parameters<FP, mio::osecirvvs::CriticalPerSevere<FP>>(params, parameter_values["CriticalPerSevere"]);
    set_parameters<FP, mio::osecirvvs::DeathsPerCritical<FP>>(params, parameter_values["DeathsPerCritical"]);

    set_parameters<FP, mio::osecirvvs::ReducExposedPartialImmunity<FP>>(params, parameter_values["ReducedExposedPartialImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducExposedImprovedImmunity<FP>>(params, parameter_values["ReducedExposedImprovedImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<FP>>(params, parameter_values["ReducedInfectedSymptomsPartialImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<FP>>(params, parameter_values["ReducedInfectedSymptomsImprovedImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<FP>>(params, parameter_values["ReducedInfectedSevereCriticalDeadPartialImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>(params, parameter_values["ReducedInfectedSevereCriticalDeadImprovedImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducTimeInfectedMild<FP>>(params, parameter_values["ReducedTimeInfectedMild"]);

    BOOST_OUTCOME_TRY(mio::Date start_date, mio::parse_date(scenario_data_run["startDate"].asString()));
    params.template get<mio::osecirvvs::StartDay>() = get_day_in_year(start_date);
    params.template get<mio::osecirvvs::Seasonality<FP>>() = parameter_values["Seasonality"];

    return mio::success();
}

template <typename FP>
mio::IOResult<void> OptimizationModel::set_contact_matrices(mio::osecirvvs::Parameters<FP>& params)
{
    auto contact_matrices = mio::ContactMatrixGroup<FP>(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(
                              (m_data_directory / "data" / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));

        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline.cast<FP>();
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic>::Zero(6, 6);
    }
    params.template get<mio::osecirvvs::ContactPatterns<FP>>() = mio::UncertainContactMatrix<FP>(contact_matrices);

    return mio::success();
}

template <typename FP>
mio::IOResult<void> OptimizationModel::set_population_data(mio::osecirvvs::Model<FP>& model)
{  

    const fs::path filename = (m_data_directory / "data"  / "compartment_initialization_2025-05-26" / "initialization_sum.txt");
    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return mio::failure(mio::StatusCode::FileNotFound, (filename).string());
    }

    try {
        std::string tp;
        size_t linenumber = 0;
        while (getline(file, tp)) {
            auto line = mio::split(tp, '\t');

            for (size_t j = 0; j < size_t(line.size()); j++) {
                model.populations[{(mio::AgeGroup)linenumber, (mio::osecirvvs::InfectionState)j}] = std::stod(line[j]);
            }
            linenumber++;
        }
    }
    catch (std::runtime_error& ex) {
        return mio::failure(mio::StatusCode::InvalidFileFormat, (filename).string() + ": " + ex.what());
    }
    // mio::osecirvvs::draw_sample_demographics(model);

    // std::cout << ad::value(model.populations[{mio::AgeGroup(5), mio::osecirvvs::InfectionState::SusceptibleNaive}].value()) << "\n";
    
    return mio::success();
}

/**
 * @brief set_initial_values sets the initial values of the pandemic ODE
 * @param model an instance of the pandemic model
 */
template <typename FP>
mio::IOResult<void> OptimizationModel::set_initial_values(mio::osecirvvs::Model<FP>& model)
{

    mio::osecirvvs::Parameters<FP> params(m_num_age_groups);

    BOOST_OUTCOME_TRY(set_covid_parameters<FP>(params));
    BOOST_OUTCOME_TRY(set_contact_matrices<FP>(params));
    model.parameters = params;

    std::cout << model.parameters.template get<mio::osecirvvs::StartDay>() << "\n";

    // BOOST_OUTCOME_TRY(set_synthetic_population_data(model));
    BOOST_OUTCOME_TRY(set_population_data<FP>(model));
    model.apply_constraints();

    return mio::success();
}

template <typename FP>
mio::osecirvvs::Model<FP> OptimizationModel::create_model()
{
    mio::osecirvvs::Model<FP> model(m_num_age_groups);
    auto out = set_initial_values<FP>(model);
    return model;
}

// -------------------------- //
// --- Artificial Example --- //
// -------------------------- //
// template <typename FP>
// mio::osecirvvs::Model<FP> OptimizationModel::create_synthetic_model() const
// {
//     // --------------------------- //
//     // Define the model parameters //
//     constexpr size_t num_age_groups = m_num_age_groups;
//     mio::osecirvvs::Model<FP> model(num_age_groups);
//     auto& params = model.parameters;

//     params.template get<mio::osecirvvs::ICUCapacity<FP>>()                    = 100;
//     params.template get<mio::osecirvvs::TestAndTraceCapacity<FP>>()           = 0.0143;
//     params.template get<mio::osecirvvs::Seasonality<FP>>()                    = 0.2;
//     params.template get<mio::osecirvvs::DynamicNPIsImplementationDelay<FP>>() = 7;

//     size_t tmax_days = static_cast<size_t>(this->tmax());
//     params.template get<mio::osecirvvs::DailyPartialVaccinations<FP>>().resize(mio::SimulationDay(tmax_days + 1));
//     params.template get<mio::osecirvvs::DailyFullVaccinations<FP>>().resize(mio::SimulationDay(tmax_days + 1));

//     size_t daily_vaccinations = 10;
//     for (size_t age_group = 0; age_group < num_age_groups; age_group++) {
//         for (size_t current_day = 0; current_day <= tmax_days; current_day++) {
//             mio::Index<mio::AgeGroup, mio::SimulationDay> index{mio::AgeGroup(age_group),
//                                                                 mio::SimulationDay(current_day)};
//             FP num_vaccinations = static_cast<FP>(current_day * daily_vaccinations);
//             params.template get<mio::osecirvvs::DailyPartialVaccinations<FP>>()[index] = num_vaccinations;
//             params.template get<mio::osecirvvs::DailyFullVaccinations<FP>>()[index]    = num_vaccinations;
//         }
//     }

//     auto set_all_groups = [&](auto Tag, FP value) {
//         auto& age_group_params = params.template get<decltype(Tag)>();
//         for (size_t age_group = 0; age_group < num_age_groups; age_group++) {
//             age_group_params[mio::AgeGroup(age_group)] = value;
//         }
//     };

//     // — times (all groups same value)
//     set_all_groups(mio::osecirvvs::TimeExposed<FP>{}, 3.33);
//     set_all_groups(mio::osecirvvs::TimeInfectedNoSymptoms<FP>{}, 1.87);
//     set_all_groups(mio::osecirvvs::TimeInfectedSymptoms<FP>{}, 7);
//     set_all_groups(mio::osecirvvs::TimeInfectedSevere<FP>{}, 6);
//     set_all_groups(mio::osecirvvs::TimeInfectedCritical<FP>{}, 7);
//     // — probabilities (all groups same value)
//     set_all_groups(mio::osecirvvs::TransmissionProbabilityOnContact<FP>{}, 0.15);
//     set_all_groups(mio::osecirvvs::RelativeTransmissionNoSymptoms<FP>{}, 0.5);
//     set_all_groups(mio::osecirvvs::RiskOfInfectionFromSymptomatic<FP>{}, 0.0);
//     set_all_groups(mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<FP>{}, 0.4);
//     set_all_groups(mio::osecirvvs::RecoveredPerInfectedNoSymptoms<FP>{}, 0.2);
//     set_all_groups(mio::osecirvvs::SeverePerInfectedSymptoms<FP>{}, 0.1);
//     set_all_groups(mio::osecirvvs::CriticalPerSevere<FP>{}, 0.1);
//     set_all_groups(mio::osecirvvs::DeathsPerCritical<FP>{}, 0.1);
//     // — immunity reductions (all groups same value)
//     set_all_groups(mio::osecirvvs::ReducExposedPartialImmunity<FP>{}, 0.8);
//     set_all_groups(mio::osecirvvs::ReducExposedImprovedImmunity<FP>{}, 0.331);
//     set_all_groups(mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<FP>{}, 0.65);
//     set_all_groups(mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<FP>{}, 0.243);
//     set_all_groups(mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<FP>{}, 0.1);
//     set_all_groups(mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>{}, 0.091);
//     set_all_groups(mio::osecirvvs::ReducTimeInfectedMild<FP>{}, 0.9);

//     // --------------------------- //
//     // Define the contact matrices //
//     int num_contact_locations = 4;
//     Eigen::Matrix<FP, num_age_groups, num_age_groups> baseline_home;
//     Eigen::Matrix<FP, num_age_groups, num_age_groups> baseline_school_pf_eig;
//     Eigen::Matrix<FP, num_age_groups, num_age_groups> baseline_work;
//     Eigen::Matrix<FP, num_age_groups, num_age_groups> baseline_other;

//     baseline_home << 0.4413, 0.4504, 1.2383, 0.8033, 0.0494, 0.0017, 0.0485, 0.7616, 0.6532, 1.1614, 0.0256, 0.0013,
//         0.1800, 0.1795, 0.8806, 0.6413, 0.0429, 0.0032, 0.0495, 0.2639, 0.5189, 0.8277, 0.0679, 0.0014, 0.0087, 0.0394,
//         0.1417, 0.3834, 0.7064, 0.0447, 0.0292, 0.0648, 0.1248, 0.4179, 0.3497, 0.1544;

//     baseline_school_pf_eig << 2.9964, 0.2501, 0.9132, 0.2509, 0.0000, 0.0000, 0.2210, 1.9155, 0.2574, 0.2863, 0.0070,
//         0.0000, 0.0324, 0.3598, 1.2613, 0.1854, 0.0041, 0.0000, 0.1043, 0.4794, 1.1886, 0.1836, 0.0052, 0.0022, 0.0000,
//         0.1150, 0.0000, 0.0000, 0.0168, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000;

//     baseline_other << 0.5170, 0.3997, 0.7957, 0.9958, 0.3239, 0.0428, 0.0632, 0.9121, 0.3254, 0.4731, 0.2355, 0.0148,
//         0.0336, 0.1604, 1.7529, 0.8622, 0.1440, 0.0077, 0.0204, 0.1444, 0.5738, 1.2127, 0.3433, 0.0178, 0.0371, 0.0393,
//         0.4171, 0.9666, 0.7495, 0.0257, 0.0791, 0.0800, 0.3480, 0.5588, 0.2769, 0.0180;

//     baseline_work << 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
//         0.0000, 0.0127, 1.7570, 1.6050, 0.0133, 0.0000, 0.0000, 0.0020, 1.0311, 2.3166, 0.0098, 0.0000, 0.0000, 0.0002,
//         0.0194, 0.0325, 0.0003, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000;

//     Eigen::Matrix<FP, num_age_groups, num_age_groups> minimum_home;
//     Eigen::Matrix<FP, num_age_groups, num_age_groups> minimum_school_pf_eig;
//     Eigen::Matrix<FP, num_age_groups, num_age_groups> minimum_work;
//     Eigen::Matrix<FP, num_age_groups, num_age_groups> minimum_other;

//     minimum_home.setZero();
//     minimum_school_pf_eig.setZero();
//     minimum_work.setZero();
//     minimum_other.setZero();

//     mio::ContactMatrixGroup<FP> contact_matrices(num_contact_locations, num_age_groups);

//     contact_matrices[0].get_baseline() = baseline_home;
//     contact_matrices[1].get_baseline() = baseline_school_pf_eig;
//     contact_matrices[2].get_baseline() = baseline_work;
//     contact_matrices[3].get_baseline() = baseline_other;

//     contact_matrices[0].get_minimum() = minimum_home;
//     contact_matrices[1].get_minimum() = minimum_school_pf_eig;
//     contact_matrices[2].get_minimum() = minimum_work;
//     contact_matrices[3].get_minimum() = minimum_other;

//     params.template get<mio::osecirvvs::ContactPatterns<FP>>() =
//         mio::UncertainContactMatrix<FP>(std::move(contact_matrices));

//     // ----------------------------- //
//     // Define the initial population //
//     size_t num_infection_states = static_cast<size_t>(mio::osecirvvs::InfectionState::Count);
//     Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic> population(num_age_groups, num_infection_states);
//     population << 2753030.497430, 0.000000, 11.708820, 0.000000, 2.212815, 7.578866, 0.000000, 4.084619, 0.000000,
//         0.000000, 0.000000, 22.853384, 0.000000, 2.146616, 0.000000, 0.000000, 0.000000, 0.150378, 0.000000, 0.001865,
//         0.150441, 0.000000, 0.005732, 1051258.609035, 64.000000, 0.000000, 0.000000, 3326384.779416, 0.000000, 5.015423,
//         0.000000, 3.287046, 2.945401, 0.000000, 5.907620, 0.000000, 0.000000, 0.000000, 5.168372, 0.000000, 1.831628,
//         0.000000, 0.000000, 0.000000, 0.023957, 0.000000, 0.001460, 0.138450, 0.000000, 0.017723, 4506996.883505,
//         45.000000, 0.000000, 0.000000, 7474316.734280, 0.000000, 15.185091, 0.000000, 9.963930, 10.119921, 0.000000,
//         19.347768, 0.000000, 0.000000, 0.000000, 27.662806, 0.000000, 9.908623, 0.000000, 0.000000, 0.000000, 0.839570,
//         0.000000, 0.044404, 0.349355, 0.000000, 0.046282, 11247005.297970, 510.000000, 0.000000, 0.000000,
//         13079641.709818, 0.000000, 28.087697, 0.000000, 14.523549, 23.028794, 0.000000, 36.771572, 0.000000, 0.000000,
//         0.000000, 74.950401, 0.000000, 21.335313, 0.000000, 0.000000, 0.000000, 7.139552, 0.000000, 0.287452, 2.329300,
//         0.000000, 0.231930, 15095393.104624, 8936.000000, 0.000000, 0.000000, 12370660.206757, 0.000000, 66.165784,
//         0.000000, 13.006396, 51.098574, 0.000000, 30.752489, 0.000000, 0.000000, 0.000000, 162.557513, 0.000000,
//         16.728201, 0.000000, 0.000000, 0.000000, 38.893189, 0.000000, 0.564814, 13.248394, 0.000000, 0.494789,
//         5239420.683100, 56248.571429, 0.000000, 0.000000, 5614932.848062, 0.000000, 61.140294, 0.000000, 9.149288,
//         46.298860, 0.000000, 21.109925, 0.000000, 0.000000, 0.000000, 137.630764, 0.000000, 11.369236, 0.000000,
//         0.000000, 0.000000, 59.433799, 0.000000, 0.673454, 24.245028, 0.000000, 0.742577, 1695860.958713, 121825.428571,
//         0.000000, 0.000000;

//     for (Eigen::Index column = 0; column < population.cols(); column++) {
//         for (Eigen::Index row = 0; row < population.rows(); row++) {
//             mio::Index<mio::AgeGroup, mio::osecirvvs::InfectionState> index{mio::AgeGroup(row),
//                                                                             mio::osecirvvs::InfectionState(column)};
//             model.populations[index] = population(row, column);
//         }
//     }

//     model.apply_constraints();

//     return model;
// }