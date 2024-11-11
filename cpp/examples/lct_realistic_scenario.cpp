/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Lena Ploetzke
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
#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"
#include "lct_secir/parameters_io.h"

#include "memilio/config.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/io/epi_data.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>
#include <iostream>
#include <vector>

namespace params
{
// Necessary because num_subcompartments is used as a template argument and has to be a constexpr.
constexpr int num_subcompartments = NUM_SUBCOMPARTMENTS;
constexpr size_t num_groups       = 6;

// Parameters
std::map<std::string, ScalarType> simulation_parameter = {{"RelativeTransmissionNoSymptoms", 1},
                                                          {"RiskOfInfectionFromSymptomatic", 0.3},
                                                          {"Seasonality", 0.},
                                                          {"scale_confirmed_cases", 1.},
                                                          {"scale_contacts", 1.},
                                                          {"npi_size", 371 * 14 / (45 * 401.)},
                                                          {"tmax", 45}};
mio::Date start_date(2021, 01, 01);
const ScalarType dt                = 0.01;
const ScalarType age_group_sizes[] = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};
const ScalarType total_population  = 83155031.0;

const ScalarType TransmissionProbabilityOnContact[] = {0.03, 0.06, 0.06, 0.06, 0.09, 0.175};

const ScalarType TimeExposed[]            = {3.335, 3.335, 3.335, 3.335, 3.335, 3.335};
const ScalarType TimeInfectedNoSymptoms[] = {2.74, 2.74, 2.565, 2.565, 2.565, 2.565};
const ScalarType TimeInfectedSymptoms[]   = {7.02625, 7.02625, 7.0665, 6.9385, 6.835, 6.775};
const ScalarType TimeInfectedSevere[]     = {5., 5., 5.925, 7.55, 8.5, 11.};
const ScalarType TimeInfectedCritical[]   = {6.95, 6.95, 6.86, 17.36, 17.1, 11.6};

const ScalarType RecoveredPerInfectedNoSymptoms[] = {1 - 0.75, 1 - 0.75, 1 - 0.8, 1 - 0.8, 1 - 0.8, 1 - 0.8};
const ScalarType SeverePerInfectedSymptoms[]      = {0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225};
const ScalarType CriticalPerSevere[]              = {0.075, 0.075, 0.075, 0.15, 0.3, 0.4};
const ScalarType DeathsPerCritical[]              = {0.05, 0.05, 0.14, 0.14, 0.4, 0.6};

/**
 * @brief Indices of contact matrix corresponding to locations where contacts occur.
 */
enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count
};

// Map the ContactLocation%s to file names.
static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}};
} // namespace params

mio::TimeSeries<ScalarType> add_age_groups(mio::TimeSeries<ScalarType> ageres)
{
    using namespace params;
    size_t infstatecount = (size_t)mio::lsecir::InfectionState::Count;
    mio::TimeSeries<ScalarType> noage(infstatecount);

    for (Eigen::Index timepoint = 0; timepoint < ageres.get_num_time_points(); ++timepoint) {
        Eigen::VectorXd result = Eigen::VectorXd::Zero(infstatecount);
        for (size_t infstate = 0; infstate < infstatecount; infstate++) {
            for (size_t group = 0; group < num_groups; group++) {
                result[infstate] += ageres.get_value(timepoint)[group * infstatecount + infstate];
            }
        }
        noage.add_time_point(ageres.get_time(timepoint), result);
    }

    return noage;
}

/**
 * @brief Add NPIs to a given contact matrix from 01/10/2020 on.
 *
 * NPIs from the Paper "Assessment of effective mitigation ..." (doi: 10.1016/j.mbs.2021.108648) are used with slight 
 * modifications for a period of 45 days from 01/10/2020 on.
 * 
 * @param[in,out] contact_matrices The contact matrices where the NPIs should be added to.
 * @param[in] start_date Start date of the simulation used for setting the NPIs.
 */
void set_npi_october(mio::ContactMatrixGroup& contact_matrices)
{
    using namespace params;
    // ---------------------24/10/2020--------------------------------
    auto offset_npi = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 10, 25), start_date));
    ScalarType contact_reduction_all_locations = simulation_parameter["npi_size"];
    for (auto&& contact_location : contact_locations) {
        contact_matrices[size_t(contact_location.first)].add_damping(
            Eigen::MatrixXd::Constant(num_groups, num_groups, contact_reduction_all_locations), offset_npi);
    }
}

void set_npi_july(mio::ContactMatrixGroup& contact_matrices)
{
    using namespace params;
    // -----------------------------------------------------
    auto offset_npi = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 07, 07), start_date));
    ScalarType contact_reduction_all_locations = simulation_parameter["npi_size"];
    for (auto&& contact_location : contact_locations) {
        contact_matrices[size_t(contact_location.first)].add_damping(
            Eigen::MatrixXd::Constant(num_groups, num_groups, contact_reduction_all_locations), offset_npi);
    }
}

/**
 * @brief Set the contact pattern of parameters for a Model without division in age groups.
 *
 * The contacts are calculated using contact matrices from files in the data directory for different locations.
 * Also set Non-pharmaceutical Interventions influencing the ContactPatterns used for simulation in the period from start_date to end_date.
 * 
 * @param[in] data_dir Directory to files with minimum and baseline contact matrices.
 * @returns Any io errors that happen during reading of the input files.
 */
mio::IOResult<mio::UncertainContactMatrix<ScalarType>> get_contact_matrix(const fs::path& contact_data_dir)
{
    using namespace params;
    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), num_groups);

    // Load and set baseline contacts for each contact location.
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(
            auto&& baseline,
            mio::read_mobility_plain(
                (contact_data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));
        contact_matrices[size_t(contact_location.first)].get_baseline() =
            simulation_parameter["scale_contacts"] * baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum() = Eigen::MatrixXd::Zero(num_groups, num_groups);
    }

    // ----- Add NPIs to the contact matrices. -----
    mio::Date end_date  = mio::offset_date_by_days(start_date, simulation_parameter["tmax"]);
    auto start_npi_july = mio::Date(2020, 7, 7);
    if ((start_npi_july < end_date) & (start_date < start_npi_july)) {
        set_npi_july(contact_matrices);
    }
    // Set of NPIs for October.
    auto start_npi_october = mio::Date(2020, 10, 25);
    if ((start_npi_october < end_date) & (start_date < start_npi_october)) {
        set_npi_october(contact_matrices);
    }
    return mio::success(mio::UncertainContactMatrix<ScalarType>(contact_matrices));

    ;
}

/**
 * @brief Performs a simulation of a real scenario with an LCT and an ODE model.
 *
 * @param[in] path Path of the RKI file that should be used to compute initial values for simulations.
 * @param[in] save_dir Specifies the directory where the results should be stored. Provide an empty string if results should not be saved.
 * @returns Any io errors that happen during reading of the RKI file or files for contact matrices or saving the results.
 */
mio::IOResult<void> simulate_other_subcompartments(std::string const& dir_to_contact_data,
                                                   std::string const& infection_data_dir, std::string save_dir = "")
{
    using namespace params;
    std::cout << "Realistic scenario with number of subcompartments is approximately the mean stay time." << std::endl;
    // ----- Initialize age resolved model. -----
    using InfState      = mio::lsecir::InfectionState;
    using LctState0_14  = mio::LctInfectionState<InfState, 1, 3, 3, 7, 5, 7, 1, 1>;
    using LctState15_34 = mio::LctInfectionState<InfState, 1, 3, 3, 7, 6, 7, 1, 1>;
    using LctState35_59 = mio::LctInfectionState<InfState, 1, 3, 3, 7, 8, 17, 1, 1>;
    using LctState60_79 = mio::LctInfectionState<InfState, 1, 3, 3, 7, 9, 17, 1, 1>;
    using LctState80    = mio::LctInfectionState<InfState, 1, 3, 3, 7, 11, 12, 1, 1>;
    using Model =
        mio::lsecir::Model<LctState0_14, LctState0_14, LctState15_34, LctState35_59, LctState60_79, LctState80>;
    Model model;

    // Define parameters used for simulation and initialization.
    for (size_t group = 0; group < num_groups; group++) {
        model.parameters.get<mio::lsecir::TimeExposed>()[group]            = TimeExposed[group];
        model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[group] = TimeInfectedNoSymptoms[group];
        model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[group]   = TimeInfectedSymptoms[group];
        model.parameters.get<mio::lsecir::TimeInfectedSevere>()[group]     = TimeInfectedSevere[group];
        model.parameters.get<mio::lsecir::TimeInfectedCritical>()[group]   = TimeInfectedCritical[group];
        model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[group] =
            TransmissionProbabilityOnContact[group];

        model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[group] =
            simulation_parameter["RelativeTransmissionNoSymptoms"];
        model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[group] =
            simulation_parameter["RiskOfInfectionFromSymptomatic"];

        model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[group] =
            RecoveredPerInfectedNoSymptoms[group];
        model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[group] = SeverePerInfectedSymptoms[group];
        model.parameters.get<mio::lsecir::CriticalPerSevere>()[group]         = CriticalPerSevere[group];
        model.parameters.get<mio::lsecir::DeathsPerCritical>()[group]         = DeathsPerCritical[group];
    }

    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix(dir_to_contact_data));

    model.parameters.get<mio::lsecir::ContactPatterns>() = contact_matrix;
    model.parameters.get<mio::lsecir::Seasonality>()     = simulation_parameter["Seasonality"];
    model.parameters.get<mio::lsecir::StartDay>()        = mio::get_day_in_year(start_date);

    BOOST_OUTCOME_TRY(auto&& rki_data, mio::read_confirmed_cases_data(infection_data_dir));
    auto init = mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesDataEntry>(
        rki_data, model.populations, model.parameters, start_date,
        std::vector<ScalarType>(age_group_sizes, age_group_sizes + num_groups),
        std::vector<ScalarType>(num_groups, simulation_parameter["scale_confirmed_cases"]));
    if (!init) {
        printf("%s\n", init.error().formatted_message().c_str());
        return init;
    }

    // Perform simulation.
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max so that we have a fixed time step and can compare to the result with one group.
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);
    mio::TimeSeries<ScalarType> result =
        mio::simulate<ScalarType, Model>(0, simulation_parameter["tmax"], dt, model, integrator);

    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations                 = model.calculate_compartments(result);
    mio::TimeSeries<ScalarType> populations_accumulated_age = add_age_groups(populations);

    if (!save_dir.empty()) {
        std::string filename = save_dir + "real_" + std::to_string(start_date.year) + "-" +
                               std::to_string(start_date.month) + "-" + std::to_string(start_date.day) + "_var";
        // Age-resolved.
        mio::IOResult<void> save_result_status = mio::save_result({populations}, {0}, 1, filename + "_ageres.h5");
        // Not resolved by age.
        save_result_status = mio::save_result({populations_accumulated_age}, {0}, 1, filename + "_accumulated.h5");
    }
    // Print commands to get the number of new infections on the first day of simulation. Could be used to scale the contacts.
    std::cout << "Number of new infections on the first day of simulation: " << std::endl;
    std::cout << std::fixed << std::setprecision(1)
              << (populations_accumulated_age[0][0] - populations_accumulated_age[1][0]) /
                     (populations.get_time(1) - populations.get_time(0))
              << std::endl;

    return mio::success();
}

/**
 * @brief Performs a simulation of a real scenario with an LCT and an ODE model.
 *
 * @param[in] path Path of the RKI file that should be used to compute initial values for simulations.
 * @param[in] save_dir Specifies the directory where the results should be stored. Provide an empty string if results should not be saved.
 * @returns Any io errors that happen during reading of the RKI file or files for contact matrices or saving the results.
 */
mio::IOResult<void> simulate(std::string const& dir_to_contact_data, std::string const& infection_data_dir,
                             std::string save_dir = "")
{
    using namespace params;
    std::cout << "Realistic scenario with " << num_subcompartments << " subcompartments." << std::endl;
    // ----- Initialize age resolved model. -----
    using InfState = mio::lsecir::InfectionState;
    using LctState = mio::LctInfectionState<InfState, 1, num_subcompartments, num_subcompartments, num_subcompartments,
                                            num_subcompartments, num_subcompartments, 1, 1>;
    using Model    = mio::lsecir::Model<LctState, LctState, LctState, LctState, LctState, LctState>;
    Model model;

    // Define parameters used for simulation and initialization.
    for (size_t group = 0; group < num_groups; group++) {
        model.parameters.get<mio::lsecir::TimeExposed>()[group]            = TimeExposed[group];
        model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[group] = TimeInfectedNoSymptoms[group];
        model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[group]   = TimeInfectedSymptoms[group];
        model.parameters.get<mio::lsecir::TimeInfectedSevere>()[group]     = TimeInfectedSevere[group];
        model.parameters.get<mio::lsecir::TimeInfectedCritical>()[group]   = TimeInfectedCritical[group];
        model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[group] =
            TransmissionProbabilityOnContact[group];

        model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[group] =
            simulation_parameter["RelativeTransmissionNoSymptoms"];
        model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[group] =
            simulation_parameter["RiskOfInfectionFromSymptomatic"];

        model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[group] =
            RecoveredPerInfectedNoSymptoms[group];
        model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[group] = SeverePerInfectedSymptoms[group];
        model.parameters.get<mio::lsecir::CriticalPerSevere>()[group]         = CriticalPerSevere[group];
        model.parameters.get<mio::lsecir::DeathsPerCritical>()[group]         = DeathsPerCritical[group];
    }

    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix(dir_to_contact_data));

    model.parameters.get<mio::lsecir::ContactPatterns>() = contact_matrix;
    model.parameters.get<mio::lsecir::Seasonality>()     = simulation_parameter["Seasonality"];
    model.parameters.get<mio::lsecir::StartDay>()        = mio::get_day_in_year(start_date);

    BOOST_OUTCOME_TRY(auto&& rki_data, mio::read_confirmed_cases_data(infection_data_dir));
    auto init = mio::lsecir::set_initial_data_from_confirmed_cases<Model::Populations, mio::ConfirmedCasesDataEntry>(
        rki_data, model.populations, model.parameters, start_date,
        std::vector<ScalarType>(age_group_sizes, age_group_sizes + num_groups),
        std::vector<ScalarType>(num_groups, simulation_parameter["scale_confirmed_cases"]));
    if (!init) {
        printf("%s\n", init.error().formatted_message().c_str());
        return init;
    }

    // Perform simulation.
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max so that we have a fixed time step and can compare to the result with one group.
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);
    mio::TimeSeries<ScalarType> result =
        mio::simulate<ScalarType, Model>(0, simulation_parameter["tmax"], dt, model, integrator);

    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations                 = model.calculate_compartments(result);
    mio::TimeSeries<ScalarType> populations_accumulated_age = add_age_groups(populations);

    if (!save_dir.empty()) {
        std::string filename = save_dir + "real_" + std::to_string(start_date.year) + "-" +
                               std::to_string(start_date.month) + "-" + std::to_string(start_date.day) + "_" +
                               std::to_string(num_subcompartments);
        mio::IOResult<void> save_result_status = mio::save_result({populations}, {0}, 1, filename + "_ageres.h5");
        if (!save_result_status) {
            return save_result_status;
        }
        save_result_status = mio::save_result({populations_accumulated_age}, {0}, 1, filename + "_accumulated.h5");
        if (!save_result_status) {
            return save_result_status;
        }
    }
    // Print commands to get the number of new infections on the first day of simulation. Could be used to scale the contacts.
    std::cout << "Number of new infections on the first day of simulation: " << std::endl;
    std::cout << std::fixed << std::setprecision(1)
              << (populations_accumulated_age[0][0] - populations_accumulated_age[1][0]) /
                     (populations.get_time(1) - populations.get_time(0))
              << std::endl;

    return mio::success();
}

int main(int argc, char** argv)
{
    std::string dir_to_data    = "../../data";
    bool other_subcompartments = false;
    if (argc > 1) {
        other_subcompartments = std::stoi(argv[1]);
    }
    if (argc > 5) {
        params::start_date = mio::Date(std::stod(argv[2]), std::stod(argv[3]), std::stod(argv[4]));
        dir_to_data        = argv[5];
    }

    if (argc > 9) {
        params::simulation_parameter["RelativeTransmissionNoSymptoms"] = std::stod(argv[6]);
        params::simulation_parameter["RiskOfInfectionFromSymptomatic"] = std::stod(argv[7]);
        params::simulation_parameter["scale_contacts"]                 = std::stod(argv[8]);
        params::simulation_parameter["npi_size"]                       = std::stod(argv[9]);
    }
    std::string save_dir           = dir_to_data + "/simulation_lct_real/";
    std::string infection_data_dir = dir_to_data + "/pydata/Germany/cases_all_age_ma7.json";
    if (other_subcompartments) {
        auto result = simulate_other_subcompartments(dir_to_data, infection_data_dir, save_dir);
        if (!result) {
            printf("%s\n", result.error().formatted_message().c_str());
            return -1;
        }
    }
    else {
        auto result = simulate(dir_to_data, infection_data_dir, save_dir);
        if (!result) {
            printf("%s\n", result.error().formatted_message().c_str());
            return -1;
        }
    }

    return 0;
}
