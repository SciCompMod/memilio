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
constexpr int num_subcompartments = 10;
constexpr size_t num_groups       = 6;

// Parameters
std::map<std::string, ScalarType> simulation_parameter = {{"RelativeTransmissionNoSymptoms", 1},
                                                          {"RiskOfInfectionFromSymptomatic", 0.3},
                                                          {"Seasonality", 0.},
                                                          {"scale_confirmed_cases", 1.},
                                                          {"scale_contacts", 1.},
                                                          {"lockdown_hard", 371 * 14 / (45 * 401.)},
                                                          {"tmax", 45}};
mio::Date start_date(2020, 9, 01);
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

/**
 * @brief Different types of NPI, used as DampingType.
 */
enum class Intervention
{
    Home,
    SchoolClosure,
    HomeOffice,
    GatheringBanFacilitiesClosure,
    PhysicalDistanceAndMasks,
    Count
};

/**
 * @brief Different level of NPI, used as DampingLevel.
 */
enum class InterventionLevel
{
    Main,
    PhysicalDistanceAndMasks,
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
 * @param[in] lockdown_hard Proportion of counties for which a hard lockdown is implemented.
 */
void set_npi_october(mio::ContactMatrixGroup& contact_matrices)
{
    using namespace params;
    // ---------------------01/10/2020--------------------------------
    // NPIs from Paper for october.
    auto offset_npi = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 10, 1), start_date));
    // For the beginning of the time period, we assume only half of the defined proportion of counties is in a hard lockdown.
    ScalarType lockdown_hard = simulation_parameter["lockdown_hard"];
    lockdown_hard            = lockdown_hard / 2;
    // Contact reduction at home.
    ScalarType v = 0.3 * (1 - lockdown_hard) + lockdown_hard * 0.5;
    contact_matrices[size_t(ContactLocation::Home)].add_damping(Eigen::MatrixXd::Constant(1, 1, v),
                                                                mio::DampingLevel(int(InterventionLevel::Main)),
                                                                mio::DampingType(int(Intervention::Home)), offset_npi);
    // Home-Office + people stopped working.
    v = (0.25 + 0.025) * (1 - lockdown_hard) + lockdown_hard * (0.25 + 0.15);
    contact_matrices[size_t(ContactLocation::Work)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::Main)),
        mio::DampingType(int(Intervention::HomeOffice)), offset_npi);
    // GatheringBanFacilitiesClosure affects ContactLocation Other.
    v = 0.1 * (1 - lockdown_hard) + lockdown_hard * 0.7;
    contact_matrices[size_t(ContactLocation::Other)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::Main)),
        mio::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), offset_npi);
    // PhysicalDistanceAndMasks in all locations.
    v = 0.3 * (1 - lockdown_hard) + lockdown_hard * 0.7;
    for (auto&& contact_location : contact_locations) {
        contact_matrices[size_t(contact_location.first)].add_damping(
            Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), offset_npi);
    }
    // Remote schooling.
    v = lockdown_hard * 0.25;
    contact_matrices[size_t(ContactLocation::School)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::Main)),
        mio::DampingType(int(Intervention::SchoolClosure)), offset_npi);

    // ---------------------24/10/2020--------------------------------
    /* We assume that the stricter NPIs of november defined in the paper are beginning about a week earlier, 
    which can be seen from the RKI data.
    Moreover the lockdown value of PhysicalDistanceAndMasks in the location school is assumed to apply for all counties.*/
    offset_npi = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 10, 24), start_date));
    // For the second half of the simulation, the proportion of counties in hard lockdown is increased to compensate for the lower proportion before.
    lockdown_hard = lockdown_hard * 3;
    // Contact reduction at home.
    v = 0.5;
    contact_matrices[size_t(ContactLocation::Home)].add_damping(Eigen::MatrixXd::Constant(1, 1, v),
                                                                mio::DampingLevel(int(InterventionLevel::Main)),
                                                                mio::DampingType(int(Intervention::Home)), offset_npi);
    // Home-Office + people stopped working.
    v = (0.25 + 0.05) * (1 - lockdown_hard) + lockdown_hard * (0.25 + 0.15);
    contact_matrices[size_t(ContactLocation::Work)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::Main)),
        mio::DampingType(int(Intervention::HomeOffice)), offset_npi);
    // GatheringBanFacilitiesClosure affects ContactLocation Other.
    v = 0.7;
    contact_matrices[size_t(ContactLocation::Other)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::Main)),
        mio::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), offset_npi);
    // PhysicalDistanceAndMasks in ContactLocation%s Home.
    v = 0.3 * (1 - lockdown_hard) + lockdown_hard * 0.7;
    contact_matrices[size_t(ContactLocation::Home)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
        mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), offset_npi);
    // PhysicalDistanceAndMasks in ContactLocation%s School.
    v = 0.7;
    contact_matrices[size_t(ContactLocation::School)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
        mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), offset_npi);
    // PhysicalDistanceAndMasks in ContactLocation%s Work and Other.
    v = 0.5 * (1 - lockdown_hard) + lockdown_hard * 0.7;
    contact_matrices[size_t(ContactLocation::Work)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
        mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), offset_npi);
    contact_matrices[size_t(ContactLocation::Other)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
        mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), offset_npi);
    // Remote schooling.
    v = lockdown_hard * 0.25;
    contact_matrices[size_t(ContactLocation::School)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::Main)),
        mio::DampingType(int(Intervention::SchoolClosure)), offset_npi);
}

/**
 * @brief Set the contact pattern of parameters for a Model without division in age groups.
 *
 * The contacts are calculated using contact matrices from files in the data directory for different locations.
 * Also set Non-pharmaceutical Interventions influencing the ContactPatterns used for simulation in the period from start_date to end_date.
 * 
 * @param[in] data_dir Directory to files with minimum and baseline contact matrices.
 * @param[in] parameters Object that the contact pattern will be added to.
 * @param[in] simulation_parameters Map with parameters necessary for the calculation of contacts an NPIs which can be different for different start dates.
 *      Function uses the values to define start and end date, lockdown_hard and scale_contacts.
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
        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::MatrixXd::Zero(num_groups, num_groups);
    }

    // ----- Add NPIs to the contact matrices. -----
    mio::Date end_date = mio::offset_date_by_days(start_date, simulation_parameter["tmax"]);
    // Set of NPIs for October.
    auto start_npi_october = mio::Date(2020, 10, 1);
    if (start_npi_october < end_date) {
        set_npi_october(contact_matrices);
    }
    return mio::success(mio::UncertainContactMatrix<ScalarType>(contact_matrices));

    ;
}

/**
 * @brief Performs a simulation of a real scenario with an LCT and an ODE model.
 *
 * @param[in] path Path of the RKI file that should be used to compute initial values for simulations.
 * @param[in] simulation_parameters Map with parameters necessary for the simulation which can be different for different start dates.
 *        Provide the parameters "start_month", "start_day","seasonality" (parameter k for the seasonality of the models), 
 *        "RelativeTransmissionNoSymptoms", "RiskOfInfectionFromSymptomatic", "scale_confirmed_cases" (to scale the RKI data while computing an initialization vector),
 *        "lockdown_hard" (Proportion of counties for which a hard lockdown is implemented) and 
 *        "scale_contacts" (scales contacts per hand to match the new infections in the RKI data).
 *        The assumption regarding the number of subcompartments of the LCT model can be controlled via the parameter "num_subcompartments".
 * @param[in] save_dir Specifies the directory where the results should be stored. Provide an empty string if results should not be saved.
 * @param[in] print_result Specifies if the results should be printed.
 * @returns Any io errors that happen during reading of the RKI file or files for contact matrices or saving the results.
 */
mio::IOResult<void> simulate(std::string const& infection_data_dir, std::string save_dir = "")
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

    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix("../../data"));

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

int main()
{
    std::string save_dir           = "../../data/simulation_lct_real/";
    std::string infection_data_dir = "../../data/pydata/Germany/cases_all_age_ma7.json";

    /* Values for "RelativeTransmissionNoSymptoms", "RelativeTransmissionNoSymptoms" and "seasonality" are suitable values based on doi: 10.1016/j.mbs.2021.108648.
    "scale_confirmed_cases" are values directly from this paper. 
    "lockdown_hard" is based on the number of beginnings of a strict lockdown for 14 days of the 45 days simulation period in the lockdown. Value is not used for 01/06/2020 and 
    scaled in october to assume that in the beginning of the period are less counties in lockdown and more in the second half. 
    "lockdown_hard" should give the average percentage of counties that are in a hard lockdown on a simulation day. 
    "scale_contacts" is used to match the predicted number of new infections on the first simulation day to the RKI data. */

    // Paths are valid if file is executed eg in memilio/build/bin
    // Simulation with start date 01.06.2020 with 1 subcompartment.
    auto result = simulate(infection_data_dir, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}
