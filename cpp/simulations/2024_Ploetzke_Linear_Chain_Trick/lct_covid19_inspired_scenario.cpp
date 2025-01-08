/* 
* Copyright (C) 2020-2025 MEmilio
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
#include <type_traits>

namespace params
{
// num_subcompartments is used as a template argument and has to be a constexpr.
// If NUM_SUBCOMPARTMENTS is set to zero, an LCT model is used with numbers of subcompartments such that
// each corresponds to the approximate stay time in the compartment.
constexpr int num_subcompartments = NUM_SUBCOMPARTMENTS;
constexpr size_t num_groups       = 6;

// Define (age-resolved) parameters.
mio::Date start_date(2021, 01, 01);
const ScalarType tmax                     = 45;
const ScalarType dt                       = 0.01;
const ScalarType seasonality              = 0.;
ScalarType relativeTransmissionNoSymptoms = 1.;
ScalarType riskOfInfectionFromSymptomatic = 0.3;
const ScalarType age_group_sizes[]        = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};

const ScalarType transmissionProbabilityOnContact[] = {0.03, 0.06, 0.06, 0.06, 0.09, 0.175};

const ScalarType timeExposed[]            = {3.335, 3.335, 3.335, 3.335, 3.335, 3.335};
const ScalarType timeInfectedNoSymptoms[] = {2.74, 2.74, 2.565, 2.565, 2.565, 2.565};
const ScalarType timeInfectedSymptoms[]   = {7.02625, 7.02625, 7.0665, 6.9385, 6.835, 6.775};
const ScalarType timeInfectedSevere[]     = {5., 5., 5.925, 7.55, 8.5, 11.};
const ScalarType timeInfectedCritical[]   = {6.95, 6.95, 6.86, 17.36, 17.1, 11.6};

const ScalarType recoveredPerInfectedNoSymptoms[] = {1 - 0.75, 1 - 0.75, 1 - 0.8, 1 - 0.8, 1 - 0.8, 1 - 0.8};
const ScalarType severePerInfectedSymptoms[]      = {0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225};
const ScalarType criticalPerSevere[]              = {0.075, 0.075, 0.075, 0.15, 0.3, 0.4};
const ScalarType deathsPerCritical[]              = {0.05, 0.05, 0.14, 0.14, 0.4, 0.6};

// Define scalings that can be changed via command line.
ScalarType scale_confirmed_cases = 1.; ///< Scale confirmed case data to incorporate a detection ratio.
// Scale the contact data so that the simulation results for the number of daily new transmissions align with the
// extrapolated RKI data (for aggregated age groups).
ScalarType scale_contacts = 1.;
ScalarType npi_size       = 0.; ///< Size of the NPI to be implemented from 25/20/2020.

// Define contact locations, which are used to define the contact matrix and to define location-resolved NPIs.
/// @brief Indices of contact matrix corresponding to locations where contacts occur.
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

/** 
* @brief Function to transform an age-resolved simulation result into a result without age resolution.
*
* Sums up the values in the age groups to transform the simulation result into a result without age resolution. 
* To provide a clear overview, we use non-age-resolved results for visualizations.
* This implementation is only valid if the simulation is run with equal LctStates for all groups or if the result 
* does not contain any subcompartments.
*   
* @param[in] ageresolved_result TimeSeries with an age-resolved simulation result.
* @returns TimeSeries with the result where the values of the age groups are summed up.
*/
mio::TimeSeries<ScalarType> sum_age_groups(const mio::TimeSeries<ScalarType> ageresolved_result)
{
    using namespace params;
    size_t infstatecount = size_t((ScalarType)ageresolved_result.get_num_elements() / (ScalarType)num_groups);
    mio::TimeSeries<ScalarType> nonageresolved_result(infstatecount);

    // For each time point, calculate the result without age resolution and add the time point
    // to the non-age-resolved result.
    for (Eigen::Index timepoint = 0; timepoint < ageresolved_result.get_num_time_points(); ++timepoint) {
        Eigen::VectorX<ScalarType> result = Eigen::VectorX<ScalarType>::Zero(infstatecount);
        for (size_t infstate = 0; infstate < infstatecount; infstate++) {
            for (size_t group = 0; group < num_groups; group++) {
                result[infstate] += ageresolved_result.get_value(timepoint)[group * infstatecount + infstate];
            }
        }
        nonageresolved_result.add_time_point(ageresolved_result.get_time(timepoint), result);
    }

    return nonageresolved_result;
}

/**
 * @brief Add NPIs to a given contact matrix from 25/10/2020 onwards.
 *
 * Add NPIs to the contact matrix contact_matrices from 25/10/2020 onwards. 
 * The size of the NPI is specified in the variable npi_size. The NPI is applied to all contact locations equally.
 * 
 * @param[in,out] contact_matrices The contact matrices where the NPIs should be added to.
 */
void set_npi_october(mio::ContactMatrixGroup& contact_matrices)
{
    using namespace params;
    auto offset_npi = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 10, 25), start_date));
    for (auto&& contact_location : contact_locations) {
        contact_matrices[size_t(contact_location.first)].add_damping(
            Eigen::MatrixXd::Constant(num_groups, num_groups, npi_size), offset_npi);
    }
}

/**
 * @brief Set the contact pattern from data files.
 *
 * The contacts are set using contact matrices from files in the data directory for different locations.
 * Non-pharmaceutical interventions (NPIs) influencing the ContactPatterns are set.
 * The contact matrices are scaled using the parameter scale_contacts such that the simulation results for the number
 *  of daily new transmissions align with the extrapolated RKI data.
 * 
 * @param[in] contact_data_dir Directory to files with minimum and baseline contact matrices.
 * @returns Any io errors that happen during reading of the input files.
 */
mio::IOResult<mio::UncertainContactMatrix<ScalarType>> get_contact_matrix(std::string contact_data_dir)
{
    using namespace params;
    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), num_groups);

    // Load and set baseline contacts for each contact location.
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(contact_data_dir + "baseline_" + contact_location.second + ".txt"));
        contact_matrices[size_t(contact_location.first)].get_baseline() = scale_contacts * baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::MatrixXd::Zero(num_groups, num_groups);
    }

    // Add NPIs to the contact matrices.
    mio::Date end_date     = mio::offset_date_by_days(start_date, (int)tmax);
    auto start_npi_october = mio::Date(2020, 10, 25);
    if ((int)(start_npi_october < end_date) & (int)(start_date < start_npi_october)) {
        set_npi_october(contact_matrices);
    }
    return mio::success(mio::UncertainContactMatrix<ScalarType>(contact_matrices));
}

/**
 * @brief Perform simulation for a Covid-19 inspired scenario in Germany.
 *
 *   The simulation uses LCT models with Covid-19 inspired parameters and a contact rate for Germany. The initial values
 *   are set using reported data provided by the RKI (and the DIVI data set).
 *   The scaling factors scale_confirmed_cases, scale_contacts and npi_size are used for their defined reasons.
 *   The simulation results are stored with age resolution in a file with the suffix "_ageresolved" and without 
 *   age resolution (age-resolved results are accumulated) in a file with suffix "_accumulated".
 *   
 *   The LCT model is constructed with NUM_SUBCOMPARTMENT subcompartments for all compartments where subcompartments 
 *   make sense (so all except Susceptibles, Recovered and Dead) for all age groups.
 *   If NUM_SUBCOMPARTMENTS is set to zero, an LCT model is used with numbers of subcompartments such that 
 *   each corresponds to the approximate stay time in the compartment.
 *
 * @param[in] contact_data_dir Directory to the contact data.
 * @param[in] infection_data_dir Directory to infection data provided by the RKI for Germany.
 * @param[in] divi_data_dir Directory to DIVI data regarding the number of patients in intensive care units.
 * @param[in] save_dir Specifies the directory where the results should be stored. 
 *   Provide an empty string if the results should not be saved.
 * @returns Any io errors that happen during reading of the RKI file or files for contact matrices or saving the results.
 */
mio::IOResult<void> simulate(std::string const& contact_data_dir, std::string const& infection_data_dir,
                             std::string const& divi_data_dir, std::string save_dir = "")
{
    using namespace params;
    std::cout << "Realistic scenario with " << num_subcompartments << " subcompartments." << std::endl;
    // Initialize (age-resolved) model.
    using InfState = mio::lsecir::InfectionState;
    // Define appropriate LCT model type: 1.) An LCT model with NUM_SUBCOMPARTMENTS subcompartments for all compartments
    // and age groups if NUM_SUBCOMPOARTMENTS if greater than zero or 2.) if the value is zero, an LCT model
    // with numbers of subcompartments such that each corresponds to the approximate stay time in the compartment.
    // For 1.) : Define single LctState.
    // If NUM_SUBCOMPOARTMENTS=0, define LctState as a random other LctInfectionState as the template
    // arguments for LctInfectionState have to be greater than zero.
    using LctState = std::conditional<
        (num_subcompartments == 0), mio::LctInfectionState<InfState, 1, 1, 1, 1, 1, 1, 1, 1>,
        mio::LctInfectionState<InfState, 1, num_subcompartments, num_subcompartments, num_subcompartments,
                               num_subcompartments, num_subcompartments, 1, 1>>::type;
    // Define LctStates for 2.): We need to define a separate LctState for each age group.
    using LctState0_14  = mio::LctInfectionState<InfState, 1, 3, 3, 7, 5, 7, 1, 1>;
    using LctState15_34 = mio::LctInfectionState<InfState, 1, 3, 3, 7, 6, 7, 1, 1>;
    using LctState35_59 = mio::LctInfectionState<InfState, 1, 3, 3, 7, 8, 17, 1, 1>;
    using LctState60_79 = mio::LctInfectionState<InfState, 1, 3, 3, 7, 9, 17, 1, 1>;
    using LctState80    = mio::LctInfectionState<InfState, 1, 3, 3, 7, 11, 12, 1, 1>;
    // Decide which LCT model type should be used.
    using Model = std::conditional<
        (num_subcompartments == 0),
        mio::lsecir::Model<LctState0_14, LctState0_14, LctState15_34, LctState35_59, LctState60_79, LctState80>,
        mio::lsecir::Model<LctState, LctState, LctState, LctState, LctState, LctState>>::type;
    Model model;

    // Define parameters used for simulation and initialization.
    for (size_t group = 0; group < num_groups; group++) {
        model.parameters.get<mio::lsecir::TimeExposed>()[group]            = timeExposed[group];
        model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()[group] = timeInfectedNoSymptoms[group];
        model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()[group]   = timeInfectedSymptoms[group];
        model.parameters.get<mio::lsecir::TimeInfectedSevere>()[group]     = timeInfectedSevere[group];
        model.parameters.get<mio::lsecir::TimeInfectedCritical>()[group]   = timeInfectedCritical[group];
        model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>()[group] =
            transmissionProbabilityOnContact[group];

        model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>()[group] = relativeTransmissionNoSymptoms;
        model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>()[group] = riskOfInfectionFromSymptomatic;

        model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>()[group] =
            recoveredPerInfectedNoSymptoms[group];
        model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()[group] = severePerInfectedSymptoms[group];
        model.parameters.get<mio::lsecir::CriticalPerSevere>()[group]         = criticalPerSevere[group];
        model.parameters.get<mio::lsecir::DeathsPerCritical>()[group]         = deathsPerCritical[group];
    }

    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix(contact_data_dir));

    model.parameters.get<mio::lsecir::ContactPatterns>() = contact_matrix;
    model.parameters.get<mio::lsecir::Seasonality>()     = seasonality;
    model.parameters.get<mio::lsecir::StartDay>()        = mio::get_day_in_year(start_date);

    // Set initial values using reported data.
    BOOST_OUTCOME_TRY(auto&& rki_data, mio::read_confirmed_cases_data(infection_data_dir));
    BOOST_OUTCOME_TRY(auto&& divi_data, mio::read_divi_data(divi_data_dir));
    auto init = mio::lsecir::set_initial_values_from_reported_data<Model::Populations, mio::ConfirmedCasesDataEntry>(
        rki_data, model.populations, model.parameters, start_date,
        std::vector<ScalarType>(age_group_sizes, age_group_sizes + num_groups),
        std::vector<ScalarType>(num_groups, scale_confirmed_cases), divi_data);
    if (!init) {
        printf("%s\n", init.error().formatted_message().c_str());
        return init;
    }

    // Set integrator of fifth order with fixed step size and perform simulation.
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max to get a fixed step size.
    integrator->set_dt_min(dt);
    integrator->set_dt_max(dt);
    mio::TimeSeries<ScalarType> result = mio::simulate<ScalarType, Model>(0, tmax, dt, model, integrator);

    // Calculate result without division in subcompartments and without division in age groups.
    mio::TimeSeries<ScalarType> populations                 = model.calculate_compartments(result);
    mio::TimeSeries<ScalarType> populations_accumulated_age = sum_age_groups(populations);

    if (!save_dir.empty()) {
        std::string filename = save_dir + "lct_" + std::to_string(start_date.year) + "-" +
                               std::to_string(start_date.month) + "-" + std::to_string(start_date.day) + "_subcomp" +
                               std::to_string(num_subcompartments);
        mio::IOResult<void> save_result_status = mio::save_result({populations}, {0}, 1, filename + "_ageresolved.h5");
        if (!save_result_status) {
            return save_result_status;
        }
        save_result_status = mio::save_result({populations_accumulated_age}, {0}, 1, filename + "_accumulated.h5");
        if (!save_result_status) {
            return save_result_status;
        }
    }
    // Print the predicted number of daily new transmissions at the start of the simulation.
    // Could be used to compare the result with the extrapolated reported data and, afterwards, to scale the
    // contacts using the variable scale_contacts.
    std::cout << "Predicted number of daily new transmissions at the start of the simulation.: " << std::endl;
    std::cout << std::fixed << std::setprecision(1)
              << (populations_accumulated_age[0][0] - populations_accumulated_age[1][0]) /
                     (populations.get_time(1) - populations.get_time(0))
              << std::endl;

    return mio::success();
}

/** 
* Usage: lct_covid19_inspired scenario <data_dir> <save_dir> <start_date-year> <start_date-month> <start_date-day> 
*    <RelativeTransmissionNoSymptoms> <RiskOfInfectionFromSymptomatic> <scale_contacts> 
*    <scale_confirmed_cases> <npi_size> 
*
*   Specifying these command line arguments has the following implications:
*   - <data_dir> Directory with contact data, divi data and rki data used for the initialization. 
*       See the README for a specification of what data is required and how to download it.
*   - <save_dir> Directory where the simulation results should be stored.
*   - <start_date-year> <start_date-month> <start_date-day> (integers): Specification of the start date of the simulation.
*   - <RelativeTransmissionNoSymptoms> <RiskOfInfectionFromSymptomatic> Parameter values regarding the isolation of different infection states. 
*   - <scale_contacts> (double): Scale confirmed case data to incorporate a detection ratio.
*   - <scale_confirmed_cases> (double): Scale the contact data so that the simulation results for the number of daily new transmissions align with the 
*       extrapolated RKI data (for aggregated age groups).
*   - <npi_size> (double): Size of the NPI to be implemented from 25/20/2020. 0 means no NPIs, 1 means no contacts.
*   
*   All command line arguments are optional but it is beneficial to specify at least the first argument <data_dir> 
*   as the default is just an educated guess.
*
*   The numbers of subcompartments used in the LCT model is determined by the preprocessor macro NUM_SUBCOMPARTMENTS.
*   You can set the number via the flag -DNUM_SUBCOMPARTMENTS=... . 
*   If NUM_SUBCOMPARTMENTS is set to zero, an LCT model is used with numbers of subcompartments such that 
*   each corresponds to the approximate stay time in the compartment.
*/
int main(int argc, char** argv)
{
    std::string data_dir = "../../data";
    std::string save_dir = "";

    switch (argc) {
    case 11:
        params::npi_size = std::stod(argv[10]);
        [[fallthrough]];
    case 10:
        params::scale_confirmed_cases = std::stod(argv[9]);
        [[fallthrough]];
    case 9:
        params::scale_contacts = std::stod(argv[8]);
        [[fallthrough]];
    case 8:
        params::riskOfInfectionFromSymptomatic = std::stod(argv[7]);
        [[fallthrough]];
    case 7:
        params::relativeTransmissionNoSymptoms = std::stod(argv[6]);
        [[fallthrough]];
    case 6:
        params::start_date = mio::Date(std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]));
        [[fallthrough]];
    case 3:
        save_dir = argv[2];
        [[fallthrough]];
    case 2:
        data_dir = argv[1];
    }

    const std::string contact_data_dir   = data_dir + "/contacts/";
    const std::string infection_data_dir = data_dir + "/pydata/Germany/cases_all_age_ma7.json";
    const std::string divi_data_dir      = data_dir + "/pydata/Germany/germany_divi_all_dates.json";

    auto result = simulate(contact_data_dir, infection_data_dir, divi_data_dir, save_dir);
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}
