/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Anna Wendler, Lena Ploetzke
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
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/utils/time_series.h"

#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/parameters.h"
#include "ide_secir/parameters_io.h"
#include "ide_secir/simulation.h"

#include "memilio/epidemiology/contact_matrix.h"
#include "ode_secir/model.h"
#include "ode_secir/infection_state.h"
#include "ode_secir/parameters.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "boost/filesystem.hpp"
#include <iomanip>
#include <string>
#include <map>
#include <iostream>

using Vector = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;

// Used parameters.
std::map<std::string, ScalarType> simulation_parameter = {{"t0", 0.},
                                                          {"dt", 0.01},
                                                          {"total_population", 83155031.},
                                                          {"total_confirmed_cases", 0.}, // set by RKI data
                                                          {"deaths", 0.}, // set by RKI data
                                                          {"TimeExposed", 4.5},
                                                          {"TimeInfectedNoSymptoms", 2.527617},
                                                          {"TimeInfectedSymptoms", 7.889900},
                                                          {"TimeInfectedSevere", 15.225278},
                                                          {"TimeInfectedCritical", 15.230258},
                                                          {"TransmissionProbabilityOnContact", 0.0733271},
                                                          {"RelativeTransmissionNoSymptoms", 1},
                                                          {"RiskOfInfectionFromSymptomatic", 0.3},
                                                          {"Seasonality", 0.},
                                                          {"InfectedSymptomsPerInfectedNoSymptoms", 0.793099},
                                                          {"SeverePerInfectedSymptoms", 0.078643},
                                                          {"CriticalPerSevere", 0.173176},
                                                          {"DeathsPerCritical", 0.387803},
                                                          {"scale_confirmed_cases", 1.},
                                                          {"scale_contacts", 1.},
                                                          {"lockdown_hard", 371 * 14 / (45 * 401.)}};

/**
 * Indices of contact matrix corresponding to locations where contacts occur.
 */
enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count,
};

/**
 * Different types of NPIs, used as DampingType.
 */
enum class Intervention
{
    Home,
    SchoolClosure,
    HomeOffice,
    GatheringBanFacilitiesClosure,
    PhysicalDistanceAndMasks,
    Count,
};

/**
 * Different levels of NPIs, used as DampingLevel.
 */
enum class InterventionLevel
{
    Main,
    PhysicalDistanceAndMasks,
    Holidays,
    Count,
};

// Map the ContactLocation%s to file names.
static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}};

/**
 * @brief Add NPIs to a given contact matrix from 01/10/2020 on.
 *
 * NPIs from the Paper "Assessment of effective mitigation ..." (doi: 10.1016/j.mbs.2021.108648) are used with slight 
 * modifications for a period of 45 days from 01/10/2020 on.
 * 
 * @param[in] contact_matrices The contact matrices to which the NPIs are applied to.
 * @param[in] start_date Start date of the simulation.
 * @param[in] lockdown_hard Proportion of counties for which a hard lockdown is implemented.
 */
void set_npi_october(mio::ContactMatrixGroup& contact_matrices, mio::Date start_date, ScalarType lockdown_hard)
{
    // ---------------------01/10/2020--------------------------------
    // NPIs from paper for October.
    auto offset_npi = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 10, 1), start_date));
    // For the beginning of the time period, we assume only half of the defined proportion of counties is in a hard lockdown.
    lockdown_hard = lockdown_hard / 2;
    // Contact reduction at home.
    ScalarType v = 0.2 * (1 - lockdown_hard) + lockdown_hard * 0.4;
    contact_matrices[size_t(ContactLocation::Home)].add_damping(Eigen::MatrixXd::Constant(1, 1, v),
                                                                mio::DampingLevel(int(InterventionLevel::Main)),
                                                                mio::DampingType(int(Intervention::Home)), offset_npi);
    // Home-Office + people stopped working.
    v = (0.2 + 0.025) * (1 - lockdown_hard) + lockdown_hard * (0.3 + 0.025);
    contact_matrices[size_t(ContactLocation::Work)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::Main)),
        mio::DampingType(int(Intervention::HomeOffice)), offset_npi);
    // GatheringBanFacilitiesClosure affects ContactLocation Other.
    v = 0. * (1 - lockdown_hard) + lockdown_hard * 0.2;
    contact_matrices[size_t(ContactLocation::Other)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::Main)),
        mio::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), offset_npi);
    // PhysicalDistanceAndMasks in all locations.
    v = 0.2 * (1 - lockdown_hard) + lockdown_hard * 0.4;
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
    /* We assume that the stricter NPIs of November defined in the paper are beginning about a week earlier, 
    which can be observed in the RKI data.
    Moreover the lockdown value of PhysicalDistanceAndMasks in the location school is assumed to apply for all counties.*/
    offset_npi = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 10, 24), start_date));
    // For the second half of the simulation, the proportion of counties in hard lockdown is increased to compensate for the lower proportion before.
    lockdown_hard = lockdown_hard * 5;
    // Contact reduction at home.
    v = 0.4 * (1 - lockdown_hard) + 0.6 * lockdown_hard;
    contact_matrices[size_t(ContactLocation::Home)].add_damping(Eigen::MatrixXd::Constant(1, 1, v),
                                                                mio::DampingLevel(int(InterventionLevel::Main)),
                                                                mio::DampingType(int(Intervention::Home)), offset_npi);
    // Home-Office + people stopped working.
    v = (0.2 + 0.05) * (1 - lockdown_hard) + lockdown_hard * (0.3 + 0.05);
    contact_matrices[size_t(ContactLocation::Work)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::Main)),
        mio::DampingType(int(Intervention::HomeOffice)), offset_npi);
    // GatheringBanFacilitiesClosure affects ContactLocation Other.
    v = 0.6 * (1 - lockdown_hard) + 0.8 * lockdown_hard;
    contact_matrices[size_t(ContactLocation::Other)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::Main)),
        mio::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), offset_npi);
    // PhysicalDistanceAndMasks in ContactLocation%s Home.
    v = 0.2 * (1 - lockdown_hard) + lockdown_hard * 0.4;
    contact_matrices[size_t(ContactLocation::Home)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
        mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), offset_npi);
    // PhysicalDistanceAndMasks in ContactLocation%s School.
    v = 0.2 * (1 - lockdown_hard) + lockdown_hard * 0.4;
    contact_matrices[size_t(ContactLocation::School)].add_damping(
        Eigen::MatrixXd::Constant(1, 1, v), mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
        mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), offset_npi);
    // PhysicalDistanceAndMasks in ContactLocation%s Work and Other.
    v = 0.4 * (1 - lockdown_hard) + lockdown_hard * 0.6;
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
 * @brief Set the contact pattern of parameters for a model without division in age groups.
 *
 * The contacts are calculated using contact matrices from files in the data directory for different locations.
 * Also set nonpharmaceutical interventions influencing the ContactPatterns used for simulation in the timeframe from start_date to end_date.
 * 
 * @param[in] data_dir Directory to files with minimum and baseline contact matrices.
 * @param[in] start_date Start date of the simulation.
 * @returns Any io errors that happen during reading of the input files.
 */
mio::IOResult<mio::ContactMatrixGroup> define_contact_matrices(const boost::filesystem::path& data_dir,
                                                               mio::Date start_date)
{
    // Files in data_dir are containing contact matrices with 6 agegroups. We use this to compute a contact pattern without division of age groups.
    // Age group sizes are calculated using table number 12411-04-02-4-B from www.regionalstatistik.de for the date 31.12.2020.
    const ScalarType age_group_sizes[] = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};
    const ScalarType total             = simulation_parameter["total_population"];
    const int numagegroups             = 6;

    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), 1);
    // Load and set minimum and baseline contacts for each contact location.
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));
        BOOST_OUTCOME_TRY(auto&& minimum,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("minimum_" + contact_location.second + ".txt")).string()));
        ScalarType base = 0;
        ScalarType min  = 0;
        for (int i = 0; i < numagegroups; i++) {
            for (int j = 0; j < numagegroups; j++) {
                // Calculate a weighted average according to the age group sizes of the total contacts.
                base += age_group_sizes[i] / total * baseline(i, j);
                min += age_group_sizes[i] / total * minimum(i, j);
            }
        }
        contact_matrices[size_t(contact_location.first)].get_baseline() =
            simulation_parameter["scale_contacts"] * Eigen::MatrixXd::Constant(1, 1, base);
        contact_matrices[size_t(contact_location.first)].get_minimum() =
            simulation_parameter["scale_contacts"] * Eigen::MatrixXd::Constant(1, 1, min);
    }

    // ----- Add NPIs to the contact matrices. -----
    // Set of NPIs for October.
    if (start_date == mio::Date(2020, 10, 1)) {
        set_npi_october(contact_matrices, start_date, simulation_parameter["lockdown_hard"]);
    }

    return mio::success(contact_matrices);
}

/**
 * @brief Set the contact pattern of parameters for a model without division in age groups without using the 
 * age-resolved contact_matrices.
 *
 * In case of only one age group the contact matrix reduces to a 1x1 matrix. 
 * Instead of using contact matrices from files in the data directory for different locations as in the function
 * define_contact_matrices(), we set the contact frequency to the value that we obtained using the function 
 * define_contact_matrices() above. Accordingly, we set the damping as in the function above to model the implementation
 * of NPIs on Oct 24, 2020. 
 *
 * @param[in] start_date Start date of the simulation.
 * @returns Any io errors that happen during reading of the input files.
 */
mio::IOResult<mio::ContactMatrixGroup> define_contact_matrices_simplified(mio::Date start_date)
{
    // Set of NPIs for October.
    auto contact_matrices  = mio::ContactMatrixGroup(1, 1);
    auto start_npi_october = start_date;
    if (start_npi_october == mio::Date(2020, 10, 1)) {

        contact_matrices[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 7.69129));

        auto offset_npi = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 10, 24), start_npi_october));
        contact_matrices[0].add_damping(0., mio::SimulationTime(0.1));
        contact_matrices[0].add_damping(1 - 3.51782 / 7.69129, offset_npi);
    }

    return mio::success(contact_matrices);
}

/**
* @brief Simulates using an IDE model. 
*
* @param[in] start_date Start date of the simulation 
* @param[in] simulation_time Duration of the simulation.
* @param[in] contact_matrices Contact matrices used. 
* @param[in] data_dir Directory to files with minimum and baseline contact matrices and reported data from RKI.
* @param[in] save_dir Directory where simulation results will be stored. 
* @returns Any IO errros that happen.
*/
mio::IOResult<std::vector<mio::TimeSeries<ScalarType>>>
simulate_ide_model(mio::Date start_date, ScalarType simulation_time, mio::ContactMatrixGroup contact_matrices,
                   const boost::filesystem::path& data_dir, std::string save_dir = "")

{
    // Initialize model.
    size_t num_agegroups = 1;
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_population =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups),
                                                         simulation_parameter["total_population"]);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), simulation_parameter["deaths"]);
    mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_confirmed_cases =
        mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups),
                                                         simulation_parameter["total_confirmed_cases"]);

    mio::isecir::Model model_ide(mio::TimeSeries<ScalarType>((int)mio::isecir::InfectionTransition::Count),
                                 total_population, deaths, num_agegroups, total_confirmed_cases);

    // Set working parameters.
    // Set TransitionDistributions.
    mio::ConstantFunction initialfunc(0);
    mio::StateAgeFunctionWrapper delaydistributioninit(initialfunc);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib((int)mio::isecir::InfectionTransition::Count,
                                                               delaydistributioninit);
    // ExposedToInfectedNoSymptoms
    mio::LognormSurvivalFunction survivalExposedToInfectedNoSymptoms(0.32459285, 0, 4.26907484);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms].set_state_age_function(
        survivalExposedToInfectedNoSymptoms);
    // InfectedNoSymptomsToInfectedSymptoms
    mio::LognormSurvivalFunction survivalInfectedNoSymptomsToInfectedSymptoms(0.71587510, 0, 0.85135303);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
        .set_state_age_function(survivalInfectedNoSymptomsToInfectedSymptoms);
    // InfectedNoSymptomsToRecovered
    mio::LognormSurvivalFunction survivalInfectedNoSymptomsToRecovered(0.24622068, 0, 7.7611400);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered].set_state_age_function(
        survivalInfectedNoSymptomsToRecovered);
    // InfectedSymptomsToInfectedSevere
    mio::LognormSurvivalFunction survivalInfectedSymptomsToInfectedSevere(0.66258947, 0, 5.29920733);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere].set_state_age_function(
        survivalInfectedSymptomsToInfectedSevere);
    // InfectedSymptomsToRecovered
    mio::LognormSurvivalFunction survivalInfectedSymptomsToRecovered(0.24622068, 0, 7.76114000);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered].set_state_age_function(
        survivalInfectedSymptomsToRecovered);
    // InfectedSevereToInfectedCritical
    mio::LognormSurvivalFunction survivalInfectedSevereToInfectedCritical(1.010767652595, 0, 0.90000000);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical].set_state_age_function(
        survivalInfectedSevereToInfectedCritical);
    // InfectedSevereToRecovered
    mio::LognormSurvivalFunction survivalInfectedSevereToRecovered(0.33816427, 0, 17.09411753);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered].set_state_age_function(
        survivalInfectedSevereToRecovered);
    // InfectedCriticalToDead
    mio::LognormSurvivalFunction survivalInfectedCriticalToDead(0.42819924, 0, 9.76267505);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead].set_state_age_function(
        survivalInfectedCriticalToDead);
    // InfectedCriticalToRecovered
    mio::LognormSurvivalFunction survivalInfectedCriticalToRecovered(0.33816427, 0, 17.09411753);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_state_age_function(
        survivalInfectedCriticalToRecovered);

    model_ide.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // Set other parameters.
    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 1.);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
        simulation_parameter["InfectedSymptomsPerInfectedNoSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
        1 - simulation_parameter["InfectedSymptomsPerInfectedNoSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
        simulation_parameter["SeverePerInfectedSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)] =
        1 - simulation_parameter["SeverePerInfectedSymptoms"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
        simulation_parameter["CriticalPerSevere"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)] =
        1 - simulation_parameter["CriticalPerSevere"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] =
        simulation_parameter["DeathsPerCritical"];
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)] =
        1 - simulation_parameter["DeathsPerCritical"];

    model_ide.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    model_ide.parameters.get<mio::isecir::ContactPatterns>() =
        mio::UncertainContactMatrix<ScalarType>(contact_matrices);

    mio::ConstantFunction constfunc(simulation_parameter["TransmissionProbabilityOnContact"]);
    mio::StateAgeFunctionWrapper StateAgeFunctionWrapperide(constfunc);
    model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(StateAgeFunctionWrapperide);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RelativeTransmissionNoSymptoms"]);
    model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(StateAgeFunctionWrapperide);
    StateAgeFunctionWrapperide.set_distribution_parameter(simulation_parameter["RiskOfInfectionFromSymptomatic"]);
    model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(StateAgeFunctionWrapperide);

    model_ide.set_tol_for_support_max(1e-6);

    // Set initial flows according to RKI data.
    std::string path_rki = mio::path_join((data_dir / "pydata" / "Germany").string(), "cases_all_germany_ma7.json");

    mio::IOResult<void> init_flows = set_initial_flows(model_ide, simulation_parameter["dt"], path_rki, start_date,
                                                       simulation_parameter["scale_confirmed_cases"]);

    model_ide.check_constraints(simulation_parameter["dt"]);

    // Simulate.
    mio::isecir::Simulation sim(model_ide, simulation_parameter["dt"]);

    sim.advance(simulation_time);

    // Save results.
    if (!save_dir.empty()) {
        std::string tmax_string  = std::to_string(simulation_time);
        std::string dt_string    = std::to_string(simulation_parameter["dt"]);
        std::string filename_ide = save_dir + "ide_" + std::to_string(start_date.year) + "-" +
                                   std::to_string(start_date.month) + "-" + std::to_string(start_date.day) + "_" +
                                   tmax_string.substr(0, tmax_string.find(".")) + "_" +
                                   dt_string.substr(0, dt_string.find(".") + 5);

        std::string filename_ide_flows = filename_ide + "_flows.h5";
        mio::IOResult<void> save_result_status_f =
            mio::save_result({sim.get_transitions()}, {0}, 1, filename_ide_flows);
        std::string filename_ide_compartments = filename_ide + "_compartments.h5";
        mio::IOResult<void> save_result_status_c =
            mio::save_result({sim.get_result()}, {0}, 1, filename_ide_compartments);
    }

    // Return vector with simulation results.
    std::vector<mio::TimeSeries<ScalarType>> result = {sim.get_result(), sim.get_transitions()};
    return mio::success(result);
}

/**
* @brief Simulates using an ODE model. 
*
* We need intial values for the compartments as input. With this, we can make the starting conditions equivalent to an
* the previously simulated results using an IDE model. 
*
* @param[in] start_date Start date of the simulation 
* @param[in] simulation_time Duration of the simulation.
* @param[in] init_compartments Vector containing initial values for compartments. 
* @param[in] contact_matrices Contact matrices used. 
* @param[in] save_dir Directory where simulation results will be stored. 
* @returns Any IO errros that happen.
*/
mio::IOResult<void> simulate_ode_model(mio::Date start_date, ScalarType simulation_time, Vector init_compartments,
                                       mio::ContactMatrixGroup contact_matrices, std::string save_dir = "")
{
    // Use ODE FlowModel.
    mio::osecir::Model model_ode(1);

    // Set working parameters.
    model_ode.parameters.get<mio::osecir::TimeExposed<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeExposed"];
    model_ode.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeInfectedNoSymptoms"];
    model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeInfectedSymptoms"];
    model_ode.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeInfectedSevere"];
    model_ode.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeInfectedCritical"];

    // Set probabilities that determine proportion between compartments.
    model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        1 - simulation_parameter["InfectedSymptomsPerInfectedNoSymptoms"];
    model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["SeverePerInfectedSymptoms"];
    model_ode.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["CriticalPerSevere"];
    model_ode.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["DeathsPerCritical"];

    // Further model parameters.
    model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TransmissionProbabilityOnContact"];
    model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["RelativeTransmissionNoSymptoms"];
    model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["RiskOfInfectionFromSymptomatic"];
    // Choose TestAndTraceCapacity very large so that riskFromInfectedSymptomatic = RiskOfInfectionFromSymptomatic.
    model_ode.parameters.get<mio::osecir::TestAndTraceCapacity<ScalarType>>() = std::numeric_limits<ScalarType>::max();
    // Choose ICUCapacity very large so that CriticalPerSevereAdjusted = CriticalPerSevere and deathsPerSevereAdjusted = 0.
    model_ode.parameters.get<mio::osecir::ICUCapacity<ScalarType>>() = std::numeric_limits<ScalarType>::max();

    // Set Seasonality=0 so that cont_freq_eff is equal to contact_matrix.
    model_ode.parameters.set<mio::osecir::Seasonality<ScalarType>>(simulation_parameter["Seasonality"]);

    model_ode.parameters.get<mio::osecir::ContactPatterns<ScalarType>>() =
        mio::UncertainContactMatrix<ScalarType>(contact_matrices);

    // Use mio::isecir::InfectionState when accessing init_compartments since this is computed using the IDE model.
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible}] =
        init_compartments[int(mio::isecir::InfectionState::Susceptible)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}] =
        init_compartments[int(mio::isecir::InfectionState::Exposed)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}] =
        init_compartments[int(mio::isecir::InfectionState::InfectedNoSymptoms)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] =
        init_compartments[int(mio::isecir::InfectionState::InfectedSymptoms)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}] = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}] =
        init_compartments[int(mio::isecir::InfectionState::InfectedSevere)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}] =
        init_compartments[int(mio::isecir::InfectionState::InfectedCritical)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}] =
        init_compartments[int(mio::isecir::InfectionState::Recovered)];
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}] =
        init_compartments[int(mio::isecir::InfectionState::Dead)];

    model_ode.check_constraints();

    // Set integrator and fix step size.
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    integrator->set_dt_min(simulation_parameter["dt"]);
    integrator->set_dt_max(simulation_parameter["dt"]);

    // Simulate.
    std::vector<mio::TimeSeries<ScalarType>> results_ode = mio::osecir::simulate_flows<ScalarType>(
        simulation_parameter["t0"], simulation_time, simulation_parameter["dt"], model_ode, integrator);

    // Save results.
    if (!save_dir.empty()) {
        std::string simulation_time_string = std::to_string(simulation_time);
        std::string dt_string              = std::to_string(simulation_parameter["dt"]);
        std::string filename_ode           = save_dir + "ode_" + std::to_string(start_date.year) + "-" +
                                   std::to_string(start_date.month) + "-" + std::to_string(start_date.day) + "_" +
                                   simulation_time_string.substr(0, simulation_time_string.find(".")) + "_" +
                                   dt_string.substr(0, dt_string.find(".") + 5);

        std::string filename_ode_flows           = filename_ode + "_flows.h5";
        mio::IOResult<void> save_result_status_f = mio::save_result({results_ode[1]}, {0}, 1, filename_ode_flows);
        std::string filename_ode_compartments    = filename_ode + "_compartments.h5";
        mio::IOResult<void> save_result_status_c =
            mio::save_result({results_ode[0]}, {0}, 1, filename_ode_compartments);
    }

    return mio::success();
}

int main(int argc, char** argv)
{
    // Paths are valid if file is executed e.g. in memilio/build/bin.
    std::string data_dir_string = "../../data";
    std::string save_dir        = "../../data/simulation_results/covid_scenario/";

    mio::Date start_date(2020, 10, 01);

    ScalarType simulation_time = 45;

    if (argc == 9) {

        data_dir_string                        = argv[1];
        save_dir                               = argv[2];
        start_date                             = mio::Date(std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]));
        simulation_time                        = std::stod(argv[6]);
        simulation_parameter["dt"]             = std::stod(argv[7]);
        simulation_parameter["scale_contacts"] = std::stod(argv[8]);

        std::cout << std::setprecision(10) << "Contact scaling: " << simulation_parameter["scale_contacts"]
                  << std::endl;
    }

    // Make folder if not existent yet.
    boost::filesystem::path dir(save_dir);
    boost::filesystem::create_directories(dir);

    const boost::filesystem::path data_dir = data_dir_string;

    // Set contact matrices.
    mio::ContactMatrixGroup contact_matrices = define_contact_matrices(data_dir, start_date).value();
    // mio::ContactMatrixGroup contact_matrices = define_contact_matrices_simplified(start_date).value();

    // Run IDE simulation.
    auto result_ide = simulate_ide_model(start_date, simulation_time, contact_matrices, data_dir, save_dir);
    if (!result_ide) {
        printf("%s\n", result_ide.error().formatted_message().c_str());
        return -1;
    }

    // Use results from IDE simulation as initialization for ODE model.
    Vector init_compartments = result_ide.value()[0].get_value(0);

    // Run ODE simulation.
    auto result_ode = simulate_ode_model(start_date, simulation_time, init_compartments, contact_matrices, save_dir);
    if (!result_ode) {
        printf("%s\n", result_ode.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}
