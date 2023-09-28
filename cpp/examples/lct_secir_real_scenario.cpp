/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
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

#include "boost/outcome/try.hpp"
#include "lct_secir/model.h"
#include "lct_secir/infection_state.h"
#include "lct_secir/parameters.h"
#include "lct_secir/simulation.h"
#include "lct_secir/parameters_io.h"

#include "memilio/config.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/utils/time_series.h"
#include "memilio/math/eigen.h"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "memilio/utils/uncertain_value.h"
#include <iostream>

#include "boost/filesystem.hpp"

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

/**
 * @brief Set the contact pattern of parameters for a Model without division in age groups.
 *
 *  The contacts are calculated using contact matrices from files in the data directory for different locations.
 * 
 * @param data_dir Directory to files with contact matrices.
 * @param parameters Object that the contact pattern will be added to.
 * @returns Any io errors that happen during reading of the files.
 */
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::lsecir::Parameters& parameters)
{
    // Files in data_dir are containing contact matrices with 6 agegroups. We use this to compute a contact pattern without division of age groups.
    // Age group sizes are calculated using table number 12411-04-02-4-B from www.regionalstatistik.de for the date 31.12.2020.
    const double age_group_sizes[] = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};
    const int total                = 83155031.0;
    const int numagegroups         = 6;

    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), 1);
    // Load and set minimum and baseline contacts for each contact location.
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(baseline,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));
        BOOST_OUTCOME_TRY(minimum,
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
        contact_matrices[size_t(contact_location.first)].get_baseline() = Eigen::MatrixXd::Constant(1, 1, base);
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::MatrixXd::Constant(1, 1, min);
    }
    // Set ContactPatterns in parameters.
    parameters.get<mio::lsecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrices);

    return mio::success();
}

/**
 * @brief Set Nonpharmaceutical Interventions influencing the ContactPatterns used for simulation in a specified timeframe.
 * @param start_date Start date of the simulation.
 * @param end_date End date of the simulation.
 * @param parameters Object that the NPIs will be added to.
 */
void set_npis(mio::Date start_date, mio::Date end_date, mio::lsecir::Parameters& parameters)
{
    auto& contacts         = parameters.get<mio::lsecir::ContactPatterns>();
    auto& contact_dampings = contacts.get_dampings();

    auto v = mio::UncertainValue(0.1);

    // Set of NPIs for June.
    auto start_npi_june = mio::Date(2020, 6, 1);
    if (start_npi_june < end_date) {
        auto offset_june = mio::SimulationTime(mio::get_offset_in_days(start_npi_june, start_date));

        // Contact reduction at home.
        v = 0.1;
        contact_dampings.push_back(mio::DampingSampling(
            v, mio::DampingLevel(int(InterventionLevel::Main)), mio::DampingType(int(Intervention::Home)), offset_june,
            {size_t(ContactLocation::Home)}, Eigen::VectorXd::Constant(1, 1.0)));
        // Home-Office + people stopped working.
        v = 0.25 + 0.025;
        contact_dampings.push_back(mio::DampingSampling(
            v, mio::DampingLevel(int(InterventionLevel::Main)), mio::DampingType(int(Intervention::HomeOffice)),
            offset_june, {size_t(ContactLocation::Work)}, Eigen::VectorXd::Constant(1, 1.0)));
        // GatheringBanFacilitiesClosure affects ContactLocation Other.
        v = 0.1;
        contact_dampings.push_back(
            mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::Main)),
                                 mio::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), offset_june,
                                 {size_t(ContactLocation::Other)}, Eigen::VectorXd::Constant(1, 1.0)));
        // PhysicalDistanceAndMasks in all locations.
        v = 0.1;
        contact_dampings.push_back(
            mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                 mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), offset_june,
                                 {size_t(ContactLocation::Home), size_t(ContactLocation::School),
                                  size_t(ContactLocation::Work), size_t(ContactLocation::Other)},
                                 Eigen::VectorXd::Constant(1, 1.0)));
        // Schools closed.
        v = 0.5;
        contact_dampings.push_back(mio::DampingSampling(
            v, mio::DampingLevel(int(InterventionLevel::Main)), mio::DampingType(int(Intervention::SchoolClosure)),
            offset_june, {size_t(ContactLocation::School)}, Eigen::VectorXd::Constant(1, 1.0)));
        auto school_reopen_time = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 6, 15), start_date));
        // School fully reopened.
        v = 0.0;
        contact_dampings.push_back(mio::DampingSampling(
            v, mio::DampingLevel(int(InterventionLevel::Main)), mio::DampingType(int(Intervention::SchoolClosure)),
            school_reopen_time, {size_t(ContactLocation::School)}, Eigen::VectorXd::Constant(1, 1.0)));
    }

    // Set of NPIs for October.
    auto start_npi_october = mio::Date(2020, 10, 1);
    if (start_npi_october < end_date) {
        auto start_autumn = mio::SimulationTime(mio::get_offset_in_days(start_npi_october, start_date));

        // Contact reduction at home.
        v = 0.3;
        contact_dampings.push_back(mio::DampingSampling(
            v, mio::DampingLevel(int(InterventionLevel::Main)), mio::DampingType(int(Intervention::Home)), start_autumn,
            {size_t(ContactLocation::Home)}, Eigen::VectorXd::Constant(1, 1.0)));
        // PhysicalDistanceAndMasks in all locations.
        v = 0.3;
        contact_dampings.push_back(
            mio::DampingSampling(v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                 mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), start_autumn,
                                 {size_t(ContactLocation::Home), size_t(ContactLocation::School),
                                  size_t(ContactLocation::Work), size_t(ContactLocation::Other)},
                                 Eigen::VectorXd::Constant(1, 1.0)));
        // Schools closed because of autumn break (use values for NRW).
        v                       = 1.;
        auto autumn_break_begin = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 10, 12), start_date));
        contact_dampings.push_back(mio::DampingSampling(
            v, mio::DampingLevel(int(InterventionLevel::Main)), mio::DampingType(int(Intervention::SchoolClosure)),
            autumn_break_begin, {size_t(ContactLocation::School)}, Eigen::VectorXd::Constant(1, 1.0)));
        // School fully reopened because of ending autumn break.
        auto autumn_break_end = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 10, 24), start_date));
        v                     = 0.0;
        contact_dampings.push_back(mio::DampingSampling(
            v, mio::DampingLevel(int(InterventionLevel::Main)), mio::DampingType(int(Intervention::SchoolClosure)),
            autumn_break_end, {size_t(ContactLocation::School)}, Eigen::VectorXd::Constant(1, 1.0)));
    }

    //Set of NPIs for November.
    auto start_npi_november = mio::Date(2020, 11, 1);
    if (start_npi_november < end_date) {

        auto start_autumn_lockdown = mio::SimulationTime(mio::get_offset_in_days(start_npi_november, start_date));
        // Contact reduction at home.
        v = 0.5;
        contact_dampings.push_back(mio::DampingSampling(
            v, mio::DampingLevel(int(InterventionLevel::Main)), mio::DampingType(int(Intervention::Home)),
            start_autumn_lockdown, {size_t(ContactLocation::Home)}, Eigen::VectorXd::Constant(1, 1.0)));
        // Home-Office + people stopped working.
        v = 0.25 + 0.05;
        contact_dampings.push_back(mio::DampingSampling(
            v, mio::DampingLevel(int(InterventionLevel::Main)), mio::DampingType(int(Intervention::HomeOffice)),
            start_autumn_lockdown, {size_t(ContactLocation::Work)}, Eigen::VectorXd::Constant(1, 1.0)));
        // GatheringBanFacilitiesClosure affects ContactLocation Other.
        v = 0.7;
        contact_dampings.push_back(mio::DampingSampling(
            v, mio::DampingLevel(int(InterventionLevel::Main)),
            mio::DampingType(int(Intervention::GatheringBanFacilitiesClosure)), start_autumn_lockdown,
            {size_t(ContactLocation::Other)}, Eigen::VectorXd::Constant(1, 1.0)));
        // PhysicalDistanceAndMasks in ContactLocation%s School and Home.
        v = 0.3;
        contact_dampings.push_back(mio::DampingSampling(
            v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), start_autumn_lockdown,
            {size_t(ContactLocation::Home), size_t(ContactLocation::School)}, Eigen::VectorXd::Constant(1, 1.0)));
        // PhysicalDistanceAndMasks in ContactLocation%s Work and Other.
        v = 0.5;
        contact_dampings.push_back(mio::DampingSampling(
            v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)), start_autumn_lockdown,
            {size_t(ContactLocation::Work), size_t(ContactLocation::Other)}, Eigen::VectorXd::Constant(1, 1.0)));
    }
    parameters.get<mio::lsecir::ContactPatterns>().make_matrix();
}

/**
 * @brief Performs a simulation of a real scenario with an LCT and an ODE model.
 *
 * @param path Path of the RKI file that should be used to compute initial values for simulations.
 * @param save_result Specifies if the results should be saved.
 * @param print_result Specifies if the results should be printed.
 * @returns Any io errors that happen during reading of the RKI file or files for contact matrices.
 */
mio::IOResult<void> simulate(std::string const& path, bool save_result = true, bool print_result = false)
{
    // Set values needed for initialization.
    ScalarType total_population = 83155031.0;
    auto start_date             = mio::Date(2020, 10, 15);
    auto end_date               = mio::Date(2020, 11, 15);

    // Define parameters used for simulation and initialization.
    mio::lsecir::Parameters parameters;
    parameters.get<mio::lsecir::TimeExposed>()                      = 3.335;
    parameters.get<mio::lsecir::TimeInfectedNoSymptoms>()           = 3.31331;
    parameters.get<mio::lsecir::TimeInfectedSymptoms>()             = 6.94547;
    parameters.get<mio::lsecir::TimeInfectedSevere>()               = 11.634346;
    parameters.get<mio::lsecir::TimeInfectedCritical>()             = 17.476959;
    parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.0733271;

    parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 1;
    parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.3;
    parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.206901;
    parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.0786429;
    parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.173176;
    parameters.get<mio::lsecir::DeathsPerCritical>()              = 0.217177;

    auto status = set_contact_matrices("../../data", parameters);
    set_npis(start_date, end_date, parameters);
    // TODO: if simulated cases are higher than expected maybe the autumn break is missing in end of october.

    // Define number of subcompartments.
    std::vector<int> vec_subcompartments((int)mio::lsecir::InfectionStateBase::Count, 1);
    // Use subcompartments with a soujourn time of approximately one day in each subcompartment.
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::Exposed] =
        round(parameters.get<mio::lsecir::TimeExposed>());
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] =
        round(parameters.get<mio::lsecir::TimeInfectedNoSymptoms>());
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSymptoms] =
        round(parameters.get<mio::lsecir::TimeInfectedSymptoms>());
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedSevere] =
        round(parameters.get<mio::lsecir::TimeInfectedSevere>());
    // Both realistic distributions for times corresponding to InfectedCritical of the IDE model are exponential distributions.
    vec_subcompartments[(int)mio::lsecir::InfectionStateBase::InfectedCritical] = 1;
    mio::lsecir::InfectionState infectionState(vec_subcompartments);

    // Calculate initial value vector for subcompartments with RKI data.
    BOOST_OUTCOME_TRY(init_subcompartments,
                      mio::lsecir::get_initial_data_from_file(path, start_date, infectionState, std::move(parameters),
                                                              total_population, 2.));
    // Sum subcompartments to get initial values ​​for the ODE model.
    Eigen::VectorXd init_base((int)mio::lsecir::InfectionStateBase::Count);
    for (int i = 0; i < (int)mio::lsecir::InfectionStateBase::Count; i++) {
        init_base[i] =
            init_subcompartments
                .segment(Eigen::Index(infectionState.get_firstindex(i)), Eigen::Index(infectionState.get_number(i)))
                .sum();
    }

    // Initialize LCT model and perform simulation.
    mio::lsecir::Model model_lct(std::move(init_subcompartments), infectionState, std::move(parameters));
    mio::TimeSeries<ScalarType> result_lct = mio::lsecir::simulate(
        0, mio::get_offset_in_days(end_date, start_date), 0.1, model_lct,
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>(1e-10, 1e-5, 0,
                                                                                                         0.1));
    // Calculate result without division in subcompartments.
    mio::TimeSeries<ScalarType> populations_lct = model_lct.calculate_populations(result_lct);

    // Initialize ODE model and perform simulation.
    mio::lsecir::InfectionState infectionState_ode(std::vector<int>((int)mio::lsecir::InfectionStateBase::Count, 1));
    mio::lsecir::Model model_ode(std::move(init_base), infectionState_ode, std::move(parameters));
    mio::TimeSeries<ScalarType> result_ode = mio::lsecir::simulate(
        0, mio::get_offset_in_days(end_date, start_date), 0.1, model_ode,
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>(1e-10, 1e-5, 0,
                                                                                                         0.1));

    if (print_result) {
        // Print results.
        std::cout << "Result LCT model:" << std::endl;
        mio::lsecir::print_TimeSeries(populations_lct, model_lct.get_heading_CompartmentsBase());
        std::cout << "Result ODE model:" << std::endl;
        mio::lsecir::print_TimeSeries(result_ode, model_ode.get_heading_CompartmentsBase());
    }
    if (save_result) {
        // Save results.
        auto save_result_status = mio::save_result({populations_lct}, {0}, 1, "real_lct_2020_10_15.h5");
        save_result_status      = mio::save_result({result_ode}, {0}, 1, "real_ode_2020_10_15.h5");
    }
    return mio::success();
}

int main()
{
    // Path is valid if file is executed eg in memilio/build/bin
    auto result = simulate("../../data/pydata/Germany/cases_all_germany_ma7.json");
    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }
    return 0;
}