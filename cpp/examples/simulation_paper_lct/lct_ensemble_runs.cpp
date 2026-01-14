/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Lena Ploetzke, Rene Schmieding
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

namespace mio
{

template <class IOContext, class... States>
void serialize_internal(IOContext& io, const mio::lsecir::Model<States...>& m)
{
    auto obj = io.create_object("LctModel");
    obj.add_element("Parameters", m.parameters);
    obj.add_element("Populations", m.populations);
}

template <class IOContext, class... States>
IOResult<mio::lsecir::Model<States...>> deserialize_internal(IOContext& io, Tag<mio::lsecir::Model<States...>> tag)
{
    using M     = mio::lsecir::Model<States...>;
    auto obj    = io.expect_object("LctModel");
    auto params = obj.expect_element("Parameters", Tag<typename M::ParameterSet>{});
    auto pop    = obj.expect_element("Populations", Tag<typename M::Populations>{});
    return mio::apply(
        io,
        [](auto&& params_, auto&& pop_) {
            return M{pop_, params_};
        },
        params, pop);
}

template <class IOContext, class... States>
void serialize_internal(IOContext& io, const mio::LctPopulations<States...>& pop)
{
    auto obj = io.create_object("LctPopulations");
    obj.add_element("Populations", pop.array());
}

template <class IOContext, class... States>
IOResult<mio::LctPopulations<States...>> deserialize_internal(IOContext& io, Tag<mio::LctPopulations<States...>> tag)
{
    auto obj = io.expect_object("LctPopulations");
    auto pop = obj.expect_element("Populations", Tag<Eigen::Array<UncertainValue<ScalarType>, Eigen::Dynamic, 1>>{});
    return mio::apply(
        io,
        [](auto&& pop_) {
            mio::LctPopulations<States...> p;
            p.array() = pop_;
            return p;
        },
        pop);
}

} // namespace mio

#include "memilio/config.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/epidemiology/lct_infection_state.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/compartments/simulation.h"
#include "memilio/data/analyze_result.h"
#include "memilio/utils/base_dir.h"
#include "memilio/io/cli.h"
#include "memilio/io/result_io.h"
#include "memilio/io/io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/io/epi_data.h"

#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include <string>
#include <iostream>
#include <vector>
#include <type_traits>
#include <omp.h>
#include <mpi.h>

#include "memilio/utils/miompi.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/compartments/parameter_studies.h"

mio::UncertainValue<ScalarType> uncertain(ScalarType v)
{
    const double var = .1;
    return mio::UncertainValue<ScalarType>(v, mio::ParameterDistributionUniform(v * (1 - var), v * (1 + var)));
}

namespace params
{
constexpr int num_subcompartments = 5;
constexpr size_t num_groups       = 6;
int num_processes                 = 1;

// Define (age-resolved) parameters.
mio::Date start_date(2021, 01, 01);
const ScalarType tmax = 30;
const ScalarType dt   = 0.1;
const ScalarType t0   = 0;

const ScalarType seasonality              = 0.;
ScalarType relativeTransmissionNoSymptoms = 1.;
ScalarType riskOfInfectionFromSymptomatic = 0.3;
const ScalarType age_group_sizes[]        = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};
const ScalarType total_population         = 83155031.0;

const ScalarType transmissionProbabilityOnContact[] = {0.03, 0.06, 0.06, 0.06, 0.09, 0.175};

// Mean stay times
const ScalarType timeExposed[]            = {3.335, 3.335, 3.335, 3.335, 3.335, 3.335};
const ScalarType timeInfectedNoSymptoms[] = {2.74, 2.74, 2.565, 2.565, 2.565, 2.565};
const ScalarType timeInfectedSymptoms[]   = {7.02625, 7.02625, 7.0665, 6.9385, 6.835, 6.775};
const ScalarType timeInfectedSevere[]     = {5., 5., 5.925, 7.55, 8.5, 11.};
const ScalarType timeInfectedCritical[]   = {6.95, 6.95, 6.86, 17.36, 17.1, 11.6};
// Transition probabilities
const ScalarType recoveredPerInfectedNoSymptoms[] = {1 - 0.75, 1 - 0.75, 1 - 0.8, 1 - 0.8, 1 - 0.8, 1 - 0.8};
const ScalarType severePerInfectedSymptoms[]      = {0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225};
const ScalarType criticalPerSevere[]              = {0.075, 0.075, 0.075, 0.15, 0.3, 0.4};
const ScalarType deathsPerCritical[]              = {0.05, 0.05, 0.14, 0.14, 0.4, 0.6};

// Define default scalings that can be changed via command line.
ScalarType scale_confirmed_cases = 1.; ///< Scale confirmed case data to incorporate a detection ratio.
ScalarType scale_contacts =
    1.; ///< Scale the contact data so that the simulation results for the number of daily new transmissions align with the extrapolated RKI data (for aggregated age groups).
ScalarType npi_size = 0.; ///< Effect on the contacts of the NPI to be implemented from 25/10/2020 on.

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
 * @brief Add NPIs to a given contact matrix from 25/10/2020 onwards.
 *
 * Add NPIs to the contact matrix contact_matrices from 25/10/2020 onwards. 
 * The size of the NPI is specified in the variable npi_size. The NPI is applied to all contact locations equally.
 * 
 * @param[in,out] contact_matrices The contact matrices where the NPIs should be added to.
 */
// void set_npi_october(mio::ContactMatrixGroup& contact_matrices)
// {
//     using namespace params;
//     auto offset_npi = mio::SimulationTime(mio::get_offset_in_days(mio::Date(2020, 10, 25), start_date));
//     for (auto&& contact_location : contact_locations) {
//         contact_matrices[size_t(contact_location.first)].add_damping(
//             Eigen::MatrixXd::Constant(num_groups, num_groups, npi_size),
//             offset_npi); // no uncertain: internal type is MatrixXd
//     }
// }

/**
 * @brief Set the contact pattern from data files.
 *
 * The contacts are set using contact matrices from files in the data directory for different locations.
 * Non-pharmaceutical interventions (NPIs) influencing the ContactPatterns are set.
 * The contact matrices are scaled using the parameter scale_contacts such that the simulation results for the number
 *  of daily new transmissions aligns with the extrapolated RKI data.
 * 
 * @param[in] contact_data_dir Directory to files with minimum and baseline contact matrices.
 * @returns The contact matrix or any io errors that happen during reading of the input files.
 */
mio::IOResult<mio::UncertainContactMatrix<ScalarType>> get_contact_matrix(std::string contact_data_dir)
{
    using namespace params;
    auto contact_matrices = mio::ContactMatrixGroup<ScalarType>(contact_locations.size(), num_groups);

    // Load and set baseline contacts for each contact location.
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline, mio::read_mobility_plain(mio::path_join(
                                               contact_data_dir, ("baseline_" + contact_location.second + ".txt"))));
        contact_matrices[size_t(contact_location.first)].get_baseline() =
            scale_contacts * baseline; // no uncertain: is matrix
        contact_matrices[size_t(contact_location.first)].get_minimum() =
            Eigen::MatrixXd::Zero(num_groups, num_groups); // no uncertain: is 0
    }

    // Add NPIs to the contact matrices.
    mio::Date end_date     = mio::offset_date_by_days(start_date, (int)tmax);
    auto start_npi_october = mio::Date(2020, 10, 11);
    if ((int)(start_npi_october < end_date) & (int)(start_date < start_npi_october)) {
        auto offset_npi = mio::SimulationTime<ScalarType>(mio::get_offset_in_days(start_npi_october, start_date));
        for (auto&& contact_location : contact_locations) {
            contact_matrices[size_t(contact_location.first)].add_damping(
                Eigen::MatrixXd::Constant(num_groups, num_groups, 0.4),
                offset_npi); // no uncertain: internal type is MatrixXd
        }
    }

    // Add NPIs to the contact matrices.
    start_npi_october = mio::Date(2020, 10, 21);
    if ((int)(start_npi_october < end_date) & (int)(start_date < start_npi_october)) {
        auto offset_npi = mio::SimulationTime<ScalarType>(mio::get_offset_in_days(start_npi_october, start_date));
        for (auto&& contact_location : contact_locations) {
            contact_matrices[size_t(contact_location.first)].add_damping(
                Eigen::MatrixXd::Constant(num_groups, num_groups, 0.1),
                offset_npi); // no uncertain: internal type is MatrixXd
        }
    }
    return mio::success(mio::UncertainContactMatrix<ScalarType>(contact_matrices));
}

template <class Model>
Model draw_sample(const Model& model)
{
    auto copy = model;
    for (size_t group = 0; group < params::num_groups; group++) {
        copy.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[group].draw_sample();
        copy.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[group].draw_sample();
        copy.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[group].draw_sample();
        copy.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[group].draw_sample();
        copy.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[group].draw_sample();
        copy.parameters.template get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[group].draw_sample();

        copy.parameters.template get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[group].draw_sample();
        copy.parameters.template get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[group].draw_sample();
        copy.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[group].draw_sample();
        copy.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[group].draw_sample();
        copy.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[group].draw_sample();
        copy.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[group].draw_sample();
    }

    return copy;
}

template <class Model>
mio::IOResult<Model> initialize_lsecir(std::string data_dir)
{

    using namespace params;
    Model model;

    // Define parameters used for simulation and initialization.
    for (size_t group = 0; group < num_groups; group++) {
        model.parameters.template get<mio::lsecir::TimeExposed<ScalarType>>()[group] = uncertain(timeExposed[group]);
        model.parameters.template get<mio::lsecir::TimeInfectedNoSymptoms<ScalarType>>()[group] =
            uncertain(timeInfectedNoSymptoms[group]);
        model.parameters.template get<mio::lsecir::TimeInfectedSymptoms<ScalarType>>()[group] =
            uncertain(timeInfectedSymptoms[group]);
        model.parameters.template get<mio::lsecir::TimeInfectedSevere<ScalarType>>()[group] =
            uncertain(timeInfectedSevere[group]);
        model.parameters.template get<mio::lsecir::TimeInfectedCritical<ScalarType>>()[group] =
            uncertain(timeInfectedCritical[group]);
        model.parameters.template get<mio::lsecir::TransmissionProbabilityOnContact<ScalarType>>()[group] =
            uncertain(transmissionProbabilityOnContact[group]);

        model.parameters.template get<mio::lsecir::RelativeTransmissionNoSymptoms<ScalarType>>()[group] =
            uncertain(relativeTransmissionNoSymptoms);
        model.parameters.template get<mio::lsecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[group] =
            uncertain(riskOfInfectionFromSymptomatic);
        model.parameters.template get<mio::lsecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[group] =
            uncertain(recoveredPerInfectedNoSymptoms[group]);
        model.parameters.template get<mio::lsecir::SeverePerInfectedSymptoms<ScalarType>>()[group] =
            uncertain(severePerInfectedSymptoms[group]);
        model.parameters.template get<mio::lsecir::CriticalPerSevere<ScalarType>>()[group] =
            uncertain(criticalPerSevere[group]);
        model.parameters.template get<mio::lsecir::DeathsPerCritical<ScalarType>>()[group] =
            uncertain(deathsPerCritical[group]);
    }

    BOOST_OUTCOME_TRY(auto&& contact_matrix, get_contact_matrix(mio::path_join(data_dir, "Germany/contacts")));

    model.parameters.template get<mio::lsecir::ContactPatterns<ScalarType>>() =
        contact_matrix; // no uncertain: set by get_contact_matrix
    model.parameters.template get<mio::lsecir::Seasonality<ScalarType>>() = seasonality; // no uncertain: is 0
    model.parameters.template get<mio::lsecir::StartDay<ScalarType>>() =
        mio::get_day_in_year(start_date); // no uncertain: is date

    // Set initial values using reported data.
    BOOST_OUTCOME_TRY(auto&& rki_data, mio::read_confirmed_cases_data(
                                           mio::path_join(data_dir, "Germany/pydata/cases_all_age_ma7.json")));
    BOOST_OUTCOME_TRY(auto&& divi_data,
                      mio::read_divi_data(mio::path_join(data_dir, "Germany/pydata/germany_divi.json")));
    auto init =
        mio::lsecir::set_initial_values_from_reported_data<typename Model::Populations, mio::ConfirmedCasesDataEntry>(
            rki_data, model.populations, model.parameters, start_date,
            std::vector<ScalarType>(age_group_sizes, age_group_sizes + num_groups),
            std::vector<ScalarType>(num_groups, scale_confirmed_cases), divi_data);

    return mio::success(model);
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
 *   the number of subcompartments corresponds to the approximate stay time in the compartment.
 *
 * @param[in] data_dir Directory to data.
 * @param[in] save_dir Specifies the directory where the results should be stored. 
 *   Provide an empty string if the results should not be saved.
 * @returns Any io errors that happen during reading of the RKI file or files for contact matrices or saving the results.
 */
mio::IOResult<void> simulate(std::string save_dir, std::string data_dir, size_t num_ensemble_runs)
{
    mio::set_log_level(mio::LogLevel::off);

    using namespace params;
    if (mio::mpi::is_root()) {
        std::cout << "Realistic scenario with " << num_subcompartments << " subcompartments." << std::endl;
    }
    // Initialize (age-resolved) model.
    using InfState = mio::lsecir::InfectionState;
    // For 1.) : Define single LctState.
    using LctState = mio::LctInfectionState<ScalarType, InfState, 1, num_subcompartments, num_subcompartments,
                                            num_subcompartments, num_subcompartments, num_subcompartments, 1, 1>;
    // Decide which LCT model type should be used.
    using Model = mio::lsecir::Model<ScalarType, LctState, LctState, LctState, LctState, LctState, LctState>;

    Model model = initialize_lsecir<Model>(data_dir).value();
    // Set integrator of fifth order with fixed step size and perform simulation.
    // auto integrator =
    //     std::make_unique<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // // Choose dt_min = dt_max to get a fixed step size.
    // integrator->set_dt_min(dt);
    // integrator->set_dt_max(dt);

    mio::ParameterStudy parameter_study(model, params::t0, params::tmax, params::dt, num_ensemble_runs);
    ScalarType total_time = 0;
    if (mio::mpi::is_root()) {
        total_time -= omp_get_wtime();
    }
    parameter_study.run(
        [](auto&& params_model, ScalarType t0, ScalarType dt, size_t) {
            auto copy = params_model;
            return mio::Simulation<ScalarType, Model>(draw_sample<Model>(copy), t0, dt);
        },
        [&](auto results_model, auto&& run_idx) {
            auto interpolated_results = mio::interpolate_simulation_result(results_model.get_result());
            mio::unused(interpolated_results, run_idx);
        });
    if (mio::mpi::is_root()) {
        total_time += omp_get_wtime();
    }

    if (mio::mpi::is_root()) {
        std::cout << "{ \"Subcompartments\": " << num_subcompartments << ", " << std::endl;
        std::cout << "\"Processes\": " << num_processes << "," << std::endl;
        std::cout << "\"Num_ensemble_runs\": " << num_ensemble_runs << "," << std::endl;
        std::cout << "\"Time\": " << total_time << "\n}," << std::endl;
    }
    mio::set_log_level(mio::LogLevel::warn);
    mio::unused(save_dir);
    return mio::success();
}

/** 
* Usage: lct_covid19_inspired scenario <data_dir> <save_dir> <start_date-year> <start_date-month> <start_date-day> 
*    <RelativeTransmissionNoSymptoms> <RiskOfInfectionFromSymptomatic> <scale_contacts> 
*    <scale_confirmed_cases> <npi_size> 
*
*   Specifying these command line arguments has the following implications:
*   - <infection_data_dir> Directory with divi data and rki data used for the initialization. 
*       See the README for a specification of what data is required and how to download it.
*   - <contact_data_dir> Directory with contact data. Should be provided by the download of the MEmilio repository.
*   - <save_dir> Directory where the simulation results should be stored.
*   - <start_date-year> <start_date-month> <start_date-day> (integers): Specification of the start date of the simulation.
*   - <RelativeTransmissionNoSymptoms> <RiskOfInfectionFromSymptomatic> Parameter values regarding the isolation of different infection states. 
*   - <scale_contacts> (double): Scale confirmed case data to incorporate a detection ratio.
*   - <scale_confirmed_cases> (double): Scale the contact data so that the simulation results for the number of daily new transmissions align with the 
*       extrapolated RKI data (for aggregated age groups).
*   - <npi_size> (double): Effect of the NPI to be implemented from 25/10/2020 on, sensible values lie between 0 and 1 where 0 means no NPI and 1 means no contacts.
*   
*   All command line arguments are optional but it is beneficial to specify at least the first argument <data_dir> 
*   as the default is just an educated guess and save_dir as results will not be saved otherwise.
*
*   The numbers of subcompartments used in the LCT model is determined by the preprocessor macro NUM_SUBCOMPARTMENTS.
*   You can set the number via the flag -DNUM_SUBCOMPARTMENTS=... . 
*   If NUM_SUBCOMPARTMENTS is set to zero, an LCT model is used with numbers of subcompartments such that 
*   this number corresponds to the approximate stay time in the compartment.
*/
int main(int argc, char** argv)
{
    auto cli_parameters = mio::cli::ParameterSetBuilder()
                              .add<"ResultDirectory">(
                                  mio::path_join(mio::base_dir(), "cpp/examples/simulation_paper_lct/results_ensemble"))
                              .add<"DataDirectory">(mio::path_join(mio::base_dir(), "data"))
                              .add<"NumberEnsembleRuns">(100, {.alias = "nRun"})
                              .build();

    auto cli_result = mio::command_line_interface(argv[0], argc, argv, cli_parameters, {"ResultDirectory"});
    if (!cli_result) {
        std::cout << cli_result.error().message();
        return cli_result.error().code().value();
    }

    boost::filesystem::path res_dir(cli_parameters.get<"ResultDirectory">());
    boost::filesystem::create_directories(res_dir);

    mio::mpi::init();
    int size;
    MPI_Comm_size(mio::mpi::get_world(), &size);
    params::num_processes = size;

    auto result = simulate(cli_parameters.get<"ResultDirectory">(), cli_parameters.get<"DataDirectory">(),
                           cli_parameters.get<"NumberEnsembleRuns">());
    if (!result) {
        if (mio::mpi::is_root()) {
            printf("%s\n", result.error().formatted_message().c_str());
        }
        mio::mpi::finalize();
        return -1;
    }
    mio::mpi::finalize();

    return 0;
}
