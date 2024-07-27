/*
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Wadim Koslow, Martin Kühn, David Kerkmann
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

#include "memilio/compartments/parameter_studies.h"
#include "memilio/geography/regions.h"
#include "memilio/io/epi_data.h"
#include "memilio/io/result_io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/miompi.h"
#include "memilio/utils/random_number_generator.h"
#include "ode_secir/parameters_io.h"
#include "ode_secir/parameter_space.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/stl_util.h"
#include "boost/filesystem.hpp"
#include <cstdio>
#include <iomanip>

namespace fs = boost::filesystem;

/**
 * indices of contact matrix corresponding to locations where contacts occur.
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
 * Set a value and distribution of an UncertainValue.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution.
 * @param p uncertain value to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void assign_uniform_distribution(mio::UncertainValue& p, double min, double max)
{
    p = mio::UncertainValue(0.5 * (max + min));
    p.set_distribution(mio::ParameterDistributionUniform(min, max));
}

/**
 * Set a value and distribution of an array of UncertainValues.
 * Assigns average of min[i] and max[i] as a value and UNIFORM(min[i], max[i]) as a distribution for
 * each element i of the array.
 * @param array array of UncertainValues to set.
 * @param min minimum of distribution for each element of array.
 * @param max minimum of distribution for each element of array.
 */
template <size_t N>
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>& array,
                                       const double (&min)[N], const double (&max)[N])
{
    assert(N == array.numel());
    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(N); ++i) {
        assign_uniform_distribution(array[i], min[size_t(i)], max[size_t(i)]);
    }
}

/**
 * Set a value and distribution of an array of UncertainValues.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution to every element of the array.
 * @param array array of UncertainValues to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>& array, double min,
                                       double max)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max);
    }
}

/**
 * Set epidemiological parameters of Sars-CoV-2 for a immune-naive
 * population and wild type variant.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_covid_parameters(mio::osecir::Parameters& params)
{
    //times
    const double incubationTime            = 5.2;
    const double serialIntervalMin         = 1.1;
    const double serialIntervalMax         = 1.1;
    const double timeInfectedSymptomsMin[] = {5.6255, 5.6255, 5.6646, 5.5631, 5.501, 5.465};
    const double timeInfectedSymptomsMax[] = {8.427, 8.427, 8.4684, 8.3139, 8.169, 8.085};
    const double timeInfectedSevereMin[]   = {3.925, 3.925, 4.85, 6.4, 7.2, 9.};
    const double timeInfectedSevereMax[]   = {6.075, 6.075, 7., 8.7, 9.8, 13.};
    const double timeInfectedCriticalMin[] = {4.95, 4.95, 4.86, 14.14, 14.4, 10.};
    const double timeInfectedCriticalMax[] = {8.95, 8.95, 8.86, 20.58, 19.8, 13.2};

    array_assign_uniform_distribution(params.get<mio::osecir::TimeExposed>(), incubationTime, incubationTime);
    array_assign_uniform_distribution(params.get<mio::osecir::TimeInfectedNoSymptoms>(), serialIntervalMin,
                                      serialIntervalMax);
    array_assign_uniform_distribution(params.get<mio::osecir::TimeInfectedSymptoms>(), timeInfectedSymptomsMin,
                                      timeInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecir::TimeInfectedSevere>(), timeInfectedSevereMin,
                                      timeInfectedSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecir::TimeInfectedCritical>(), timeInfectedCriticalMin,
                                      timeInfectedCriticalMax);

    //probabilities
    const double transmissionProbabilityOnContactMin[] = {0.02, 0.05, 0.05, 0.05, 0.08, 0.15};
    const double transmissionProbabilityOnContactMax[] = {0.04, 0.07, 0.07, 0.07, 0.10, 0.20};
    const double relativeTransmissionNoSymptomsMin     = 1;
    const double relativeTransmissionNoSymptomsMax     = 1;
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    const double riskOfInfectionFromSymptomaticMin    = 0.1;
    const double riskOfInfectionFromSymptomaticMax    = 0.3;
    const double maxRiskOfInfectionFromSymptomaticMin = 0.3;
    const double maxRiskOfInfectionFromSymptomaticMax = 0.5;
    const double recoveredPerInfectedNoSymptomsMin[]  = {0.5, 0.45, 0.4, 0.35, 0.175, 0.05};
    const double recoveredPerInfectedNoSymptomsMax[]  = {0.5, 0.45, 0.4, 0.35, 0.175, 0.05};
    const double severePerInfectedSymptomsMin[]       = {0.03, 0.03, 0.04, 0.17, 0.05, 0.08};
    const double severePerInfectedSymptomsMax[]       = {0.03, 0.03, 0.04, 0.17, 0.05, 0.08};
    const double criticalPerSevereMin[]               = {0.10, 0.11, 0.12, 0.14, 0.33, 0.8};
    const double criticalPerSevereMax[]               = {0.10, 0.11, 0.12, 0.14, 0.33, 0.8};
    const double deathsPerCriticalMin[]               = {0.12, 0.13, 0.15, 0.29, 0.40, 0.48};
    const double deathsPerCriticalMax[]               = {0.12, 0.13, 0.15, 0.29, 0.40, 0.48};

    array_assign_uniform_distribution(params.get<mio::osecir::TransmissionProbabilityOnContact>(),
                                      transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax);
    array_assign_uniform_distribution(params.get<mio::osecir::RelativeTransmissionNoSymptoms>(),
                                      relativeTransmissionNoSymptomsMin, relativeTransmissionNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecir::RiskOfInfectionFromSymptomatic>(),
                                      riskOfInfectionFromSymptomaticMin, riskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic>(),
                                      maxRiskOfInfectionFromSymptomaticMin, maxRiskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>(),
                                      recoveredPerInfectedNoSymptomsMin, recoveredPerInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecir::SeverePerInfectedSymptoms>(),
                                      severePerInfectedSymptomsMin, severePerInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecir::CriticalPerSevere>(), criticalPerSevereMin,
                                      criticalPerSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecir::DeathsPerCritical>(), deathsPerCriticalMin,
                                      deathsPerCriticalMax);

    //sasonality
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.3;

    assign_uniform_distribution(params.get<mio::osecir::Seasonality>(), seasonality_min, seasonality_max);

    return mio::success();
}

static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}};

/**
 * Set contact matrices.
 * Reads contact matrices from files in the data directory.
 * @param data_dir data directory.
 * @param params Object that the contact matrices will be added to.
 * @returns any io errors that happen during reading of the files.
 */
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::osecir::Parameters& params)
{
    //TODO: io error handling
    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));
        BOOST_OUTCOME_TRY(auto&& minimum,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("minimum_" + contact_location.second + ".txt")).string()));
        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = minimum;
    }
    params.get<mio::osecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrices);

    return mio::success();
}

/**
 * Create the input graph for the parameter study.
 * Reads files from the data directory.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param data_dir data directory.
 * @returns created graph or any io errors that happen during reading of the files.
 */
mio::IOResult<std::vector<mio::osecir::Model>> get_graph(mio::Date start_date, const int num_days,
                                                         const fs::path& data_dir)
{

    // global parameters
    const int num_groups = 6;
    mio::osecir::Parameters params(num_groups);
    params.get<mio::osecir::StartDay>() = mio::get_day_in_year(start_date);
    BOOST_OUTCOME_TRY(set_covid_parameters(params));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));

    auto population_data_path =
        mio::path_join((data_dir / "pydata" / "Germany").string(), "county_current_population.json");

    BOOST_OUTCOME_TRY(auto&& node_ids, mio::get_node_ids(population_data_path, true));
    std::vector<mio::osecir::Model> nodes(node_ids.size(), mio::osecir::Model(num_groups));
    for (auto& node : nodes) {
        node.parameters = params;
    }
    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 2.8);
    auto scaling_factor_icu      = std::vector<double>(size_t(params.get_num_groups()), 0.5);
    scaling_factor_icu[4]        = 25.0;
    scaling_factor_icu[5]        = 5.0;

    const auto& read_function_nodes = mio::osecir::read_input_data_county<mio::osecir::Model>;
    BOOST_OUTCOME_TRY(read_function_nodes(nodes, start_date, node_ids, scaling_factor_infected, scaling_factor_icu,
                                          data_dir.string(), num_days, false));
    return mio::success(nodes);
}

/**
 * Calls the function to generate the extrapolated data based on the reported data. The generated data is saved in the
    * data directory.
 * @param start_date start date of the simulation.
 * @param num_days number of days to extrapolate the data.
 * @param data_dir data directory.
 * @returns created graph or any io errors that happen during reading of the files.
 */
mio::IOResult<void> generate_extrapolated_data(std::vector<mio::osecir::Model> nodes, std::vector<int> node_ids,
                                               mio::Date start_date, const int num_days, const fs::path& data_dir)
{

    auto scaling_factor_infected = std::vector<double>(size_t(nodes[0].parameters.get_num_groups()), 1.0);
    auto scaling_factor_icu      = 1.0;

    BOOST_OUTCOME_TRY(export_input_data_county_timeseries(
        nodes, data_dir.string(), node_ids, start_date, scaling_factor_infected, scaling_factor_icu, num_days,
        mio::path_join(data_dir.string(), "pydata/Germany", "county_divi_ma7.json"),
        mio::path_join(data_dir.string(), "pydata/Germany", "cases_all_county_age_repdate_ma7.json"),
        mio::path_join(data_dir.string(), "pydata/Germany", "county_current_population.json")));

    return mio::success();
}

/*
int main(int argc, char** argv)
{

    mio::set_log_level(mio::LogLevel::warn);
    mio::mpi::init();

    std::string data_dir;
    mio::Date start_date = mio::Date(2023, 6, 1);
    double num_days_sim  = 90.0;

    if (argc == 1) {
        data_dir = "../../data";
    }
    else if (argc == 2) {
        data_dir = argv[1];
    }
    else if (argc == 6) {
        data_dir     = argv[1];
        start_date   = mio::Date(std::atoi(argv[2]), std::atoi(argv[3]), std::atoi(argv[4]));
        num_days_sim = std::atoi(argv[5]);
    }
    else {
        mio::mpi::finalize();
        return 0;
    }

    //mio::thread_local_rng().seed(
    //    {114381446, 2427727386, 806223567, 832414962, 4121923627, 1581162203}); //set seeds, e.g., for debugging
    mio::thread_local_rng().synchronize();
    if (mio::mpi::is_root()) {
        printf("Seeds: ");
        for (auto s : mio::thread_local_rng().get_seeds()) {
            printf("%u, ", s);
        }
        printf("\n");
    }

    //create or load graph
    auto create_extrapolated_data = get_graph(start_date, num_days_sim, data_dir);

    mio::mpi::finalize();
    return 0;
}
*/