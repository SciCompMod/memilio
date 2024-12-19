/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "ode_seir/model.h"
#include "ode_seir/infection_state.h"
#include "ode_seir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/math/euler.h"
#include "memilio/utils/time_series.h"

#include "memilio/utils/time_series.h"
#include "memilio/io/io.h"
#include "memilio/io/result_io.h"
#include "memilio/io/epi_data.h"

/**
 * Set epidemiological parameters of Sars-CoV-2 for a immune-naive
 * population and wild type variant.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_covid_parameters(mio::oseir::Parameters<double>& params, bool synthetic_population)
{
    params.template set<mio::oseir::TimeExposed<>>(3.335);
    if (!synthetic_population) {
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(0)] = 8.0096875;
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(1)] = 8.0096875;
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(2)] = 8.2182;
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(3)] = 8.1158;
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(4)] = 8.033;
        params.get<mio::oseir::TimeInfected<>>()[mio::AgeGroup(5)] = 7.985;

        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(0)] = 0.03;
        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(1)] = 0.06;
        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(2)] = 0.06;
        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(3)] = 0.06;
        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(4)] = 0.09;
        params.get<mio::oseir::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(5)] = 0.175;
    }
    else {
        params.template set<mio::oseir::TimeInfected<>>(8.097612257);

        params.template set<mio::oseir::TransmissionProbabilityOnContact<>>(0.07333);
    }

    printf("Setting epidemiological parameters successful.\n");
    return mio::success();
}

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
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::oseir::Parameters<double>& params,
                                         bool synthetic_population)
{
    if (!synthetic_population) {
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
        params.get<mio::oseir::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);
    }
    else {
        mio::ContactMatrixGroup& contact_matrix = params.get<mio::oseir::ContactPatterns<>>().get_cont_freq_mat();
        contact_matrix[0].get_baseline().setConstant(7.95 / (size_t)params.get_num_groups());
    }

    printf("Setting contact matrices successful.\n");
    return mio::success();
}

template <typename FP = ScalarType>
mio::IOResult<void> set_population_data(mio::oseir::Model<FP>& model, const fs::path& data_dir)
{
    BOOST_OUTCOME_TRY(
        auto&& node_ids,
        mio::get_node_ids((data_dir / "pydata" / "Germany" / "county_current_population_nrw.json").string(), true,
                          true));

    BOOST_OUTCOME_TRY(const auto&& population_data,
                      mio::read_population_data(
                          (data_dir / "pydata" / "Germany" / "county_current_population_nrw.json").string(), true));

    for (auto&& entry : population_data) {
        auto it = std::find_if(node_ids.begin(), node_ids.end(), [&entry](auto r) {
            return r == 0 ||
                   (entry.county_id && mio::regions::StateId(r) == mio::regions::get_state_id(int(*entry.county_id))) ||
                   (entry.county_id && mio::regions::CountyId(r) == *entry.county_id) ||
                   (entry.district_id && mio::regions::DistrictId(r) == *entry.district_id);
        });
        if (it != node_ids.end()) {
            for (size_t age = 0; age < (size_t)model.parameters.get_num_groups(); age++) {
                model.populations[{mio::AgeGroup(age), mio::oseir::InfectionState::Susceptible}] +=
                    entry.population[mio::AgeGroup(age)];
            }
        }
    }

    printf("Setting population data successful.\n");
    return mio::success();
}
template <typename FP = ScalarType>
mio::IOResult<void> set_parameters_and_population(mio::oseir::Model<FP>& model, const fs::path& data_dir,
                                                  bool synthetic_population)
{
    auto& populations = model.populations;
    auto& parameters  = model.parameters;

    size_t number_age_groups = (size_t)parameters.get_num_groups();

    if (synthetic_population) {
        printf("Data is not compatible, using synthetic population instead.\n");
        for (size_t j = 0; j < number_age_groups; j++) {
            model.populations[{mio::AgeGroup(j), mio::oseir::InfectionState::Exposed}]     = 100;
            model.populations[{mio::AgeGroup(j), mio::oseir::InfectionState::Susceptible}] = 999900;
        }
    }
    else {
        BOOST_OUTCOME_TRY(set_population_data(model, data_dir));
        populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] -= 100;
        populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}] += 100;
    }

    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, parameters, synthetic_population))

    BOOST_OUTCOME_TRY(set_covid_parameters(parameters, synthetic_population));

    return mio::success();
}

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    ScalarType t0   = 0.;
    ScalarType tmax = 50.;
    ScalarType dt   = 0.1;

    ScalarType number_age_groups = 6;
    bool synthetic_population    = false;
    if (number_age_groups != 6) {
        synthetic_population = true;
    }

    mio::log_info("Simulating ODE SEIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    const std::string& data_dir = "";

    mio::oseir::Model<ScalarType> model(number_age_groups);
    auto result_prepare_simulation = set_parameters_and_population(model, data_dir, synthetic_population);

    std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = std::make_shared<mio::EulerIntegratorCore<>>();

    auto seir = simulate(t0, tmax, dt, model, integrator);

    auto reproduction_numbers = model.get_reproduction_numbers(seir);
    std::cout << "\nbasis reproduction number: " << reproduction_numbers[0] << "\n";

    // seir.print_table({"S", "E", "I", "R"});
    // std::cout << "\nnumber total: " << seir.get_last_value().sum() << "\n";
}
