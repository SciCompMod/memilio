/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Carlotta Gerstein
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

#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/custom_index_array.h"
#include "models/ode_metapop/infection_state.h"
#include "models/ode_metapop/model.h"
#include "models/ode_metapop/parameters.h"
#include "models/ode_metapop/regions.h"

#include <omp.h>

bool age_groups = false;

template <typename FP>
void set_contact_matrix(mio::oseirmetapop::Model<FP>& model)
{
    if (age_groups) {
        Eigen::MatrixXd contact_matrix_eigen(6, 6);
        contact_matrix_eigen << 3.9547, 1.1002, 2.9472, 2.05, 0.3733, 0.0445, 0.3327, 3.5892, 1.236, 1.9208, 0.2681,
            0.0161, 0.246, 0.7124, 5.6518, 3.2939, 0.2043, 0.0109, 0.1742, 0.8897, 3.3124, 4.5406, 0.4262, 0.0214,
            0.0458, 0.1939, 0.5782, 1.3825, 1.473, 0.0704, 0.1083, 0.1448, 0.4728, 0.9767, 0.6266, 0.1724;
        mio::ContactMatrixGroup& contact_matrix =
            model.parameters.template get<mio::oseirmetapop::ContactPatterns<>>().get_cont_freq_mat();
        contact_matrix[0].get_baseline() = contact_matrix_eigen;
    }
    {
        mio::ContactMatrixGroup& contact_matrix =
            model.parameters.template get<mio::oseirmetapop::ContactPatterns<>>().get_cont_freq_mat();
        contact_matrix[0].get_baseline().setConstant(7.95);
    }
}

/**
 * Set epidemiological parameters of Sars-CoV-2 for a immune-naive
 * population and wild type variant.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
void set_covid_parameters(mio::oseirmetapop::Parameters<double>& params)
{
    params.template set<mio::oseirmetapop::TimeExposed<>>(3.335);

    if (age_groups) {
        params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(0)] = 8.0096875;
        params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(1)] = 8.0096875;
        params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(2)] = 8.2182;
        params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(3)] = 8.1158;
        params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(4)] = 8.033;
        params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(5)] = 7.985;

        params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(0)] = 0.03;
        params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(1)] = 0.06;
        params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(2)] = 0.06;
        params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(3)] = 0.06;
        params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(4)] = 0.09;
        params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(5)] = 0.175;
    }
    else {
        params.get<mio::oseirmetapop::TransmissionProbabilityOnContact<>>()[mio::AgeGroup(0)] = 0.07333;
        params.get<mio::oseirmetapop::TimeInfected<>>()[mio::AgeGroup(0)]                     = 8.097612257;
    }
}

template <typename FP = ScalarType>
void set_mobility_weights(mio::oseirmetapop::Model<FP>& model)
{
    size_t number_regions    = (size_t)model.parameters.get_num_regions();
    double fraction_commuter = 1. / (2 * number_regions);
    Eigen::MatrixXd mobility_data_commuter =
        Eigen::MatrixXd::Constant(number_regions, number_regions, fraction_commuter) -
        fraction_commuter *
            Eigen::MatrixXd::Identity(number_regions, number_regions); // Ensure that the diagonal is zero
    for (size_t county_idx_i = 0; county_idx_i < number_regions; ++county_idx_i) {
        mobility_data_commuter(county_idx_i, county_idx_i) = 1 - mobility_data_commuter.rowwise().sum()(county_idx_i);
    }
    model.parameters.template get<mio::oseirmetapop::CommutingStrengths<>>().get_cont_freq_mat()[0].get_baseline() =
        mobility_data_commuter;
}

template <typename FP = ScalarType>
void set_parameters_and_population(mio::oseirmetapop::Model<FP>& model)
{
    auto& populations = model.populations;
    auto& parameters  = model.parameters;

    size_t number_regions    = (size_t)parameters.get_num_regions();
    size_t number_age_groups = (size_t)parameters.get_num_agegroups();
    for (size_t j = 0; j < number_age_groups; j++) {
        for (size_t i = 0; i < number_regions; i++) {
            populations[{mio::oseirmetapop::Region(i), mio::AgeGroup(j),
                         mio::oseirmetapop::InfectionState::Susceptible}] = 60000;
        }
    }
    populations[{mio::oseirmetapop::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Exposed}] += 100;
    populations[{mio::oseirmetapop::Region(0), mio::AgeGroup(0), mio::oseirmetapop::InfectionState::Susceptible}] -=
        100;
    set_mobility_weights(model);

    set_contact_matrix(model);

    set_covid_parameters(parameters);

    mio::ContactMatrixGroup& commuting_strengths =
        parameters.template get<mio::oseirmetapop::CommutingStrengths<>>().get_cont_freq_mat();

    auto& population_after_commuting = model.m_population_after_commuting;
    for (size_t region_n = 0; region_n < number_regions; ++region_n) {
        for (size_t age = 0; age < number_age_groups; ++age) {
            double population_n = 0;
            for (size_t state = 0; state < (size_t)mio::oseirmetapop::InfectionState::Count; state++) {
                population_n += populations[{mio::oseirmetapop::Region(region_n), mio::AgeGroup(age),
                                             mio::oseirmetapop::InfectionState(state)}];
            }
            population_after_commuting[{mio::oseirmetapop::Region(region_n), mio::AgeGroup(age)}] += population_n;
            for (size_t region_m = 0; region_m < number_regions; ++region_m) {
                population_after_commuting[{mio::oseirmetapop::Region(region_n), mio::AgeGroup(age)}] -=
                    commuting_strengths[0].get_baseline()(region_n, region_m) * population_n;
                population_after_commuting[{mio::oseirmetapop::Region(region_m), mio::AgeGroup(age)}] +=
                    commuting_strengths[0].get_baseline()(region_n, region_m) * population_n;
            }
        }
    }
}

void simulate(size_t num_warm_up_runs, size_t num_runs, size_t number_regions, ScalarType tmax)
{
    mio::set_log_level(mio::LogLevel::off);
    ScalarType t0                = 0.;
    ScalarType dt                = 0.1;
    ScalarType number_age_groups = 1;
    if (age_groups) {
        number_age_groups = 6;
    }

    mio::oseirmetapop::Model<ScalarType> model(number_regions, number_age_groups);
    set_parameters_and_population(model);

    std::shared_ptr<mio::IntegratorCore<ScalarType>> integrator = std::make_shared<mio::EulerIntegratorCore<>>();

    std::cout << "{ \"Regions\": " << number_regions << ", " << std::endl;

    // Warm up runs.
    for (size_t i = 0; i < num_warm_up_runs; i++) {
        simulate(t0, tmax, dt, model, integrator);
    }
    auto result = simulate(t0, tmax, dt, model, integrator);

    // Runs with timing.
    ScalarType total = 0;
    for (size_t i = 0; i < num_runs; i++) {
        double runtime = simulate_runtime(t0, tmax, dt, model, integrator);
        total += runtime;
    }
    std::cout << "\"Time\": " << total / num_runs << "\n}," << std::endl;
}

int main(int argc, char** argv)
{
    const ScalarType tmax = 20;
    size_t warm_up        = 10;
    size_t num_runs       = 100;
    size_t num_regions    = 10;
    if (argc > 3) {
        warm_up     = std::stod(argv[1]);
        num_runs    = std::stod(argv[2]);
        num_regions = std::stod(argv[3]);
    }
    simulate(warm_up, num_runs, num_regions, tmax);
    return 0;
}
