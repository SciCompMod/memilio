/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Rene Schmieding
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
#ifndef SECIR_AGERES_SETUPS_H_
#define SECIR_AGERES_SETUPS_H_

#include "memilio/compartments/simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/parameter_studies.h"
#include "models/ode_secir/model.h"
#include "models/ode_secir/parameter_space.h"

namespace mio
{
namespace benchmark
{
namespace detail
{
/**
         * @brief Helper function to create a secir model with consistent setup for use in benchmarking.
         */
mio::osecir::Model make_model(int num)
{

    double cont_freq = 10; // see Polymod study

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model model(num);
    auto nb_groups = model.parameters.get_num_groups();
    double fact    = 1.0 / (double)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::ICUCapacity>(std::numeric_limits<double>::max());
    params.set<mio::osecir::StartDay>(0);
    params.set<mio::osecir::Seasonality>(0);

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::IncubationTime>()[i]       = 5.2;
        params.get<mio::osecir::TimeInfectedSymptoms>()[i] = 6.;
        params.get<mio::osecir::SerialInterval>()[i]       = 4.2;
        params.get<mio::osecir::TimeInfectedSevere>()[i]   = 12;
        params.get<mio::osecir::TimeInfectedCritical>()[i] = 8;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * nb_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * nb_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * nb_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * nb_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * nb_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * nb_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * nb_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact>()[i] = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[i]   = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[i]   = 0.25;
        params.get<mio::osecir::SeverePerInfectedSymptoms>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical>()[i]                = 0.3;
    }

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                               mio::SimulationTime(30.));

    model.apply_constraints();

    return model;
}
} // namespace detail

namespace model
{
/**
         * @brief Secir model with consistent setup for use in benchmarking.
         */
mio::osecir::Model SecirAgeres(size_t num_agegroups)
{
    mio::osecir::Model model = mio::benchmark::detail::make_model(num_agegroups);

    auto nb_groups   = model.parameters.get_num_groups();
    double cont_freq = 10, fact = 1.0 / (double)(size_t)nb_groups;
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                               mio::SimulationTime(30.));

    return model;
}
/**
         * @brief Secir model with consistent setup for use in benchmarking with added dampings.
         */
mio::osecir::Model SecirAgeresDampings(size_t num_agegroups)
{
    mio::osecir::Model model = mio::benchmark::detail::make_model(num_agegroups);

    auto nb_groups   = model.parameters.get_num_groups();
    double cont_freq = 10, fact = 1.0 / (double)(size_t)nb_groups;
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                               mio::SimulationTime(25.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.3),
                               mio::SimulationTime(40.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.8),
                               mio::SimulationTime(60.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.5),
                               mio::SimulationTime(75.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 1.0),
                               mio::SimulationTime(95.));

    return model;
}
/**
         * @brief Secir model with consistent setup for use in benchmarking with added dampings.
         * Dampings are set up to challenge the integrator, not to be realistic.
         */
mio::osecir::Model SecirAgeresAbsurdDampings(size_t num_agegroups)
{
    mio::osecir::Model model = mio::benchmark::detail::make_model(num_agegroups);

    auto nb_groups   = model.parameters.get_num_groups();
    double cont_freq = 10, fact = 1.0 / (double)(size_t)nb_groups;
    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.8),
                               mio::SimulationTime(10.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.5),
                               mio::SimulationTime(11.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.2),
                               mio::SimulationTime(12.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.1),
                               mio::SimulationTime(13.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.9),
                               mio::SimulationTime(30.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.2),
                               mio::SimulationTime(30.5));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                               mio::SimulationTime(31.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.2),
                               mio::SimulationTime(31.5));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.8),
                               mio::SimulationTime(32.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.1),
                               mio::SimulationTime(40.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.001),
                               mio::SimulationTime(44.));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.9),
                               mio::SimulationTime(46.));

    return model;
}
} // namespace model
} // namespace benchmark

} // namespace mio

#endif
