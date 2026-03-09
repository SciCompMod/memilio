/*
* Copyright (C) 2020-2026 MEmilio
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
namespace details
{
/**
 * @brief Helper function to create a secir model with consistent setup for use in benchmarking.
 */
mio::osecir::Model<ScalarType> make_model(int num)
{

    ScalarType cont_freq = 10; // see Polymod study

    ScalarType nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
               nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model<ScalarType> model(num);
    auto nb_groups  = model.parameters.get_num_groups();
    ScalarType fact = 1.0 / (ScalarType)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<mio::osecir::ICUCapacity<ScalarType>>(std::numeric_limits<ScalarType>::max());
    params.set<mio::osecir::StartDay<ScalarType>>(0);
    params.set<mio::osecir::Seasonality<ScalarType>>(0);

    for (auto i = mio::AgeGroup(0); i < nb_groups; i++) {
        params.get<mio::osecir::TimeExposed<ScalarType>>()[i]            = 3.2;
        params.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[i] = 2.0;
        params.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[i]   = 6.0;
        params.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[i]     = 12;
        params.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[i]   = 8;

        model.populations[{i, mio::osecir::InfectionState::Exposed}]            = fact * nb_exp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedNoSymptoms}] = fact * nb_car_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSymptoms}]   = fact * nb_inf_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedSevere}]     = fact * nb_hosp_t0;
        model.populations[{i, mio::osecir::InfectionState::InfectedCritical}]   = fact * nb_icu_t0;
        model.populations[{i, mio::osecir::InfectionState::Recovered}]          = fact * nb_rec_t0;
        model.populations[{i, mio::osecir::InfectionState::Dead}]               = fact * nb_dead_t0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::osecir::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        params.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[i] = 0.05;
        params.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[i]   = 0.67;
        params.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[i]   = 0.09;
        params.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[i]   = 0.25;
        params.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[i]        = 0.2;
        params.get<mio::osecir::CriticalPerSevere<ScalarType>>()[i]                = 0.25;
        params.get<mio::osecir::DeathsPerCritical<ScalarType>>()[i]                = 0.3;
    }

    mio::ContactMatrixGroup<ScalarType>& contact_matrix = params.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0]                                   = mio::ContactMatrix<ScalarType>(
        Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                               mio::SimulationTime<ScalarType>(30.0));

    model.apply_constraints();

    return model;
}
} // namespace details

namespace model
{
/**
 * @brief Secir model with consistent setup for use in benchmarking.
 */
mio::osecir::Model<ScalarType> SecirAgeres(size_t num_agegroups)
{
    mio::osecir::Model<ScalarType> model = mio::benchmark::details::make_model(num_agegroups);

    auto nb_groups       = model.parameters.get_num_groups();
    ScalarType cont_freq = 10, fact = 1.0 / (ScalarType)(size_t)nb_groups;
    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(
        Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                               mio::SimulationTime<ScalarType>(30.0));

    return model;
}
/**
 * @brief Secir model with consistent setup for use in benchmarking with added dampings.
 */
mio::osecir::Model<ScalarType> SecirAgeresDampings(size_t num_agegroups)
{
    mio::osecir::Model<ScalarType> model = mio::benchmark::details::make_model(num_agegroups);

    auto nb_groups       = model.parameters.get_num_groups();
    ScalarType cont_freq = 10, fact = 1.0 / (ScalarType)(size_t)nb_groups;
    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(
        Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                               mio::SimulationTime<ScalarType>(25.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.3),
                               mio::SimulationTime<ScalarType>(40.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.8),
                               mio::SimulationTime<ScalarType>(60.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.5),
                               mio::SimulationTime<ScalarType>(75.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 1.0),
                               mio::SimulationTime<ScalarType>(95.0));

    return model;
}
/**
 * @brief Secir model with consistent setup for use in benchmarking with added dampings.
 * Dampings are set up to challenge the integrator, not to be realistic.
 */
mio::osecir::Model<ScalarType> SecirAgeresAbsurdDampings(size_t num_agegroups)
{
    mio::osecir::Model<ScalarType> model = mio::benchmark::details::make_model(num_agegroups);

    auto nb_groups       = model.parameters.get_num_groups();
    ScalarType cont_freq = 10, fact = 1.0 / (ScalarType)(size_t)nb_groups;
    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] = mio::ContactMatrix<ScalarType>(
        Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.8),
                               mio::SimulationTime<ScalarType>(10.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.5),
                               mio::SimulationTime<ScalarType>(11.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.2),
                               mio::SimulationTime<ScalarType>(12.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.1),
                               mio::SimulationTime<ScalarType>(13.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.9),
                               mio::SimulationTime<ScalarType>(30.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.2),
                               mio::SimulationTime<ScalarType>(30.5));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                               mio::SimulationTime<ScalarType>(31.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.2),
                               mio::SimulationTime<ScalarType>(31.5));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.8),
                               mio::SimulationTime<ScalarType>(32.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.1),
                               mio::SimulationTime<ScalarType>(40.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.001),
                               mio::SimulationTime<ScalarType>(44.0));
    contact_matrix.add_damping(Eigen::MatrixX<ScalarType>::Constant((size_t)nb_groups, (size_t)nb_groups, 0.9),
                               mio::SimulationTime<ScalarType>(46.0));

    return model;
}
} // namespace model
} // namespace benchmark

} // namespace mio

#endif
