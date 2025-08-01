/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Henrik Zunker
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
#include "ode_secirts/analyze_result.h"
#include "ode_secirts/model.h"
#include "ode_secirts/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/data/analyze_result.h"

int main()
{
    // This example demonstrates how to simulate a SECIRTS model.
    // The SECIRTS model is an extension of the SECIRVVS model that includes waning and  temporary immunity.
    // After the simulation, the aggregated size of the temporary immunity states are printed.
    mio::set_log_level(mio::LogLevel::warn);

    double t0   = 0;
    double tmax = 30;
    double dt   = 0.1;

    mio::log_info("Simulating SECIRTS; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::osecirts::Model<double> model(3);
    auto nb_groups = model.parameters.get_num_groups();

    for (mio::AgeGroup i = 0; i < nb_groups; i++) {
        // population
        model.populations[{i, mio::osecirts::InfectionState::ExposedNaive}]                                = 20;
        model.populations[{i, mio::osecirts::InfectionState::ExposedImprovedImmunity}]                     = 20;
        model.populations[{i, mio::osecirts::InfectionState::ExposedPartialImmunity}]                      = 20;
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaive}]                     = 30;
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaiveConfirmed}]            = 0;
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunity}]           = 30;
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]  = 0;
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunity}]          = 30;
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 0;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaive}]                       = 40;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaiveConfirmed}]              = 0;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunity}]             = 40;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]    = 0;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunity}]            = 40;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]   = 0;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSevereNaive}]                         = 30;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSevereImprovedImmunity}]              = 30;
        model.populations[{i, mio::osecirts::InfectionState::InfectedSeverePartialImmunity}]               = 30;
        model.populations[{i, mio::osecirts::InfectionState::InfectedCriticalNaive}]                       = 20;
        model.populations[{i, mio::osecirts::InfectionState::InfectedCriticalPartialImmunity}]             = 20;
        model.populations[{i, mio::osecirts::InfectionState::InfectedCriticalImprovedImmunity}]            = 20;
        model.populations[{i, mio::osecirts::InfectionState::SusceptibleNaive}]                            = 1000;
        model.populations[{i, mio::osecirts::InfectionState::SusceptiblePartialImmunity}]                  = 1200;
        model.populations[{i, mio::osecirts::InfectionState::SusceptibleImprovedImmunity}]                 = 1000;
        model.populations[{i, mio::osecirts::InfectionState::TemporaryImmunePartialImmunity}]              = 60;
        model.populations[{i, mio::osecirts::InfectionState::TemporaryImmuneImprovedImmunity}]             = 70;
        model.populations[{i, mio::osecirts::InfectionState::DeadNaive}]                                   = 0;
        model.populations[{i, mio::osecirts::InfectionState::DeadPartialImmunity}]                         = 0;
        model.populations[{i, mio::osecirts::InfectionState::DeadImprovedImmunity}]                        = 0;

        // parameters
        //times
        model.parameters.get<mio::osecirts::TimeExposed<double>>()[i]                = 3.33;
        model.parameters.get<mio::osecirts::TimeInfectedNoSymptoms<double>>()[i]     = 1.87;
        model.parameters.get<mio::osecirts::TimeInfectedSymptoms<double>>()[i]       = 7;
        model.parameters.get<mio::osecirts::TimeInfectedSevere<double>>()[i]         = 6;
        model.parameters.get<mio::osecirts::TimeInfectedCritical<double>>()[i]       = 7;
        model.parameters.get<mio::osecirts::TimeTemporaryImmunityPI<double>>()[i]    = 60;
        model.parameters.get<mio::osecirts::TimeTemporaryImmunityII<double>>()[i]    = 60;
        model.parameters.get<mio::osecirts::TimeWaningPartialImmunity<double>>()[i]  = 180;
        model.parameters.get<mio::osecirts::TimeWaningImprovedImmunity<double>>()[i] = 180;

        //probabilities
        model.parameters.get<mio::osecirts::TransmissionProbabilityOnContact<double>>()[i]  = 0.15;
        model.parameters.get<mio::osecirts::RelativeTransmissionNoSymptoms<double>>()[i]    = 0.5;
        model.parameters.get<mio::osecirts::RiskOfInfectionFromSymptomatic<double>>()[i]    = 0.0;
        model.parameters.get<mio::osecirts::MaxRiskOfInfectionFromSymptomatic<double>>()[i] = 0.4;
        model.parameters.get<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>()[i]    = 0.2;
        model.parameters.get<mio::osecirts::SeverePerInfectedSymptoms<double>>()[i]         = 0.1;
        model.parameters.get<mio::osecirts::CriticalPerSevere<double>>()[i]                 = 0.1;
        model.parameters.get<mio::osecirts::DeathsPerCritical<double>>()[i]                 = 0.1;

        model.parameters.get<mio::osecirts::ReducExposedPartialImmunity<double>>()[i]                     = 0.8;
        model.parameters.get<mio::osecirts::ReducExposedImprovedImmunity<double>>()[i]                    = 0.331;
        model.parameters.get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<double>>()[i]            = 0.65;
        model.parameters.get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<double>>()[i]           = 0.243;
        model.parameters.get<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[i]  = 0.1;
        model.parameters.get<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[i] = 0.091;
        model.parameters.get<mio::osecirts::ReducTimeInfectedMild<double>>()[i]                           = 0.9;
    }

    model.parameters.get<mio::osecirts::ICUCapacity<double>>()          = 100;
    model.parameters.get<mio::osecirts::TestAndTraceCapacity<double>>() = 0.0143;
    const size_t daily_vaccinations                                     = 10;
    const size_t num_days                                               = 300;
    model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirts::DailyFullVaccinations<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>().resize(mio::SimulationDay(num_days));
    for (size_t i = 0; i < num_days; ++i) {
        for (mio::AgeGroup j = 0; j < nb_groups; ++j) {
            auto num_vaccinations = static_cast<double>(i * daily_vaccinations);
            model.parameters.get<mio::osecirts::DailyPartialVaccinations<double>>()[{j, mio::SimulationDay(i)}] =
                num_vaccinations;
            model.parameters.get<mio::osecirts::DailyFullVaccinations<double>>()[{j, mio::SimulationDay(i)}] =
                num_vaccinations;
            model.parameters.get<mio::osecirts::DailyBoosterVaccinations<double>>()[{j, mio::SimulationDay(i)}] =
                num_vaccinations;
        }
    }

    mio::ContactMatrixGroup<double>& contact_matrix = model.parameters.get<mio::osecirts::ContactPatterns<double>>();
    const double cont_freq                          = 10;
    const double fact                               = 1.0 / (double)(size_t)nb_groups;
    contact_matrix[0] =
        mio::ContactMatrix<double>(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                               mio::SimulationTime<double>(30.));

    model.parameters.get<mio::osecirts::Seasonality<double>>() = 0.2;

    model.apply_constraints();

    mio::TimeSeries<double> result = mio::osecirts::simulate<double>(t0, tmax, dt, model);

    bool print_to_terminal = true;

    if (print_to_terminal) {
        auto result_interpolated = mio::interpolate_simulation_result(result);
        for (auto t_indx = 0; t_indx < result_interpolated.get_num_time_points(); t_indx++) {
            double timm_pi = 0.0;
            double timm_ii = 0.0;
            for (mio::AgeGroup i = 0; i < nb_groups; i++) {
                timm_pi += result_interpolated.get_value(t_indx)[model.populations.get_flat_index(
                    {i, mio::osecirts::InfectionState::TemporaryImmunePartialImmunity})];
                timm_ii += result_interpolated.get_value(t_indx)[model.populations.get_flat_index(
                    {i, mio::osecirts::InfectionState::TemporaryImmuneImprovedImmunity})];
            }
            printf("t=%i, timm_pi=%f, timm_ii=%f\n", int(result_interpolated.get_time(t_indx)), timm_pi, timm_ii);
        }
    }
}
