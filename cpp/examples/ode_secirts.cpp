/* 
* Copyright (C) 2020-2024 MEmilio
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
#include "ode_secirvvs/analyze_result.h"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 100;
    double dt   = 0.1;

    mio::log_info("Simulating SECIRVVS; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::osecirvvs::Model model(3);
    auto nb_groups = model.parameters.get_num_groups();

    for (mio::AgeGroup i = 0; i < nb_groups; i++) {
        // population
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedNaive}]                                = 20;
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedImprovedImmunity}]                     = 20;
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedPartialImmunity}]                      = 20;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}]                     = 30;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}]            = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}]           = 30;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]  = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}]          = 30;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]                       = 40;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]              = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}]             = 40;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]    = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]            = 40;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]   = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereNaive}]                         = 30;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereImprovedImmunity}]              = 30;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSeverePartialImmunity}]               = 30;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalNaive}]                       = 20;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalPartialImmunity}]             = 20;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalImprovedImmunity}]            = 20;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptibleNaive}]                            = 1000;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}]                  = 1200;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}]                 = 1000;
        model.populations[{i, mio::osecirvvs::InfectionState::TemporaryImmunPartialImmunity}]               = 60;
        model.populations[{i, mio::osecirvvs::InfectionState::TemporaryImmunImprovedImmunity}]              = 70;
        model.populations[{i, mio::osecirvvs::InfectionState::DeadNaive}]                                   = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::DeadPartialImmunity}]                         = 0;
        model.populations[{i, mio::osecirvvs::InfectionState::DeadImprovedImmunity}]                        = 0;

        // parameters
        //times
        model.parameters.get<mio::osecirvvs::TimeExposed>()[i]                = 3.33;
        model.parameters.get<mio::osecirvvs::TimeInfectedNoSymptoms>()[i]     = 1.87;
        model.parameters.get<mio::osecirvvs::TimeInfectedSymptoms>()[i]       = 7;
        model.parameters.get<mio::osecirvvs::TimeInfectedSevere>()[i]         = 6;
        model.parameters.get<mio::osecirvvs::TimeInfectedCritical>()[i]       = 7;
        model.parameters.get<mio::osecirvvs::TimeTemporaryImmunityPI>()[i]    = 60;
        model.parameters.get<mio::osecirvvs::TimeTemporaryImmunityPI>()[i]    = 60;
        model.parameters.get<mio::osecirvvs::TimeTemporaryImmunityII>()[i]    = 60;
        model.parameters.get<mio::osecirvvs::TimeWaningPartialImmunity>()[i]  = 180;
        model.parameters.get<mio::osecirvvs::TimeWaningImprovedImmunity>()[i] = 180;

        //probabilities
        model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact>()[i]  = 0.15;
        model.parameters.get<mio::osecirvvs::RelativeTransmissionNoSymptoms>()[i]    = 0.5;
        model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic>()[i]    = 0.0;
        model.parameters.get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic>()[i] = 0.4;
        model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms>()[i]    = 0.2;
        model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms>()[i]         = 0.1;
        model.parameters.get<mio::osecirvvs::CriticalPerSevere>()[i]                 = 0.1;
        model.parameters.get<mio::osecirvvs::DeathsPerCritical>()[i]                 = 0.1;

        model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity>()[i]                     = 0.8;
        model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity>()[i]                    = 0.331;
        model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity>()[i]            = 0.65;
        model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity>()[i]           = 0.243;
        model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity>()[i]  = 0.1;
        model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity>()[i] = 0.091;
        model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild>()[i]                           = 0.9;
    }

    model.parameters.get<mio::osecirvvs::ICUCapacity>()          = 100;
    model.parameters.get<mio::osecirvvs::TestAndTraceCapacity>() = 0.0143;
    const size_t daily_vaccinations                              = 10;
    const size_t num_days                                        = 300;
    model.parameters.get<mio::osecirvvs::DailyPartialVaccination>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirvvs::DailyFullVaccination>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirvvs::DailyBoosterVaccination>().resize(mio::SimulationDay(num_days));
    for (size_t i = 0; i < num_days; ++i) {
        auto num_vaccinations = static_cast<double>(i * daily_vaccinations);
        model.parameters.get<mio::osecirvvs::DailyPartialVaccination>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
        model.parameters.get<mio::osecirvvs::DailyFullVaccination>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
        model.parameters.get<mio::osecirvvs::DailyBoosterVaccination>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
    }

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecirvvs::ContactPatterns>();
    const double cont_freq                  = 10;
    const double fact                       = 1.0 / (double)(size_t)nb_groups;
    contact_matrix[0] =
        mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));
    contact_matrix.add_damping(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, 0.7),
                               mio::SimulationTime(30.));

    model.parameters.get<mio::osecirvvs::Seasonality>() = 0.2;

    model.apply_constraints();

    mio::TimeSeries<double> result = simulate(t0, tmax, dt, model);

    bool print_to_terminal = true;

    if (print_to_terminal) {
        auto result_interpolated = mio::interpolate_simulation_result(result);
        for (auto t_indx = 0; t_indx < result_interpolated.get_num_time_points(); t_indx++) {
            double timm_pi = 0.0;
            double timm_ii = 0.0;
            for (mio::AgeGroup i = 0; i < nb_groups; i++) {
                timm_pi += result_interpolated.get_value(t_indx)[model.populations.get_flat_index(
                    {i, mio::osecirvvs::InfectionState::TemporaryImmunPartialImmunity})];
                timm_ii += result_interpolated.get_value(t_indx)[model.populations.get_flat_index(
                    {i, mio::osecirvvs::InfectionState::TemporaryImmunImprovedImmunity})];
            }
            printf("t=%i, timm_pi=%f, timm_ii=%f\n", int(result_interpolated.get_time(t_indx)), timm_pi, timm_ii);
        }
    }
}
