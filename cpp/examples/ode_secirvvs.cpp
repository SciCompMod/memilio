/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Martin J. Kuehn
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
#include "ode_secirvvs/model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 30;
    double dt   = 0.1;

    mio::log_info("Simulating SECIRVVS; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::osecirvvs::Model<double> model(1);

    for (mio::AgeGroup i = 0; i < model.parameters.get_num_groups(); i++) {
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedNaive}]                                = 10;
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedImprovedImmunity}]                     = 11;
        model.populations[{i, mio::osecirvvs::InfectionState::ExposedPartialImmunity}]                      = 12;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaive}]                     = 13;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsNaiveConfirmed}]            = 13;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunity}]           = 14;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsPartialImmunityConfirmed}]  = 14;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunity}]          = 15;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedNoSymptomsImprovedImmunityConfirmed}] = 15;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive}]                       = 5;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsNaiveConfirmed}]              = 5;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity}]             = 6;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunityConfirmed}]    = 6;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity}]            = 7;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunityConfirmed}]   = 7;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereNaive}]                         = 8;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSevereImprovedImmunity}]              = 1;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedSeverePartialImmunity}]               = 2;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalNaive}]                       = 3;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalPartialImmunity}]             = 4;
        model.populations[{i, mio::osecirvvs::InfectionState::InfectedCriticalImprovedImmunity}]            = 5;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptibleImprovedImmunity}]                 = 6;
        model.populations[{i, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}]                  = 7;
        model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::DeadNaive}]                    = 0;
        model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::DeadPartialImmunity}]          = 0;
        model.populations[{(mio::AgeGroup)0, mio::osecirvvs::InfectionState::DeadImprovedImmunity}]         = 0;
        model.populations.set_difference_from_group_total<mio::AgeGroup>(
            {i, mio::osecirvvs::InfectionState::SusceptibleNaive}, 1000);
    }

    model.parameters.get<mio::osecirvvs::ICUCapacity<double>>()          = 100;
    model.parameters.get<mio::osecirvvs::TestAndTraceCapacity<double>>() = 0.0143;
    const size_t daily_vaccinations                                      = 10;
    model.parameters.get<mio::osecirvvs::DailyFirstVaccination<double>>().resize(mio::SimulationDay((size_t)tmax + 1));
    model.parameters.get<mio::osecirvvs::DailyFullVaccination<double>>().resize(mio::SimulationDay((size_t)tmax + 1));
    for (size_t i = 0; i < tmax + 1; ++i) {
        auto num_vaccinations = static_cast<double>(i * daily_vaccinations);
        model.parameters
            .get<mio::osecirvvs::DailyFirstVaccination<double>>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
        model.parameters
            .get<mio::osecirvvs::DailyFullVaccination<double>>()[{(mio::AgeGroup)0, mio::SimulationDay(i)}] =
            num_vaccinations;
    }
    model.parameters.get<mio::osecirvvs::DynamicNPIsImplementationDelay<double>>() = 7;

    auto& contacts       = model.parameters.get<mio::osecirvvs::ContactPatterns<double>>();
    auto& contact_matrix = contacts.get_cont_freq_mat();
    contact_matrix[0].get_baseline().setConstant(0.5);
    contact_matrix[0].get_baseline().diagonal().setConstant(5.0);
    contact_matrix[0].add_damping(0.3, mio::SimulationTime(5.0));

    //times
    model.parameters.get<mio::osecirvvs::TimeExposed<double>>()[mio::AgeGroup(0)]            = 3.33;
    model.parameters.get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>()[mio::AgeGroup(0)] = 1.87;
    model.parameters.get<mio::osecirvvs::TimeInfectedSymptoms<double>>()[mio::AgeGroup(0)]   = 7;
    model.parameters.get<mio::osecirvvs::TimeInfectedSevere<double>>()[mio::AgeGroup(0)]     = 6;
    model.parameters.get<mio::osecirvvs::TimeInfectedCritical<double>>()[mio::AgeGroup(0)]   = 7;

    //probabilities
    model.parameters.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)] = 0.15;
    model.parameters.get<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>()[mio::AgeGroup(0)]   = 0.5;
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    model.parameters.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>()[mio::AgeGroup(0)]    = 0.0;
    model.parameters.get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<double>>()[mio::AgeGroup(0)] = 0.4;
    model.parameters.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>()[mio::AgeGroup(0)]    = 0.2;
    model.parameters.get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>()[mio::AgeGroup(0)]         = 0.1;
    model.parameters.get<mio::osecirvvs::CriticalPerSevere<double>>()[mio::AgeGroup(0)]                 = 0.1;
    model.parameters.get<mio::osecirvvs::DeathsPerCritical<double>>()[mio::AgeGroup(0)]                 = 0.1;

    model.parameters.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>()[mio::AgeGroup(0)]           = 0.8;
    model.parameters.get<mio::osecirvvs::ReducExposedImprovedImmunity<double>>()[mio::AgeGroup(0)]          = 0.331;
    model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>()[mio::AgeGroup(0)]  = 0.65;
    model.parameters.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>()[mio::AgeGroup(0)] = 0.243;
    model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[mio::AgeGroup(0)] =
        0.1;
    model.parameters.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[mio::AgeGroup(0)] =
        0.091;
    model.parameters.get<mio::osecirvvs::ReducTimeInfectedMild<double>>()[mio::AgeGroup(0)] = 0.9;

    model.parameters.get<mio::osecirvvs::Seasonality<double>>() = 0.2;

    model.apply_constraints();

    // use adaptive Runge-Kutta-Fehlberg45 scheme as integrator
    // auto integrator = std::make_shared<mio::RKIntegratorCore>();
    // integrator->set_dt_min(0.3);
    // integrator->set_dt_max(1.0);
    // integrator->set_rel_tolerance(1e-4);
    // integrator->set_abs_tolerance(1e-1);
    // mio::TimeSeries<double> result = mio::osecirvvs::simulate(t0, tmax, dt, model, integrator);

    // use default Cash-Karp adaptive integrator
    mio::TimeSeries<double> result = mio::osecirvvs::simulate<double>(t0, tmax, dt, model);

    bool print_to_terminal = true;

    if (print_to_terminal) {
        printf("\n%.14f ", result.get_last_time());
        for (size_t j = 0; j < (size_t)mio::osecirvvs::InfectionState::Count; j++) {
            printf("compartment %d: %.14f\n", (int)j, result.get_last_value()[j]);
        }
    }
}
