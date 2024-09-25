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
#include "ode_secirts/analyze_result.h"
#include "ode_secirts/model.h"
#include "ode_secirts/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

void initialize_population(mio::osecirts::Model<double>& model, double variation_percent)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(1.0 - variation_percent, 1.0 + variation_percent);

    for (mio::AgeGroup i = 0; i < model.parameters.get_num_groups(); i++) {
        model.populations[{i, mio::osecirts::InfectionState::ExposedNaive}]                = 20 * dis(gen);
        model.populations[{i, mio::osecirts::InfectionState::InfectedNoSymptomsNaive}]     = 30 * dis(gen);
        model.populations[{i, mio::osecirts::InfectionState::InfectedSymptomsNaive}]       = 40 * dis(gen);
        model.populations[{i, mio::osecirts::InfectionState::InfectedSevereNaive}]         = 30 * dis(gen);
        model.populations[{i, mio::osecirts::InfectionState::InfectedCriticalNaive}]       = 20 * dis(gen);
        model.populations[{i, mio::osecirts::InfectionState::SusceptibleNaive}]            = 1000 * dis(gen);
        model.populations[{i, mio::osecirts::InfectionState::SusceptiblePartialImmunity}]  = 1200 * dis(gen);
        model.populations[{i, mio::osecirts::InfectionState::SusceptibleImprovedImmunity}] = 1000 * dis(gen);
        model.populations[{i, mio::osecirts::InfectionState::DeadNaive}]                   = 0 * dis(gen);
    }
}

double run_model_with_params(mio::osecirts::Model<double>& model, const std::vector<double>& params, double t0,
                             double tmax, double dt)
{
    // set population
    initialize_population(model, 0.2);

    // Set parameters
    auto& param_set = model.parameters;
    int i           = 0;

    param_set.get<mio::osecirts::TimeExposed<double>>()[mio::AgeGroup(0)]                = params[i++];
    param_set.get<mio::osecirts::TimeInfectedNoSymptoms<double>>()[mio::AgeGroup(0)]     = params[i++];
    param_set.get<mio::osecirts::TimeInfectedSymptoms<double>>()[mio::AgeGroup(0)]       = params[i++];
    param_set.get<mio::osecirts::TimeInfectedSevere<double>>()[mio::AgeGroup(0)]         = params[i++];
    param_set.get<mio::osecirts::TimeInfectedCritical<double>>()[mio::AgeGroup(0)]       = params[i++];
    param_set.get<mio::osecirts::TimeTemporaryImmunityPI<double>>()[mio::AgeGroup(0)]    = params[i++];
    param_set.get<mio::osecirts::TimeTemporaryImmunityII<double>>()[mio::AgeGroup(0)]    = params[i++];
    param_set.get<mio::osecirts::TimeWaningPartialImmunity<double>>()[mio::AgeGroup(0)]  = params[i++];
    param_set.get<mio::osecirts::TimeWaningImprovedImmunity<double>>()[mio::AgeGroup(0)] = params[i++];

    param_set.get<mio::osecirts::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)] = params[i++];
    param_set.get<mio::osecirts::RelativeTransmissionNoSymptoms<double>>()[mio::AgeGroup(0)]   = params[i++];
    param_set.get<mio::osecirts::RiskOfInfectionFromSymptomatic<double>>()[mio::AgeGroup(0)]   = params[i++];
    param_set.get<mio::osecirts::RecoveredPerInfectedNoSymptoms<double>>()[mio::AgeGroup(0)]   = params[i++];
    param_set.get<mio::osecirts::SeverePerInfectedSymptoms<double>>()[mio::AgeGroup(0)]        = params[i++];
    param_set.get<mio::osecirts::CriticalPerSevere<double>>()[mio::AgeGroup(0)]                = params[i++];
    param_set.get<mio::osecirts::DeathsPerCritical<double>>()[mio::AgeGroup(0)]                = params[i++];

    param_set.get<mio::osecirts::ReducExposedPartialImmunity<double>>()[mio::AgeGroup(0)]           = params[i++];
    param_set.get<mio::osecirts::ReducExposedImprovedImmunity<double>>()[mio::AgeGroup(0)]          = params[i++];
    param_set.get<mio::osecirts::ReducInfectedSymptomsPartialImmunity<double>>()[mio::AgeGroup(0)]  = params[i++];
    param_set.get<mio::osecirts::ReducInfectedSymptomsImprovedImmunity<double>>()[mio::AgeGroup(0)] = params[i++];
    param_set.get<mio::osecirts::ReducInfectedSevereCriticalDeadPartialImmunity<double>>()[mio::AgeGroup(0)] =
        params[i++];
    param_set.get<mio::osecirts::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>()[mio::AgeGroup(0)] =
        params[i++];
    param_set.get<mio::osecirts::ReducTimeInfectedMild<double>>()[mio::AgeGroup(0)] = params[i++];

    param_set.apply_constraints();

    // simulate
    auto sim = mio::osecirts::simulate_flows(t0, tmax, dt, model);

    auto result = sim[0];
    auto flows  = sim[1];

    double new_infections = 0.0;
    for (auto t = 1; t < result.get_num_time_points(); ++t) {
        new_infections +=
            flows.get_value(
                t)[model.get_flat_flow_index<mio::osecirts::InfectionState::SusceptibleNaive,
                                             mio::osecirts::InfectionState::ExposedNaive>({mio::AgeGroup(0)})] -
            flows.get_value(
                t - 1)[model.get_flat_flow_index<mio::osecirts::InfectionState::SusceptibleNaive,
                                                 mio::osecirts::InfectionState::ExposedNaive>({mio::AgeGroup(0)})] +
            flows.get_value(t)[model.get_flat_flow_index<mio::osecirts::InfectionState::SusceptiblePartialImmunity,
                                                         mio::osecirts::InfectionState::ExposedPartialImmunity>(
                {mio::AgeGroup(0)})] -
            flows.get_value(t - 1)[model.get_flat_flow_index<mio::osecirts::InfectionState::SusceptiblePartialImmunity,
                                                             mio::osecirts::InfectionState::ExposedPartialImmunity>(
                {mio::AgeGroup(0)})] +
            flows.get_value(t)[model.get_flat_flow_index<mio::osecirts::InfectionState::SusceptibleImprovedImmunity,
                                                         mio::osecirts::InfectionState::ExposedImprovedImmunity>(
                {mio::AgeGroup(0)})] -
            flows.get_value(t - 1)[model.get_flat_flow_index<mio::osecirts::InfectionState::SusceptibleImprovedImmunity,
                                                             mio::osecirts::InfectionState::ExposedImprovedImmunity>(
                {mio::AgeGroup(0)})];
    }

    double icu = 0.0;
    for (auto t = 0; t < result.get_num_time_points(); ++t) {
        icu += result.get_value(t)[model.populations.get_flat_index(
                   {mio::AgeGroup(0), mio::osecirts::InfectionState::InfectedCriticalNaive})] +
               result.get_value(t)[model.populations.get_flat_index(
                   {mio::AgeGroup(0), mio::osecirts::InfectionState::InfectedCriticalPartialImmunity})] +
               result.get_value(t)[model.populations.get_flat_index(
                   {mio::AgeGroup(0), mio::osecirts::InfectionState::InfectedCriticalImprovedImmunity})];
    }

    return icu;
}

std::vector<double> calculate_elementary_effects(mio::osecirts::Model<double>& model,
                                                 const std::vector<double>& base_params, double relative_delta,
                                                 double t0, double tmax, double dt, int num_trajectories)
{
    std::vector<double> mean_effects(base_params.size(), 0.0);
    std::vector<double> std_effects(base_params.size(), 0.0);

    for (int traj = 0; traj < num_trajectories; ++traj) {
        // choose random parameter
        std::vector<double> params = base_params;
        for (size_t i = 0; i < base_params.size(); ++i) {
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<> dis(-relative_delta, relative_delta);
            params[i] += params[i] * dis(gen); // relative change
        }

        std::vector<double> effects(base_params.size(), 0.0);
        for (size_t i = 0; i < base_params.size(); ++i) {
            std::vector<double> params_plus  = params;
            std::vector<double> params_minus = params;

            // Verändere den Parameter in positiver und negativer Richtung
            params_plus[i] += params[i] * relative_delta;
            params_minus[i] -= params[i] * relative_delta;

            double output_plus  = run_model_with_params(model, params_plus, t0, tmax, dt);
            double output_minus = run_model_with_params(model, params_minus, t0, tmax, dt);

            // Elementary Effect
            effects[i] = (output_plus - output_minus) / (2 * params[i] * relative_delta);
        }

        for (size_t i = 0; i < base_params.size(); ++i) {
            mean_effects[i] += effects[i];
            std_effects[i] += effects[i] * effects[i];
        }
    }

    for (size_t i = 0; i < base_params.size(); ++i) {
        mean_effects[i] /= num_trajectories;
        std_effects[i] = std::sqrt((std_effects[i] / num_trajectories) - std::pow(mean_effects[i], 2));
    }

    // return std_effects;

    return mean_effects;
}

int main()
{
    mio::set_log_level(mio::LogLevel::err);

    double t0   = 0;
    double tmax = 200;
    double dt   = 0.1;

    mio::osecirts::Model<double> model(1);
    model.parameters.get<mio::osecirts::ICUCapacity<double>>()          = 100;
    model.parameters.get<mio::osecirts::TestAndTraceCapacity<double>>() = 0.0143;
    const size_t daily_vaccinations                                     = 0;
    const size_t num_days                                               = 300;
    model.parameters.get<mio::osecirts::DailyPartialVaccination<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirts::DailyFullVaccination<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirts::DailyBoosterVaccination<double>>().resize(mio::SimulationDay(num_days));
    model.parameters.get<mio::osecirts::MaxRiskOfInfectionFromSymptomatic<double>>()[mio::AgeGroup(0)] = 0.5;

    for (size_t i = 0; i < num_days; ++i) {
        for (mio::AgeGroup j = 0; j < mio::AgeGroup(1); ++j) {
            auto num_vaccinations = static_cast<double>(i * daily_vaccinations);
            model.parameters.get<mio::osecirts::DailyPartialVaccination<double>>()[{j, mio::SimulationDay(i)}] =
                num_vaccinations;
            model.parameters.get<mio::osecirts::DailyFullVaccination<double>>()[{j, mio::SimulationDay(i)}] =
                num_vaccinations;
            model.parameters.get<mio::osecirts::DailyBoosterVaccination<double>>()[{j, mio::SimulationDay(i)}] =
                num_vaccinations;
        }
    }

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::osecirts::ContactPatterns<double>>();
    const double cont_freq                  = 10;
    const double fact                       = 1.0 / 1.0;
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant((size_t)1, (size_t)1, fact * cont_freq));

    std::vector<double> base_params = {
        3.33, 1.87,  7.0,  6.0,   7.0,  60.0,  60.0, 180.0, 180.0, // Times
        0.15, 0.5,   0.1,  0.2,   0.11, 0.11,  0.11, // Probabilities
        0.8,  0.331, 0.65, 0.243, 0.11, 0.091, 0.9 // Reduction factors
    };

    double delta = 0.5;

    //  Elementary Effects
    std::vector<double> effects = calculate_elementary_effects(model, base_params, delta, t0, tmax, dt, 20000);

    //  results
    std::cout << "[";
    for (size_t i = 0; i < effects.size(); ++i) {
        std::cout << effects[i] << ", ";
        if (i == 21) {
            std::cout << effects[20] << ", ";
            std::cout << effects[21] << ", ";
        }
    }
    std::cout << "]" << std::endl;

    return 0;
}
