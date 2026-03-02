/*
* Copyright (C) 2020-2026 MEmilio
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
#include "ode_secir/model.h"
#include "memilio/compartments/feedback_simulation.h"
#include "memilio/utils/logging.h"

void initialize_model(mio::osecir::Model<double>& model, int total_population, double cont_freq)
{
    model.parameters.set<mio::osecir::StartDay<double>>(60);
    model.parameters.set<mio::osecir::Seasonality<double>>(0.2);

    // time-related parameters
    model.parameters.get<mio::osecir::TimeExposed<double>>()            = 3.2;
    model.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>() = 2.0;
    model.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()   = 5.8;
    model.parameters.get<mio::osecir::TimeInfectedSevere<double>>()     = 9.5;
    model.parameters.get<mio::osecir::TimeInfectedCritical<double>>()   = 7.1;

    // Set transmission and isolation parameters
    model.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()  = 0.05;
    model.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<double>>()    = 0.7;
    model.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()    = 0.09;
    model.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<double>>()    = 0.25;
    model.parameters.get<mio::osecir::MaxRiskOfInfectionFromSymptomatic<double>>() = 0.45;
    model.parameters.get<mio::osecir::TestAndTraceCapacity<double>>()              = 35;
    model.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()         = 0.2;
    model.parameters.get<mio::osecir::CriticalPerSevere<double>>()                 = 0.25;
    model.parameters.get<mio::osecir::DeathsPerCritical<double>>()                 = 0.3;

    // contact matrix
    mio::ContactMatrixGroup<double>& contact_matrix = model.parameters.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0] = mio::ContactMatrix<double>(Eigen::MatrixXd::Constant(1, 1, cont_freq));

    // initial population
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]                     = 40;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]          = 30;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]            = 20;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]              = 10;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]            = 5;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]                   = 20;
    model.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]                        = 0;
    model.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                total_population);
    // The function apply_constraints() ensures that all parameters are within their defined bounds.
    // Note that negative values are set to zero instead of stopping the simulation.
    model.apply_constraints();
}

void initialize_feedback(mio::FeedbackSimulation<double, mio::Simulation<double, mio::osecir::Model<double>>,
                                                 mio::osecir::ContactPatterns<double>>& feedback_simulation)
{
    // nominal ICU capacity
    feedback_simulation.get_parameters().template get<mio::NominalICUCapacity<double>>() = 10;

    // ICU occupancy in the past for memory kernel
    auto& icu_occupancy     = feedback_simulation.get_parameters().template get<mio::ICUOccupancyHistory<double>>();
    Eigen::VectorXd icu_day = Eigen::VectorXd::Constant(1, 1);
    const auto cutoff       = static_cast<int>(feedback_simulation.get_parameters().template get<mio::GammaCutOff>());
    for (int t = -cutoff; t <= 0; ++t) {
        icu_occupancy.add_time_point(t, icu_day);
    }

    // bounds for contact reduction measures
    feedback_simulation.get_parameters().template get<mio::ContactReductionMin<double>>() = {0.1};
    feedback_simulation.get_parameters().template get<mio::ContactReductionMax<double>>() = {0.8};
}

int main()
{
    // This example demonstrates the implementation of a feedback mechanism for a ODE SECIR model.
    // It shows how the perceived risk dynamically impacts contact reduction measures.
    // The feedback mechanism adjusts contact rates during simulation based on the perceived
    // risk which is calculated from the ICU occupancy using a memory kernel.
    mio::set_log_level(mio::LogLevel::warn);

    const double tmax          = 35;
    const int total_population = 1000;
    const double cont_freq     = 10;

    // create and initialize ODE model for a single age group
    mio::osecir::Model<double> model(1);
    initialize_model(model, total_population, cont_freq);

    // determine the index for the ICU state (InfectedCritical) for feedback mechanism
    auto icu_index = std::vector<size_t>{
        model.populations.get_flat_index({mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical})};

    // create simulation objects: first a secir simulation, then a feedback simulation
    auto simulation = mio::osecir::Simulation<double, mio::Simulation<double, mio::osecir::Model<double>>>(model);
    auto feedback_simulation =
        mio::FeedbackSimulation<double, mio::Simulation<double, mio::osecir::Model<double>>,
                                mio::osecir::ContactPatterns<double>>(std::move(simulation), icu_index);

    // set up the parameters for the feedback simulation
    initialize_feedback(feedback_simulation);

    // run the simulation with feedback mechanism
    feedback_simulation.advance(tmax);

    // print the perceived risk and the final total population
    auto& perceived_risk = feedback_simulation.get_perceived_risk();
    perceived_risk.print_table({"Perceived Risk"});
    return 0;
}
