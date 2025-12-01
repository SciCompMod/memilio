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

#include "ode_seirv/model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

/**
 * @brief set_initial_population sets the initial population of the model 
 * 
 * We assume that no one is initially infected or exposed at the beginning of the season.
 * Infections are seeded later via the OutsideFoI parameter.
 * The distinction into the 2 layers of susceptibility is done via the parameters.
 * @tparam FP floating point type
 * @param[in, out] model SEIRV model
 * @param[in] total_pop total population size
 */
template <class FP>
void set_initial_population(mio::oseirv::Model<FP>& model, const FP total_pop)
{
    auto& params            = model.parameters;
    auto& pop               = model.populations;
    const size_t num_groups = (size_t)params.get_num_groups();

    for (size_t i = 0; i < num_groups; ++i) {
        pop.template set_difference_from_group_total<mio::AgeGroup>(
            {mio::AgeGroup(i), mio::oseirv::InfectionState::Susceptible}, total_pop / num_groups);

        // Total population N_i as currently stored in group totals
        FP Ni = pop.get_group_total(mio::AgeGroup(i));

        FP s_age  = params.template get<mio::oseirv::SusceptibilityByAge<FP>>()[mio::AgeGroup(i)];
        FP s_frac = params.template get<mio::oseirv::SusceptibleFraction<FP>>();
        FP vc     = params.template get<mio::oseirv::VaccineCoverage<FP>>()[mio::AgeGroup(i)];
        FP ve     = params.template get<mio::oseirv::VaccineEffectiveness<FP>>()[mio::AgeGroup(i)];

        pop[{mio::AgeGroup(i), mio::oseirv::InfectionState::Exposed}]            = 0;
        pop[{mio::AgeGroup(i), mio::oseirv::InfectionState::ExposedVaccinated}]  = 0;
        pop[{mio::AgeGroup(i), mio::oseirv::InfectionState::Infected}]           = 0;
        pop[{mio::AgeGroup(i), mio::oseirv::InfectionState::InfectedVaccinated}] = 0;

        pop[{mio::AgeGroup(i), mio::oseirv::InfectionState::Susceptible}] = s_frac * s_age * (FP(1) - vc) * Ni;
        pop[{mio::AgeGroup(i), mio::oseirv::InfectionState::SusceptibleVaccinated}] =
            s_frac * s_age * (FP(1) - ve) * vc * Ni;

        pop[{mio::AgeGroup(i), mio::oseirv::InfectionState::Recovered}] = (FP(1) - s_frac * s_age) * (FP(1) - vc) * Ni;
        pop[{mio::AgeGroup(i), mio::oseirv::InfectionState::RecoveredVaccinated}] =
            (FP(1) - s_frac * s_age * (FP(1) - ve)) * vc * Ni;
    }
}

// Example usage of the SEIRV ODE model presented in https://doi.org/10.1186/s12879-017-2344-6
int main()
{
    using FP = double;
    mio::set_log_level(mio::LogLevel::debug);

    const FP t0 = 0., tmax = 42.;
    const FP dt = 0.1;

    const size_t num_groups = 1;
    const auto total_pop    = 1e5;
    mio::oseirv::Model<FP> model((int)num_groups);
    auto& parameters = model.parameters;

    FP cont_freq       = 10.0;
    FP cont_freq_group = FP(1) / FP(num_groups);

    // Healthy contacts
    mio::ContactMatrixGroup<ScalarType>& contacts_healthy = parameters.get<mio::oseirv::ContactPatternsHealthy<FP>>();
    contacts_healthy[0]                                   = mio::ContactMatrix<FP>(
        Eigen::MatrixXd::Constant((Eigen::Index)num_groups, (Eigen::Index)num_groups, cont_freq_group * cont_freq));

    // Sick contacts (here 20% reduced)
    mio::ContactMatrixGroup<ScalarType>& contacts_sick = parameters.get<mio::oseirv::ContactPatternsSick<FP>>();
    contacts_sick[0]                                   = mio::ContactMatrix<FP>(Eigen::MatrixXd::Constant(
        (Eigen::Index)num_groups, (Eigen::Index)num_groups, cont_freq_group * cont_freq * 0.8));

    // Parameters
    parameters.get<mio::oseirv::BaselineTransmissibility<FP>>()   = 1.2;
    parameters.get<mio::oseirv::TimeExposed<FP>>()                = 2.0;
    parameters.get<mio::oseirv::TimeInfected<FP>>()               = 2.0; // same duration for E->I and I->R
    parameters.get<mio::oseirv::SeasonalityAmplitude<FP>>()       = 1.0;
    parameters.get<mio::oseirv::SeasonalityShiftPerSubtype<FP>>() = 0.0;
    parameters.get<mio::oseirv::SeasonalityShiftPerSeason<FP>>()  = 0.0;
    parameters.get<mio::oseirv::OutsideFoI<FP>>()                 = 1e-6;
    parameters.get<mio::oseirv::ClusteringExponent<FP>>()         = 0.9;
    parameters.get<mio::oseirv::SickMixing<FP>>()                 = 2.0;

    for (size_t i = 0; i < num_groups; ++i) {
        parameters.get<mio::oseirv::SusceptibilityByAge<FP>>()[mio::AgeGroup(i)]  = 1.0;
        parameters.get<mio::oseirv::VaccineCoverage<FP>>()[mio::AgeGroup(i)]      = 0.3;
        parameters.get<mio::oseirv::VaccineEffectiveness<FP>>()[mio::AgeGroup(i)] = 0.5;
        model.populations.set_difference_from_group_total<mio::AgeGroup>(
            {mio::AgeGroup(i), mio::oseirv::InfectionState::Susceptible}, 1e5 / num_groups);
    }

    set_initial_population(model, total_pop);

    auto seirv = mio::simulate<FP>(t0, tmax, dt, model);

    std::vector<std::string> vars = {"S", "SV", "E", "EV", "I", "IV", "R", "RV"};
    seirv.print_table(vars);
}
