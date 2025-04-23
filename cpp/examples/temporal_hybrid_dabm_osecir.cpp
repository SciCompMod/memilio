/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julia Bicker
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

#include "d_abm/model.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "models/hybrid/temporal_hybrid_model.h"
#include "models/hybrid/infection_state.h"
#include "models/d_abm/simulation.h"
#include "models/d_abm/single_well.h"
#include "memilio/compartments/simulation.h"
#include "models/ode_secir/model.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/epidemiology/adoption_rate.h"
#include "memilio/geography/regions.h"
#include "ode_secir/infection_state.h"
#include "models/hybrid/conversion_functions.cpp"
#include <cstddef>
#include <vector>

int main()
{
    mio::set_log_level(mio::LogLevel::warn);
    // Simple example to demonstrate how to run a simulation using temporal-hybrid model combining the diffusive ABM and the ODE-SECIR-model.
    // As condition to switch between models we use a threshold of 20 infected individuals (For <20 Infected the ABM is used and for >=20 Infected the ODE-Model is used).

    using ABM = mio::dabm::Model<SingleWell<mio::hybrid::InfectionState>>;
    using ODE = mio::osecir::Model<double>;

    //Initialize ABM population
    std::vector<ABM::Agent> agents(1000);
    //Random variables used to initialize agents' position and infection state
    auto& pos_sampler  = mio::UniformDistribution<double>::get_instance();
    auto& stat_sampler = mio::DiscreteDistribution<size_t>::get_instance();
    //Infection state distribution
    std::vector<double> infection_state_dist{0.995, 0.005, 0., 0., 0., 0.};
    //Sample agents' position and infection state
    for (auto& a : agents) {
        //Agents' positions are equally distributed in [-2, 2] x [-2, 2]
        a.position = Eigen::Vector2d{pos_sampler(mio::thread_local_rng(), -2., 2.),
                                     pos_sampler(mio::thread_local_rng(), -2., 2.)};
        //Agents' infection states are sampled from infection_state_dist
        a.status =
            static_cast<mio::hybrid::InfectionState>(stat_sampler(mio::thread_local_rng(), infection_state_dist));
    }
    //Transmission parameters used for both models
    double contact_frequency = 10, trans_prob_on_contact = 0.06, time_E = 3., time_Ins = 2.5, time_Isy = 5.2,
           time_Isev = 9., time_Icri = 7.2, mu_Ins_R = 0.2, mu_Isy_Isev = 0.1, mu_Isev_Icri = 0.1, mu_Icri_D = 0.2;
    //Initialize ABM adoption rates
    std::vector<mio::AdoptionRate<mio::hybrid::InfectionState>> adoption_rates;
    //Second-order adoption rate (S->E)
    adoption_rates.push_back(
        {mio::hybrid::InfectionState::Susceptible,
         mio::hybrid::InfectionState::Exposed,
         mio::regions::Region(0),
         contact_frequency * trans_prob_on_contact,
         {{mio::hybrid::InfectionState::InfectedNoSymptoms, 1}, {mio::hybrid::InfectionState::InfectedSymptoms, 1}}});
    //First-order adoption rates
    //E->Ins
    adoption_rates.push_back({mio::hybrid::InfectionState::Exposed,
                              mio::hybrid::InfectionState::InfectedNoSymptoms,
                              mio::regions::Region(0),
                              1. / time_E,
                              {}});
    //Ins->Isy
    adoption_rates.push_back({mio::hybrid::InfectionState::InfectedNoSymptoms,
                              mio::hybrid::InfectionState::InfectedSymptoms,
                              mio::regions::Region(0),
                              (1 - mu_Ins_R) / time_Ins,
                              {}});
    //Ins->R
    adoption_rates.push_back({mio::hybrid::InfectionState::InfectedNoSymptoms,
                              mio::hybrid::InfectionState::Recovered,
                              mio::regions::Region(0),
                              mu_Ins_R / time_Ins,
                              {}});
    //Isy->Isev
    adoption_rates.push_back({mio::hybrid::InfectionState::InfectedSymptoms,
                              mio::hybrid::InfectionState::InfectedSevere,
                              mio::regions::Region(0),
                              mu_Isy_Isev / time_Isy,
                              {}});
    //Isy->R
    adoption_rates.push_back({mio::hybrid::InfectionState::InfectedSymptoms,
                              mio::hybrid::InfectionState::Recovered,
                              mio::regions::Region(0),
                              (1 - mu_Isy_Isev) / time_Isy,
                              {}});
    //Isev->Icri
    adoption_rates.push_back({mio::hybrid::InfectionState::InfectedSevere,
                              mio::hybrid::InfectionState::InfectedCritical,
                              mio::regions::Region(0),
                              mu_Isev_Icri / time_Isev,
                              {}});
    //Isev->R
    adoption_rates.push_back({mio::hybrid::InfectionState::InfectedSevere,
                              mio::hybrid::InfectionState::Recovered,
                              mio::regions::Region(0),
                              (1 - mu_Isev_Icri) / time_Isev,
                              {}});
    //Icri->R
    adoption_rates.push_back({mio::hybrid::InfectionState::InfectedCritical,
                              mio::hybrid::InfectionState::Recovered,
                              mio::regions::Region(0),
                              (1 - mu_Icri_D) / time_Icri,
                              {}});
    //Icri->D
    adoption_rates.push_back({mio::hybrid::InfectionState::InfectedCritical,
                              mio::hybrid::InfectionState::Dead,
                              mio::regions::Region(0),
                              mu_Icri_D / time_Icri,
                              {}});
    //Interaction radius and noise
    double interaction_radius = 0.4, noise = 0.5;
    ABM abm(agents, adoption_rates, interaction_radius, noise,
            {mio::hybrid::InfectionState::InfectedSevere, mio::hybrid::InfectionState::InfectedCritical,
             mio::hybrid::InfectionState::Dead});

    //As we start modeling with the ABM, we don't need to initialize the population for the ODE-model
    //Initialize ODE model parameters
    ODE ode(1);
    ode.parameters.get<mio::osecir::TimeExposed<double>>()[mio::AgeGroup(0)]            = time_E;
    ode.parameters.get<mio::osecir::TimeInfectedNoSymptoms<double>>()[mio::AgeGroup(0)] = time_Ins;
    ode.parameters.get<mio::osecir::TimeInfectedSymptoms<double>>()[mio::AgeGroup(0)]   = time_Isy;
    ode.parameters.get<mio::osecir::TimeInfectedSevere<double>>()[mio::AgeGroup(0)]     = time_Isev;
    ode.parameters.get<mio::osecir::TimeInfectedCritical<double>>()[mio::AgeGroup(0)]   = time_Icri;
    ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact<double>>()[mio::AgeGroup(0)] =
        trans_prob_on_contact;
    ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<double>>()[mio::AgeGroup(0)] = mu_Ins_R;
    ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<double>>()[mio::AgeGroup(0)]      = mu_Isy_Isev;
    ode.parameters.get<mio::osecir::CriticalPerSevere<double>>()[mio::AgeGroup(0)]              = mu_Isev_Icri;
    ode.parameters.get<mio::osecir::DeathsPerCritical<double>>()[mio::AgeGroup(0)]              = mu_Icri_D;
    ode.apply_constraints();
    mio::ContactMatrixGroup& contact_matrix = ode.parameters.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, contact_frequency));

    //Set t0 and internal dt for each model
    double t0 = 0;
    double dt = 0.1;

    //Create simulations
    auto sim_abm = mio::dabm::Simulation(abm, t0, dt);
    auto sim_ode = mio::Simulation(ode, t0, dt);

    const auto result_abm = [](const mio::dabm::Simulation<SingleWell<mio::hybrid::InfectionState>>& sim,
                               double /*t*/) {
        return sim.get_result();
    };

    const auto result_ode = [](const mio::Simulation<double, ODE>& sim, double /*t*/) {
        return sim.get_result();
    };

    //Create hybrid simulation
    double dt_switch = 0.2;
    mio::hybrid::TemporalHybridSimulation<decltype(sim_abm), decltype(sim_ode), mio::TimeSeries<double>,
                                          mio::TimeSeries<double>>
        hybrid_sim(sim_abm, sim_ode, result_abm, result_ode, true, t0, dt_switch);

    //Define switching conditiond
    const auto condition = [](const mio::TimeSeries<double>& result_abm, const mio::TimeSeries<double>& result_ode,
                              bool abm_used) {
        if (abm_used) {
            auto& last_value = result_abm.get_last_value().eval();
            if ((last_value[(int)mio::hybrid::InfectionState::Exposed] +
                 last_value[(int)mio::hybrid::InfectionState::InfectedNoSymptoms] +
                 last_value[(int)mio::hybrid::InfectionState::InfectedSymptoms] +
                 last_value[(int)mio::hybrid::InfectionState::InfectedSevere] +
                 last_value[(int)mio::hybrid::InfectionState::InfectedCritical]) > 20) {
                return true;
            }
        }
        else {
            auto& last_value = result_ode.get_last_value().eval();
            if ((last_value[(int)mio::osecir::InfectionState::Exposed] +
                 last_value[(int)mio::osecir::InfectionState::InfectedNoSymptoms] +
                 last_value[(int)mio::osecir::InfectionState::InfectedNoSymptomsConfirmed] +
                 last_value[(int)mio::osecir::InfectionState::InfectedSymptoms] +
                 last_value[(int)mio::osecir::InfectionState::InfectedSymptomsConfirmed] +
                 last_value[(int)mio::osecir::InfectionState::InfectedSevere] +
                 last_value[(int)mio::osecir::InfectionState::InfectedCritical]) <= 20) {
                return true;
            }
        }
        return false;
    };

    //Simulate for 30 days
    hybrid_sim.advance(30., condition);

    //Print result time series of both models
    auto ts_abm = hybrid_sim.get_result_model1();
    auto ts_ode = hybrid_sim.get_result_model2();

    ts_abm.print_table({"S", "E", "Ins", "Isy", "Isev", "Icri", "R", "D"});
    ts_ode.print_table({"S", "E", "Ins", "Ins_confirmed", "Isy", "Isy_confirmed", "Isev", "Icri", "R", "D"});
}
