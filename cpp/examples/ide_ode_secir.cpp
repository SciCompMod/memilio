/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Anna Wendler
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
#include "boost/numeric/odeint/stepper/controlled_runge_kutta.hpp"
#include "ode_secir/infection_state.h"
#include "ode_secir/model.h"
#include "memilio/math/adapt_rk.h"
#include "ode_secir/parameter_space.h"
#include "ode_secir/analyze_result.h"
#include "ode_secir/parameters.h"

#include "boost/fusion/functional/invocation/invoke.hpp"
#include "ide_secir/infection_state.h"
#include "ide_secir/model_ide.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"
#include "ide_secir/initialize_from_ode.h"
#include "memilio/io/result_io.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include <iostream>
#include <string>

int main()
{

    // Here we decide what exactly we want to do in the example below
    bool print_to_terminal = true;
    bool save_result       = true;
    bool ide_simulation    = false;
    int dt_ode_exponent    = 6;
    int dt_ide_exponent    = 0;
    // We use setting 2 as baseline, changes for other settings are in respective if statements
    int setting = 6;

    // General set up.
    ScalarType t0       = 0;
    ScalarType tmax     = 40.00;
    ScalarType dt_ode   = pow(10, -dt_ode_exponent);
    ScalarType dt_ide   = pow(10, -dt_ide_exponent);
    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    /**********************************
    *         ODE simulation          *
    **********************************/

    ScalarType nb_total_t0 = 10000, nb_exp_t0 = 20, nb_car_t0 = 20, nb_inf_t0 = 3, nb_hosp_t0 = 1, nb_icu_t0 = 1,
               nb_rec_t0 = 10, nb_dead_t0 = 0;

    if (setting == 10 || setting == 12 || setting == 13 || setting == 14 || setting == 19 || setting == 25 ||
        setting == 26 || setting == 28 || setting == 29 || setting == 30 || setting == 31) {
        nb_rec_t0 = 0.;
    }

    if (setting == 20) {
        nb_rec_t0 = 1000.;
    }

    if (setting == 24) {
        nb_rec_t0 = 1.;
    }

    if (setting == 27) {
        nb_rec_t0 = 0.1;
    }

    mio::osecir::Model model_ode(1);

    // Set parameters
    ScalarType cont_freq = 1.0;

    // Set Seasonality=0 so that cont_freq_eff is equal to contact_matrix
    model_ode.parameters.set<mio::osecir::Seasonality>(0.0);
    mio::ContactMatrixGroup& contact_matrix = model_ode.parameters.get<mio::osecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    model_ode.parameters.get<mio::osecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    // Parameters needed to determine transition rates
    // model_ode.parameters.get<mio::osecir::IncubationTime>()[(mio::AgeGroup)0] = 2.6;
    // model_ode.parameters.get<mio::osecir::SerialInterval>()[(mio::AgeGroup)0] = 2.0;
    model_ode.parameters.get<mio::osecir::TimeExposed>()[(mio::AgeGroup)0]            = 1.4;
    model_ode.parameters.get<mio::osecir::TimeInfectedNoSymptoms>()[(mio::AgeGroup)0] = 1.2;
    model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0]   = 0.3;
    model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]     = 0.3;
    model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0]   = 0.3;

    // Set initial values for compartments
    model_ode.populations.set_total(nb_total_t0);
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]                     = nb_exp_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]          = nb_car_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]            = nb_inf_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]              = nb_hosp_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]            = nb_icu_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]                   = nb_rec_t0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]                        = nb_dead_t0;
    model_ode.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                    nb_total_t0);

    // Set probabilities that determine proportion between compartments
    model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 0.5;
    model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 0.5;
    model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 0.5;
    model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 0.5;

    if (setting == 4) {
        // Set probabilities that determine proportion between compartments
        model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 0.;
        model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 1.;
        model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 1.;
        model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 0.5;
    }

    if (setting == 5 || setting == 12 || setting == 13 || setting == 15 || setting == 20 || setting == 21 ||
        setting == 23 || setting == 25 || setting == 26 || setting == 28 || setting == 29) {
        // Set probabilities that determine proportion between compartments
        model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 0.;
        model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 1.;
        model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 1.;
        model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 1.;
    }

    if (setting == 30) {
        // Set probabilities that determine proportion between compartments
        model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 0.;
        model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 1.;
        model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 1.;
        model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 0.99;
    }

    if (setting == 31) {
        // Set probabilities that determine proportion between compartments
        model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 0.;
        model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 1.;
        model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 1.;
        model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 0.95;
    }

    if (setting == 6 || setting == 11 || setting == 14) {
        // Set probabilities that determine proportion between compartments
        model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 0.;
        model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 1.;
        model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 1.;
        model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 0.;
    }

    if (setting == 7) {
        // Set probabilities that determine proportion between compartments
        model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0] = 1.;
        model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]      = 0.;
        model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]              = 0.;
        model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]              = 0.;
    }

    if (setting == 8) {
        // Set probabilities that determine proportion between compartments
        model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0] = 0.;
    }

    // Further model parameters
    model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0] = 1.0;
    model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]   = 1.0;
    model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]   = 1.0;
    // Choose TestAndTraceCapacity very large so that riskFromInfectedSymptomatic = RiskOfInfectionFromSymptomatic
    model_ode.parameters.get<mio::osecir::TestAndTraceCapacity>() = std::numeric_limits<ScalarType>::max();
    // Choose ICUCapacity very large so that CriticalPerSevereAdjusted = CriticalPerSevere and deathsPerSevereAdjusted = 0
    model_ode.parameters.get<mio::osecir::ICUCapacity>() = std::numeric_limits<ScalarType>::max();

    model_ode.check_constraints();

    // TODO: find out how we can change the integrator to another one from boost that doesn't use adaptive time steps
    // auto integrator = std::make_shared<mio::RKIntegratorCore>();
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // auto integrator = std::make_shared<boost::numeric::odeint::runge_kutta_cash_karp54>();
    // choose dt_min = dt_max so that we have a fixed time step and can compare to IDE
    integrator->set_dt_min(dt_ode);
    integrator->set_dt_max(dt_ode);
    // set tolerance as follows so that every time step is only computed once (found out by trying)
    integrator->set_rel_tolerance(1e-1);
    integrator->set_abs_tolerance(1e-1);

    mio::TimeSeries<ScalarType> secihurd_ode = simulate(t0, tmax, dt_ode, model_ode, integrator);

    // Compute flows from ODE result.
    // Note that we are computing \tilde{\sigma} here. To be able to compare flows between different timetspes (of ODE and IDE)
    // we need to divide by dt to get \hat{\sigma}. This is not done here but in the python scripts for the analysis of results.
    mio::TimeSeries<ScalarType> secihurd_ode_flows(num_transitions);
    mio::isecir::get_flows_from_ode_compartments(model_ode, secihurd_ode, secihurd_ode_flows, tmax, tmax - t0, dt_ode);

    if (print_to_terminal) {
        secihurd_ode.print_table();
    }

    std::cout << "\n";

    if (save_result) {
        auto save_result_status_ode =
            mio::save_result({secihurd_ode}, {0}, 1,
                             "../../results/result_ode_dt=1e-" + std::to_string(dt_ode_exponent) + "_setting" +
                                 std::to_string(setting) + ".h5");
        auto save_result_status_ode_flows =
            mio::save_result({secihurd_ode_flows}, {0}, 1,
                             "../../results/result_ode_flows_dt=1e-" + std::to_string(dt_ode_exponent) + "_setting" +
                                 std::to_string(setting) + ".h5");
    }

    /**********************************
    *         IDE simulation          *
    **********************************/

    if (ide_simulation) {
        using Vec = mio::TimeSeries<ScalarType>::Vector;

        ScalarType N = nb_total_t0;
        // TODO: Set this automatically or check if this is possible (wrt to global_max_support) with the given ODE simulation
        ScalarType t0_ide = 35.0;
        // Get number of dead individuals at time t0_ide from ODE model
        ScalarType deaths = secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt_ode -
                                         1][(int)mio::osecir::InfectionState::Dead];

        mio::TimeSeries<ScalarType> init_transitions(num_transitions);

        ScalarType total_infections = 0.;

        if (setting == 11 || setting == 13 || setting == 14 || setting == 18 || setting == 21) {
            // Compute total_infections by getting number of individuals that currently are or have been infected at time t0_ide
            total_infections =
                secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt_ode - 1]
                            [(int)mio::osecir::InfectionState::InfectedSymptoms] +
                secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt_ode - 1]
                            [(int)mio::osecir::InfectionState::InfectedSevere] +
                secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt_ode - 1]
                            [(int)mio::osecir::InfectionState::InfectedCritical] +
                secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt_ode - 1]
                            [(int)mio::osecir::InfectionState::Recovered] +
                secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt_ode - 1]
                            [(int)mio::osecir::InfectionState::Dead];
        }

        // Initialize model.
        mio::isecir::Parameters parameters;
        mio::isecir::Model model_ide(std::move(init_transitions), N, deaths, total_infections, parameters);

        // Set working parameters.

        // Contact matrix; contact_matrix was already defined for ODE
        model_ide.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        // To compare with the ODE model we use ExponentialDecay functions as TransitionDistributions
        // We set the parameters so that they correspond to the above ODE model
        mio::ExponentialSurvivalFunction exponential(10.0);
        mio::StateAgeFunctionWrapper delaydistribution(exponential);
        std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
        // ExposedToInfectedNoSymptoms
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms].set_distribution_parameter(
            1 / model_ode.parameters.get<mio::osecir::TimeExposed>()[(mio::AgeGroup)0]);
        // InfectedNoSymptomsToInfectedSymptoms
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
            .set_distribution_parameter(
                1 / model_ode.parameters.get<mio::osecir::TimeInfectedNoSymptoms>()[(mio::AgeGroup)0]);
        // InfectedNoSymptomsToRecovered
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]
            .set_distribution_parameter(
                1 / model_ode.parameters.get<mio::osecir::TimeInfectedNoSymptoms>()[(mio::AgeGroup)0]);
        // InfectedSymptomsToInfectedSevere
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]
            .set_distribution_parameter(
                1 / model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0]);
        // InfectedSymptomsToRecovered
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered].set_distribution_parameter(
            1 / model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms>()[(mio::AgeGroup)0]);
        // InfectedSevereToInfectedCritical
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]
            .set_distribution_parameter(1 /
                                        model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]);
        // InfectedSevereToRecovered
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered].set_distribution_parameter(
            1 / model_ode.parameters.get<mio::osecir::TimeInfectedSevere>()[(mio::AgeGroup)0]);
        // InfectedCriticalToDead]
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead].set_distribution_parameter(
            1 / model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0]);
        // InfectedCriticalToRecovered
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_distribution_parameter(
            1 / model_ode.parameters.get<mio::osecir::TimeInfectedCritical>()[(mio::AgeGroup)0]);

        model_ide.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

        // Set probabilities that determine proportion between compartments
        std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.);
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
            1 - model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
            model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
            model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)] =
            1 - model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
            model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)] =
            1 - model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] =
            model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)] =
            1 - model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0];
        model_ide.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

        // Set further parameters
        mio::ConstantFunction constfunc_proboncontact(
            model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact>()[(mio::AgeGroup)0]);
        mio::StateAgeFunctionWrapper proboncontact(constfunc_proboncontact);
        model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(proboncontact);

        mio::ConstantFunction constfunc_reltransnosympt(
            model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms>()[(mio::AgeGroup)0]);
        mio::StateAgeFunctionWrapper reltransnosympt(constfunc_reltransnosympt);
        model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(reltransnosympt);

        mio::ConstantFunction constfunc_riskofinf(
            model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic>()[(mio::AgeGroup)0]);
        mio::StateAgeFunctionWrapper riskofinf(constfunc_riskofinf);
        model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(riskofinf);

        // Compute initial flows from results of ODE simulation
        mio::isecir::compute_initial_flows_for_ide_from_ode(model_ode, model_ide, secihurd_ode, t0_ide, dt_ode, dt_ide);

        model_ide.check_constraints(dt_ide);

        // model_ide.set_populations_before_simulation();

        if (setting == 15 || setting == 16 || setting == 25) {
            model_ide.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Recovered] =
                secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt_ode - 1]
                            [(int)mio::osecir::InfectionState::Recovered];
        }

        if (setting == 17 || setting == 26) {
            model_ide.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Susceptible] =
                secihurd_ode[(Eigen::Index)secihurd_ode.get_num_time_points() - (tmax - t0_ide) / dt_ode - 1]
                            [(int)mio::osecir::InfectionState::Susceptible];
        }

        // Carry out simulation
        std::cout << "Simulating now \n";
        mio::isecir::Simulation sim(model_ide, dt_ide);
        sim.advance(tmax);

        mio::TimeSeries<ScalarType> secihurd_ide       = sim.get_result();
        mio::TimeSeries<ScalarType> secihurd_ide_flows = sim.get_transitions();

        if (print_to_terminal) {
            secihurd_ide.print_table();
            secihurd_ide_flows.print_table();
        }

        std::cout << "Initialization method: " << sim.get_model().get_initialization_method_compartments() << "\n";

        // std::cout << "Compartments at last time step of ODE:\n";
        // std::cout << "# time  |  S  |  E  |  C  |  I  |  H  |  U  |  R  |  D  |" << std::endl;
        // for (Eigen::Index j = 0; j < secihurd_ode.get_num_elements(); ++j) {
        //     std::cout << "  |  " << std::fixed << std::setprecision(8) << secihurd_ode.get_last_value()[j];
        // }
        // std::cout << "\n" << std::endl;

        // std::cout << "Compartments at last time step of IDE:\n";
        // std::cout << "# time  |  S  |  E  |  C  |  I  |  H  |  U  |  R  |  D  |" << std::endl;
        // for (Eigen::Index j = 0; j < secihurd_ide.get_num_elements(); ++j) {
        //     std::cout << "  |  " << std::fixed << std::setprecision(8) << secihurd_ide.get_last_value()[j];
        // }
        std::cout << "\n" << std::endl;

        // if (save_result) {

        //     auto save_result_status_ide = mio::save_result(
        //         {secihurd_ide}, {0}, 1,
        //         "../../results/result_ide_dt=1e-" + std::to_string(dt_ide_exponent) + "_init_dt_ode=1e-" +
        //             std::to_string(dt_ode_exponent) + "_setting" + std::to_string(setting) + ".h5");
        //     auto save_result_status_ide_flows = mio::save_result(
        //         {secihurd_ide_flows}, {0}, 1,
        //         "../../results/result_ide_flows_dt=1e-" + std::to_string(dt_ide_exponent) + "_init_dt_ode=1e-" +
        //             std::to_string(dt_ode_exponent) + "_setting" + std::to_string(setting) + ".h5");
        // }
    }
}
