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
#include "ode_secir/infection_state.h"
#include "ode_secir/model.h"
#include "ode_secir/parameters.h"

#include "ide_secir/infection_state.h"
#include "ide_secir/model_ide.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"
#include "ide_secir/initialize_from_ode.h"

#include "memilio/io/result_io.h"
#include "memilio/utils/time_series.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/floating_point.h"
#include <iostream>
#include <string>

mio::TimeSeries<ScalarType> remove_time_points(const mio::TimeSeries<ScalarType>& simulation_result,
                                               ScalarType saving_dt, ScalarType scale = 1)
{
    mio::TimeSeries<ScalarType> removed(simulation_result.get_num_elements());
    ScalarType time = simulation_result.get_time(0);
    removed.add_time_point(time, scale * simulation_result[0]);
    time += saving_dt;
    for (int i = 1; i < simulation_result.get_num_time_points(); i++) {
        if (mio::floating_point_greater_equal(simulation_result.get_time(i), time, 1e-8)) {
            removed.add_time_point(simulation_result.get_time(i), scale * simulation_result[i]);
            time += saving_dt;
        }
    }
    return removed;
}

int main()
{
    // Here we decide what exactly we want to do in the example below.
    bool print_to_terminal = false;
    // if save_exponent is set to 0, the simulation result is not saved.
    ScalarType save_exponent = 4;
    ScalarType saving_dt     = pow(10, -save_exponent);
    bool save_ide            = true;
    // Directory where results will be stored.
    std::string result_dir = "../../results/";

    // Decide if we want to run the IDE simulation or just the ODE simulation for comparison.
    bool ide_simulation = false;

    // General set up.
    ScalarType t0           = 0.;
    ScalarType tmax         = 70.;
    ScalarType ode_exponent = 4;
    ScalarType dt_ode       = pow(10, -ode_exponent);
    ScalarType ide_exponent = 4;
    ScalarType dt_ide       = pow(10, -ide_exponent);
    int num_transitions     = (int)mio::isecir::InfectionTransition::Count;

    /**********************************
    *         ODE simulation          *
    **********************************/

    ScalarType nb_total_t0 = 10000, nb_exp_t0 = 20, nb_car_t0 = 20, nb_inf_t0 = 3, nb_hosp_t0 = 1, nb_icu_t0 = 1,
               nb_rec_t0 = 10, nb_dead_t0 = 0;

    mio::osecir::Model<ScalarType> model_ode(1);

    // Set parameters.
    ScalarType cont_freq = 1.0;

    // Set Seasonality=0 so that cont_freq_eff is equal to contact_matrix.
    model_ode.parameters.set<mio::osecir::Seasonality<ScalarType>>(0.0);
    mio::ContactMatrixGroup& contact_matrix = model_ode.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    model_ode.parameters.get<mio::osecir::ContactPatterns<ScalarType>>() = mio::UncertainContactMatrix(contact_matrix);

    // Parameters needed to determine transition rates.
    model_ode.parameters.get<mio::osecir::TimeExposed<ScalarType>>()[(mio::AgeGroup)0]            = 1.4;
    model_ode.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] = 1.2;
    model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0]   = 0.3;
    model_ode.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[(mio::AgeGroup)0]     = 0.3;
    model_ode.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[(mio::AgeGroup)0]   = 0.3;

    // Set initial values for compartments.
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

    // Set probabilities that determine proportion between compartments.
    model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] = 0.5;
    model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0]      = 0.5;
    model_ode.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0]              = 0.5;
    model_ode.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0]              = 0.5;

    // Further model parameters.
    model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[(mio::AgeGroup)0] = 1.0;
    model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[(mio::AgeGroup)0]   = 1.0;
    model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[(mio::AgeGroup)0]   = 1.0;
    // Choose TestAndTraceCapacity very large so that riskFromInfectedSymptomatic = RiskOfInfectionFromSymptomatic.
    model_ode.parameters.get<mio::osecir::TestAndTraceCapacity<ScalarType>>() = std::numeric_limits<ScalarType>::max();
    // Choose ICUCapacity very large so that CriticalPerSevereAdjusted = CriticalPerSevere and deathsPerSevereAdjusted = 0.
    model_ode.parameters.get<mio::osecir::ICUCapacity<ScalarType>>() = std::numeric_limits<ScalarType>::max();

    // TODO: find out how we can change the integrator to another one from boost that doesn't use adaptive time steps
    auto integrator =
        std::make_shared<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max so that we have a fixed time step and can compare to IDE.
    integrator->set_dt_min(dt_ode);
    integrator->set_dt_max(dt_ode);
    // Set tolerance as follows so that every time step is only computed once (found out by trying).
    integrator->set_rel_tolerance(1e-1);
    integrator->set_abs_tolerance(1e-1);

    std::cout << "Starting simulation with ODE model. \n";
    mio::TimeSeries<ScalarType> secihurd_ode =
        mio::osecir::simulate<ScalarType>(t0, tmax, dt_ode, model_ode, integrator);

    // Compute flows from ODE result to store results.
    // Note that we are computing \tilde{\sigma} here. To be able to compare flows between different timetspes (of ODE and IDE)
    // we need to divide by dt to get \hat{\sigma}. This is not done here but in the python scripts for the analysis of results.
    mio::TimeSeries<ScalarType> secihurd_ode_flows(num_transitions);
    mio::isecir::get_flows_from_ode_compartments(model_ode, secihurd_ode, secihurd_ode_flows, tmax, tmax - t0, dt_ode);

    if (print_to_terminal) {
        secihurd_ode.print_table();
        std::cout << "\n";
    }
    if (save_exponent != 0) {
        // Create directory "results" if not existent yet.
        boost::filesystem::path res_dir(result_dir);
        boost::filesystem::create_directory(res_dir);
        auto save_result_status_ode =
            mio::save_result({remove_time_points(secihurd_ode, saving_dt)}, {0}, 1,
                             result_dir + "result_ode_dt=1e-" + fmt::format("{:.0f}", ode_exponent) + "_savefrequency" +
                                 fmt::format("{:.0f}", save_exponent) + ".h5");
        auto save_result_status_ode_flows =
            mio::save_result({remove_time_points(secihurd_ode_flows, saving_dt, 1. / dt_ode)}, {0}, 1,
                             result_dir + "result_ode_flows_dt=1e-" + fmt::format("{:.0f}", ode_exponent) +
                                 "_savefrequency" + fmt::format("{:.0f}", save_exponent) + ".h5");
        std::cout << "Successfully saved the ODE simulation results. \n";
    }

    /**********************************
    *         IDE simulation          *
    **********************************/

    if (ide_simulation) {
        using Vec = mio::TimeSeries<ScalarType>::Vector;

        ScalarType N = nb_total_t0;
        // TODO: Set this automatically or check if this is possible (wrt to global_max_support) with the given ODE simulation
        ScalarType t0_ide = 35.0;
        // Number of deaths according to the ODE model will be set in the function where also the transitions are calculated.
        ScalarType deaths = 0;

        mio::TimeSeries<ScalarType> init_transitions(num_transitions);

        ScalarType total_infections = 0.;

        // Initialize model.
        mio::isecir::Parameters parameters;
        mio::isecir::Model model_ide(std::move(init_transitions), N, deaths, total_infections, parameters);

        // Set working parameters.

        // Contact matrix; contact_matrix was already defined for ODE.
        model_ide.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

        // To compare with the ODE model we use ExponentialSurvivalFunctions functions as TransitionDistributions.
        // We set the parameters so that they correspond to the above ODE model.
        mio::ExponentialSurvivalFunction exponential(10.0);
        mio::StateAgeFunctionWrapper delaydistribution(exponential);
        std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
        // ExposedToInfectedNoSymptoms
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms].set_distribution_parameter(
            1. / model_ode.parameters.get<mio::osecir::TimeExposed<ScalarType>>()[(mio::AgeGroup)0]);
        // InfectedNoSymptomsToInfectedSymptoms
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
            .set_distribution_parameter(
                1. / model_ode.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0]);
        // InfectedNoSymptomsToRecovered
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]
            .set_distribution_parameter(
                1. / model_ode.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0]);
        // InfectedSymptomsToInfectedSevere
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]
            .set_distribution_parameter(
                1. / model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0]);
        // InfectedSymptomsToRecovered
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered].set_distribution_parameter(
            1. / model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0]);
        // InfectedSevereToInfectedCritical
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]
            .set_distribution_parameter(
                1. / model_ode.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[(mio::AgeGroup)0]);
        // InfectedSevereToRecovered
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered].set_distribution_parameter(
            1. / model_ode.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[(mio::AgeGroup)0]);
        // InfectedCriticalToDead
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead].set_distribution_parameter(
            1. / model_ode.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[(mio::AgeGroup)0]);
        // InfectedCriticalToRecovered
        vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_distribution_parameter(
            1. / model_ode.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[(mio::AgeGroup)0]);

        model_ide.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

        // Set probabilities that determine proportion between compartments.
        std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.);
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
            1 - model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
            model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
            model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)] =
            1 - model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
            model_ode.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)] =
            1 - model_ode.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] =
            model_ode.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0];
        vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)] =
            1 - model_ode.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0];
        model_ide.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

        // Set further parameters.
        mio::ConstantFunction constfunc_proboncontact(
            model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[(mio::AgeGroup)0]);
        mio::StateAgeFunctionWrapper proboncontact(constfunc_proboncontact);
        model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(proboncontact);

        mio::ConstantFunction constfunc_reltransnosympt(
            model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[(mio::AgeGroup)0]);
        mio::StateAgeFunctionWrapper reltransnosympt(constfunc_reltransnosympt);
        model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(reltransnosympt);

        mio::ConstantFunction constfunc_riskofinf(
            model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[(mio::AgeGroup)0]);
        mio::StateAgeFunctionWrapper riskofinf(constfunc_riskofinf);
        model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(riskofinf);

        // Compute initial flows from results of ODE simulation.
        mio::isecir::compute_initial_flows_for_ide_from_ode(model_ode, model_ide, secihurd_ode, t0_ide, dt_ode, dt_ide);

        model_ide.check_constraints(dt_ide);

        // Carry out simulation.
        std::cout << "Starting simulation with IDE model. \n";
        mio::isecir::Simulation sim(model_ide, dt_ide);
        sim.advance(tmax);

        mio::TimeSeries<ScalarType> secihurd_ide       = sim.get_result();
        mio::TimeSeries<ScalarType> secihurd_ide_flows = sim.get_transitions();

        if (print_to_terminal) {
            secihurd_ide.print_table();
            secihurd_ide_flows.print_table();
        }

        std::cout << "Initialization method: " << sim.get_model().get_initialization_method_compartments() << "\n";
        if (save_ide) {
            auto save_result_status_ide =
                mio::save_result({secihurd_ide}, {0}, 1,
                                 result_dir + "result_ide_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                                     "_init_dt_ode=1e-" + fmt::format("{:.0f}", ode_exponent) + ".h5");
            auto save_result_status_ide_flows =
                mio::save_result({remove_time_points(secihurd_ide_flows, dt_ide, 1. / dt_ide)}, {0}, 1,
                                 result_dir + "result_ide_flows_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
                                     "_init_dt_ode=1e-" + fmt::format("{:.0f}", ode_exponent) + ".h5");
            std::cout << "Successfully saved the IDE simulation results. \n";
        }
    }
}
