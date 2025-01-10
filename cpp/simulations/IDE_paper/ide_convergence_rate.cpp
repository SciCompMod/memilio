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
#include "memilio/io/result_io.h"
#include "memilio/utils/time_series.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/floating_point.h"

#include "ode_secir/infection_state.h"
#include "ode_secir/model.h"

#include "ide_secir/infection_state.h"
// #include "ide_secir/model.h"
// #include "ide_secir/simulation.h"

#include <iostream>
#include <string>

/** This file can be used to create simulation results in order to compare ODE and IDE models.
* This means that the parameters of the IDE model are set in such a way that the continuous version of the IDE model
* is reduced to the ODE model used. 
* The IDE model is initialized with flows from the ODE model so that both models are comparable. 
* If these results are generated for different accuracies, the convergence rate of the solution scheme of the IDE model
* can be determined numerically with the ODE model with a high accuracy as ground truth.
*/

// Used parameters.
std::map<std::string, ScalarType> simulation_parameter = {{"t0", 0.},
                                                          {"total_population", 10000.},
                                                          {"TimeExposed", 1.4},
                                                          {"TimeInfectedNoSymptoms", 1.2},
                                                          {"TimeInfectedSymptoms", 0.3},
                                                          {"TimeInfectedSevere", 0.3},
                                                          {"TimeInfectedCritical", 0.3},
                                                          {"TransmissionProbabilityOnContact", 1.},
                                                          {"RelativeTransmissionNoSymptoms", 1.},
                                                          {"RiskOfInfectionFromSymptomatic", 1.},
                                                          {"Seasonality", 0.},
                                                          {"InfectedSymptomsPerInfectedNoSymptoms", 0.5},
                                                          {"SeverePerInfectedSymptoms", 0.5},
                                                          {"CriticalPerSevere", 0.5},
                                                          {"DeathsPerCritical", 0.5},
                                                          {"cont_freq", 1.}};

/**
* @brief Takes the compartments of an ODE simulation and computes the respective flows. 
*
* With t_max and t_window, it can be determined for which time window the flows will be computed. 
* By default, we compute the flows with the same time step size.
* It is also possible to compute the corresponding flows for a bigger time step which can be given by dt_target.
*
* @param[in] model_ode ODE-SECIR model used.
* @param[in] compartments TimeSeries containing compartments from an ODE simulation. 
*           The time steps must be equidistant.
* @param[out] flows TimeSeries where the computed flows will be stored. 
* @param[in] t_max Maximal time for which the flows are computed.
* @param[in] t_window Time window before t_max for which flows will be computed.
* @param[in] dt_target Time step size used for the resulting TimeSeries flows. 
*       Default is the time step size of the ODE simulation result compartments. 
*       dt_target should be a multiple of the step size used for the simulation result in compartments.
*/
void get_flows_from_ode_compartments(mio::osecir::Model<ScalarType>& model_ode,
                                     mio::TimeSeries<ScalarType> compartments, mio::TimeSeries<ScalarType>& flows,
                                     ScalarType t_max, ScalarType t_window, ScalarType dt_target = 0.)
{
    ScalarType dt_ode   = compartments.get_time(1) - compartments.get_time(0);
    int num_transitions = (int)mio::isecir::InfectionTransition::Count;
    // Check that the TimeSeries flows is empty as expected.
    if (flows.get_num_time_points() > 0) {
        flows = mio::TimeSeries<ScalarType>(num_transitions);
    }
    // If dt_target is not set, use dt_ode.
    if (dt_target < 1e-10) {
        dt_target = dt_ode;
    }
    // scale_timesteps is used to get from index wrt ODE timestep to index wrt dt_target.
    // Here we assume that the ODE model is solved on a finer (or equal) scale than dt_target.
    ScalarType scale_timesteps = dt_target / dt_ode;

    // Compute index variables with respect to dt_target.
    Eigen::Index t_window_index = Eigen::Index(std::ceil(t_window / dt_target));
    Eigen::Index t_max_index    = Eigen::Index(std::ceil(t_max / dt_target));

    Eigen::Index flows_start_index = t_max_index - t_window_index + 1;

    // Add time points to TimeSeries flows and set flow Susceptible to Exposed.
    for (Eigen::Index i = flows_start_index; i <= t_max_index; i++) {
        flows.add_time_point(i * dt_target, mio::TimeSeries<ScalarType>::Vector::Constant(num_transitions, 0));
        flows.get_last_value()[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)] +=
            compartments[(Eigen::Index)(scale_timesteps * (i - 1))]
                        [(Eigen::Index)mio::osecir::InfectionState::Susceptible] -
            compartments[(Eigen::Index)(scale_timesteps * i)][(Eigen::Index)mio::osecir::InfectionState::Susceptible];
    }

    // --- Compute flows as combination of change in compartments and previously computed flows.
    // Flow from Exposed to InfectedNoSymptoms.
    for (Eigen::Index i = flows_start_index; i <= t_max_index; i++) {
        flows[i - flows_start_index][(Eigen::Index)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms] =
            compartments[(Eigen::Index)(scale_timesteps * (i - 1))]
                        [(Eigen::Index)mio::osecir::InfectionState::Exposed] -
            compartments[(Eigen::Index)(scale_timesteps * i)][(Eigen::Index)mio::osecir::InfectionState::Exposed] +
            flows[i - flows_start_index][(Eigen::Index)mio::isecir::InfectionTransition::SusceptibleToExposed];
    }
    ScalarType out_flow = 0;
    // Flow from InfectedNoSymptoms to InfectedSymptoms and from InfectedNoSymptoms to Recovered.
    for (Eigen::Index i = flows_start_index; i <= t_max_index; i++) {
        out_flow =
            compartments[(Eigen::Index)(scale_timesteps * (i - 1))]
                        [(Eigen::Index)mio::osecir::InfectionState::InfectedNoSymptoms] -
            compartments[(Eigen::Index)(scale_timesteps * i)]
                        [(Eigen::Index)mio::osecir::InfectionState::InfectedNoSymptoms] +
            flows[i - flows_start_index][(Eigen::Index)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms];
        flows[i -
              flows_start_index][(Eigen::Index)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] =
            (1 -
             model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0]) *
            out_flow;
        flows[i - flows_start_index][(Eigen::Index)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered] =
            model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] *
            out_flow;
    }
    // Flow from InfectedSymptoms to InfectedSevere and from InfectedSymptoms to Recovered.
    for (Eigen::Index i = flows_start_index; i <= t_max_index; i++) {
        out_flow = compartments[(Eigen::Index)(scale_timesteps * (i - 1))]
                               [(Eigen::Index)mio::osecir::InfectionState::InfectedSymptoms] -
                   compartments[(Eigen::Index)(scale_timesteps * i)]
                               [(Eigen::Index)mio::osecir::InfectionState::InfectedSymptoms] +
                   flows[i - flows_start_index]
                        [(Eigen::Index)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms];
        flows[i - flows_start_index][(Eigen::Index)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere] =
            model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0] * out_flow;
        flows[i - flows_start_index][(Eigen::Index)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered] =
            (1 - model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0]) *
            out_flow;
    }
    // Flow from InfectedSevere to InfectedCritical and from InfectedSevere to Recovered.
    for (Eigen::Index i = flows_start_index; i <= t_max_index; i++) {
        out_flow = compartments[(Eigen::Index)(scale_timesteps * (i - 1))]
                               [(Eigen::Index)mio::osecir::InfectionState::InfectedSevere] -
                   compartments[(Eigen::Index)(scale_timesteps * i)]
                               [(Eigen::Index)mio::osecir::InfectionState::InfectedSevere] +
                   flows[i - flows_start_index]
                        [(Eigen::Index)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere];
        flows[i - flows_start_index][(Eigen::Index)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical] =
            model_ode.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0] * out_flow;
        flows[i - flows_start_index][(Eigen::Index)mio::isecir::InfectionTransition::InfectedSevereToRecovered] =
            (1 - model_ode.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0]) * out_flow;
    }
    // Flow from InfectedCritical to Dead and from InfectedCritical to Recovered.
    for (Eigen::Index i = flows_start_index; i <= t_max_index; i++) {
        out_flow = compartments[(Eigen::Index)(scale_timesteps * (i - 1))]
                               [(Eigen::Index)mio::osecir::InfectionState::InfectedCritical] -
                   compartments[(Eigen::Index)(scale_timesteps * i)]
                               [(Eigen::Index)mio::osecir::InfectionState::InfectedCritical] +
                   flows[i - flows_start_index]
                        [(Eigen::Index)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical];
        flows[i - flows_start_index][(Eigen::Index)mio::isecir::InfectionTransition::InfectedCriticalToDead] =
            model_ode.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0] * out_flow;
        flows[i - flows_start_index][(Eigen::Index)mio::isecir::InfectionTransition::InfectedCriticalToRecovered] =
            (1 - model_ode.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0]) * out_flow;
    }
}

/**
* @brief Computes the inital flows that are needed for an IDE simulation given we have the compartments from an ODE 
* simulation for an adequate time window before t0_ide. 
*
* Here we assume, that the ODE and the IDE model are matching, i.e. that the parameters of the IDE model are chosen 
* such that the continous version reduces to the ODE model. This is achieved by choosing exponentially distributed 
* transitions with the corresponding mean stay times etc. 
* The time step size of ODE and IDE simulation can be chosen independently.
* However, we assume that the time step size of the IDE model is a multiple of the one of the ODE model.
*
* @param[in] model_ode ODE model that is used.
* @param[in] model_ide IDE model that is used.
* @param[in] compartments TimeSeries containing compartments from a simulation with model_ode. 
* @param[in] t0_ide Start time of IDE simulation that we want to compute initial flows for. 
* @param[in] dt_ide Time step size of IDE simulation. 
*/
// void compute_initial_flows_for_ide_from_ode(mio::osecir::Model<ScalarType>& model_ode, mio::isecir::Model& model_ide,
//                                             mio::TimeSeries<ScalarType> compartments, ScalarType t0_ide,
//                                             ScalarType dt_ide)
// {
//     std::cout << "Computing initial flows. \n";

//     // Use t_window=t0_ide to get flows from t0 onwards.
//     get_flows_from_ode_compartments(model_ode, compartments, model_ide.m_transitions, t0_ide, t0_ide, dt_ide);
//     ScalarType dt_ode = compartments.get_time(1) - compartments.get_time(0);
//     // Remove time series from previous run and set initial values in populations.
//     if (model_ide.m_populations.get_num_time_points() > 0) {
//         model_ide.m_populations = mio::TimeSeries<ScalarType>((int)mio::isecir::InfectionState::Count);
//     }
//     model_ide.m_populations.add_time_point<Eigen::VectorXd>(
//         model_ide.m_transitions.get_last_time(),
//         mio::TimeSeries<ScalarType>::Vector::Constant((int)mio::isecir::InfectionState::Count, 0));
//     model_ide.m_populations[0][Eigen::Index(mio::isecir::InfectionState::Dead)] =
//         compartments[(Eigen::Index)compartments.get_num_time_points() -
//                      (Eigen::Index)((compartments.get_last_time() - t0_ide) / dt_ode) - 1]
//                     [(Eigen::Index)mio::osecir::InfectionState::Dead];
// }

/**
* @brief Function to remove time points from some simulation results so that not every point has to be saved afterwards.
*
* @param[in] simulation_result TimeSeries containing simulation results. Can contain compartments or flows.
* @param[in] saving_dt Step size in between the time points of the TimeSeries with less time points.
*   This should be a multiple of the time step size used in simulation_results.
* @param[in] scale Factor by which the TimeSeries values should be scaled.
* @returns TimeSeries with simulation results where some time points have been removed. 
*/
mio::TimeSeries<ScalarType> remove_time_points(const mio::TimeSeries<ScalarType>& simulation_result,
                                               ScalarType saving_dt, ScalarType scale = 1.)
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

/**
* @brief Simulates from t0 until tmax using an ODE model and subsequently simulates from 
* t0_ide = (tmax-t0)/2 until tmax using a matching IDE model to determine convergence of the IDE solver. 
*
* To make both ODE and IDE model comparable, we set the parameters of the IDE model according to the parameters of the 
* ODE model. Furthermore, we determine the initial flows for the IDE model based on the ODE results so that we have
* equivalent conditions for both models at t0_ide. 

* The time step size of ODE and IDE simulation can be chosen independently. A vector containing the desired 
* ide_exponents is passed for which an IDE simulation is run. If an empty vector is passed, only an ODE simulation is 
* run, e.g., to create a ground truth
* However, we assume that the time step size of the IDE model is a multiple of the one of the ODE model.
*
* @param[in] t0 Start time of the ODE simulation. 
* @param[in] tmax Maximal time for which we simulate.
* @param[in] ode_exponent The ODE model is simulated using a fixed step size dt=10^{-ode_exponent}.
* @param[in] ide_exponents The IDE model is simulated using fixed step sizes dt=10^{-ide_exponent} for ide_exponent in 
* ide_exponents.
* @param[in] save_exponent The results of the ODE model will be saved using the step size 10^{-save_exponent}, should 
* not be larger than the maximum ide_exponent.
* @param[in] result_dir Directory where simulation results will be stored. 
* @returns Any io errors that happen. 
*/
mio::IOResult<void> simulate_ode_and_ide(ScalarType t0, ScalarType tmax, ScalarType ode_exponent,
                                         ScalarType ide_exponent, ScalarType save_exponent, std::string result_dir)
{
    mio::unused(ide_exponent);
    /**********************************
    *         ODE simulation          *
    **********************************/

    // The ODE model is simulated using a fixed step size dt=10^{-ode_exponent}.
    ScalarType dt_ode = pow(10, -ode_exponent);

    mio::osecir::Model<ScalarType> model_ode(1);

    // Set initial values for compartments.
    ScalarType nb_exp_t0 = 20, nb_car_t0 = 20, nb_inf_t0 = 3, nb_hosp_t0 = 1, nb_icu_t0 = 1, nb_rec_t0 = 10,
               nb_dead_t0 = 0;

    model_ode.populations.set_total(simulation_parameter["total_population"]);
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
                                                    simulation_parameter["total_population"]);

    // Set parameters.
    ScalarType cont_freq = simulation_parameter["cont_freq"];
    // Set Seasonality=0 so that cont_freq_eff is equal to contact_matrix.
    model_ode.parameters.set<mio::osecir::Seasonality<ScalarType>>(simulation_parameter["Seasonality"]);
    mio::ContactMatrixGroup& contact_matrix = model_ode.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    model_ode.parameters.get<mio::osecir::ContactPatterns<ScalarType>>() = mio::UncertainContactMatrix(contact_matrix);

    // Parameters needed to determine transition rates.
    model_ode.parameters.get<mio::osecir::TimeExposed<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeExposed"];
    model_ode.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeInfectedNoSymptoms"];
    model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeInfectedSymptoms"];
    model_ode.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeInfectedSevere"];
    model_ode.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TimeInfectedCritical"];

    // Set probabilities that determine proportion between compartments.
    model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        1 - simulation_parameter["InfectedSymptomsPerInfectedNoSymptoms"];
    model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["SeverePerInfectedSymptoms"];
    model_ode.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["CriticalPerSevere"];
    model_ode.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["DeathsPerCritical"];

    // Further model parameters.
    model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["TransmissionProbabilityOnContact"];
    model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["RelativeTransmissionNoSymptoms"];
    model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[(mio::AgeGroup)0] =
        simulation_parameter["RiskOfInfectionFromSymptomatic"];
    // Choose TestAndTraceCapacity very large so that parameter has no effect.
    model_ode.parameters.get<mio::osecir::TestAndTraceCapacity<ScalarType>>() = std::numeric_limits<ScalarType>::max();
    // Choose ICUCapacity very large so that parameter has no effect.
    model_ode.parameters.get<mio::osecir::ICUCapacity<ScalarType>>() = std::numeric_limits<ScalarType>::max();

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

    if (!result_dir.empty() && save_exponent > 0) {
        // Create result directory if not existent yet.
        boost::filesystem::path res_dir(result_dir);
        boost::filesystem::create_directory(res_dir);
        auto save_result_status_ode =
            mio::save_result({remove_time_points(secihurd_ode, pow(10, -save_exponent))}, {0}, 1,
                             result_dir + "result_ode_dt=1e-" + fmt::format("{:.0f}", ode_exponent) + "_savefrequency" +
                                 fmt::format("{:.0f}", save_exponent) + ".h5");

        // Compute flows from ODE result to store results.
        // Note that we are computing \tilde{\sigma} here. To be able to compare flows between different timesteps (of ODE and IDE)
        // we need to divide by dt to get \hat{\sigma}. This is done while saving the results.
        mio::TimeSeries<ScalarType> secihurd_ode_flows((int)mio::isecir::InfectionTransition::Count);
        get_flows_from_ode_compartments(model_ode, secihurd_ode, secihurd_ode_flows, tmax, tmax - t0);
        auto save_result_status_ode_flows =
            mio::save_result({remove_time_points(secihurd_ode_flows, pow(10, -save_exponent), 1. / dt_ode)}, {0}, 1,
                             result_dir + "result_ode_flows_dt=1e-" + fmt::format("{:.0f}", ode_exponent) +
                                 "_savefrequency" + fmt::format("{:.0f}", save_exponent) + ".h5");

        if (save_result_status_ode && save_result_status_ode_flows) {
            std::cout << "Successfully saved the ODE simulation results. \n\n";
        }
        else {
            return mio::failure(mio::StatusCode::InvalidValue,
                                "Error occured while saving the ODE simulation results.");
        }
    }

    // /**********************************
    // *         IDE simulation          *
    // **********************************/

    // // Start IDE model simulation at half of tmax.
    // ScalarType t0_ide = (tmax - t0) / 2.;
    // // Number of deaths will be set according to the ODE model later in the function where also the transitions are calculated.
    // ScalarType deaths_init_value = 0.;

    // // Initialize model.
    // mio::TimeSeries<ScalarType> init_transitions((int)mio::isecir::InfectionTransition::Count);
    // size_t num_agegroups = 1;
    // mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_population =
    //     mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups),
    //                                                      simulation_parameter["total_population"]);
    // mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths =
    //     mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), deaths_init_value);
    // mio::isecir::Model model_ide(std::move(init_transitions), total_population, deaths, num_agegroups);

    // // Set parameters.
    // // Contact matrix; contact_matrix was already defined for ODE.
    // model_ide.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    // // To compare with the ODE model we use ExponentialSurvivalFunctions functions as TransitionDistributions.
    // // We set the parameters so that they correspond to the above ODE model.
    // mio::ExponentialSurvivalFunction exponential(10.0);
    // mio::StateAgeFunctionWrapper delaydistribution(exponential);
    // std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib((int)mio::isecir::InfectionTransition::Count,
    //                                                            delaydistribution);
    // // ExposedToInfectedNoSymptoms
    // vec_delaydistrib[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms].set_distribution_parameter(
    //     1. / model_ode.parameters.get<mio::osecir::TimeExposed<ScalarType>>()[(mio::AgeGroup)0]);
    // // InfectedNoSymptomsToInfectedSymptoms
    // vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
    //     .set_distribution_parameter(
    //         1. / model_ode.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0]);
    // // InfectedNoSymptomsToRecovered
    // vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered].set_distribution_parameter(
    //     1. / model_ode.parameters.get<mio::osecir::TimeInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0]);
    // // InfectedSymptomsToInfectedSevere
    // vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]
    //     .set_distribution_parameter(
    //         1. / model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0]);
    // // InfectedSymptomsToRecovered
    // vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered].set_distribution_parameter(
    //     1. / model_ode.parameters.get<mio::osecir::TimeInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0]);
    // // InfectedSevereToInfectedCritical
    // vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]
    //     .set_distribution_parameter(
    //         1. / model_ode.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[(mio::AgeGroup)0]);
    // // InfectedSevereToRecovered
    // vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered].set_distribution_parameter(
    //     1. / model_ode.parameters.get<mio::osecir::TimeInfectedSevere<ScalarType>>()[(mio::AgeGroup)0]);
    // // InfectedCriticalToDead
    // vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead].set_distribution_parameter(
    //     1. / model_ode.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[(mio::AgeGroup)0]);
    // // InfectedCriticalToRecovered
    // vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered].set_distribution_parameter(
    //     1. / model_ode.parameters.get<mio::osecir::TimeInfectedCritical<ScalarType>>()[(mio::AgeGroup)0]);

    // model_ide.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    // // Set probabilities.
    // std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 0.);
    // vec_prob[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)]        = 1;
    // vec_prob[Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] = 1;
    // vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
    //     1 - model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0];
    // vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
    //     model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms<ScalarType>>()[(mio::AgeGroup)0];
    // vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
    //     model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0];
    // vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)] =
    //     1 - model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms<ScalarType>>()[(mio::AgeGroup)0];
    // vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
    //     model_ode.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0];
    // vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)] =
    //     1 - model_ode.parameters.get<mio::osecir::CriticalPerSevere<ScalarType>>()[(mio::AgeGroup)0];
    // vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] =
    //     model_ode.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0];
    // vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)] =
    //     1 - model_ode.parameters.get<mio::osecir::DeathsPerCritical<ScalarType>>()[(mio::AgeGroup)0];

    // model_ide.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    // // Set further parameters.
    // mio::ConstantFunction constfunc_proboncontact(
    //     model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[(mio::AgeGroup)0]);
    // mio::StateAgeFunctionWrapper proboncontact(constfunc_proboncontact);
    // model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(proboncontact);

    // mio::ConstantFunction constfunc_reltransnosympt(
    //     model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[(mio::AgeGroup)0]);
    // mio::StateAgeFunctionWrapper reltransnosympt(constfunc_reltransnosympt);
    // model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(reltransnosympt);

    // mio::ConstantFunction constfunc_riskofinf(
    //     model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[(mio::AgeGroup)0]);
    // mio::StateAgeFunctionWrapper riskofinf(constfunc_riskofinf);
    // model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(riskofinf);

    // // for (ScalarType ide_exponent : ide_exponents) {

    // // The IDE model is simulated using a fixed step size dt=10^{-ide_exponent}.
    // ScalarType dt_ide = pow(10, -ide_exponent);

    // // Compute initial flows from results of ODE simulation and set initial values for populations.
    // compute_initial_flows_for_ide_from_ode(model_ode, model_ide, secihurd_ode, t0_ide, dt_ide);

    // model_ide.check_constraints(dt_ide);

    // // Carry out simulation.
    // std::cout << "Starting simulation with IDE model. \n";
    // mio::isecir::Simulation sim(model_ide, dt_ide);
    // sim.advance(tmax);

    // std::cout << "Initialization method of the IDE model: " << sim.get_model().get_initialization_method_compartments()
    //           << "\n";
    // if (!result_dir.empty()) {
    //     // Save compartments.
    //     mio::TimeSeries<ScalarType> secihurd_ide = sim.get_result();
    //     auto save_result_status_ide =
    //         mio::save_result({secihurd_ide}, {0}, 1,
    //                          result_dir + "result_ide_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
    //                              "_init_dt_ode=1e-" + fmt::format("{:.0f}", ode_exponent) + ".h5");
    //     // Save flows.
    //     mio::TimeSeries<ScalarType> secihurd_ide_flows = sim.get_transitions();
    //     auto save_result_status_ide_flows =
    //         mio::save_result({remove_time_points(secihurd_ide_flows, dt_ide, 1. / dt_ide)}, {0}, 1,
    //                          result_dir + "result_ide_flows_dt=1e-" + fmt::format("{:.0f}", ide_exponent) +
    //                              "_init_dt_ode=1e-" + fmt::format("{:.0f}", ode_exponent) + ".h5");
    //     if (save_result_status_ide && save_result_status_ide_flows) {
    //         std::cout << "Successfully saved the IDE simulation results. \n\n";
    //     }
    //     else {
    //         std::cout << "Error occured while saving the IDE simulation results. \n";
    //         return mio::failure(mio::StatusCode::InvalidValue,
    //                             "Error occured while saving the IDE simulation results.");
    //     }
    // }
    return mio::success();
}

int main()
{
    // Directory where results will be stored. If this string is empty, results will not be saved.
    std::string result_dir = "../../data/simulation_results/convergence/";

    // General set up.
    ScalarType t0   = 0.;
    ScalarType tmax = 70.;
    // The ODE model will be simulated using a fixed step size dt=10^{-ode_exponent}.
    ScalarType ode_exponent = 6;
    // The results of the ODE model will be saved using the step size 10^{-save_exponent}
    // as for very small step sizes used for the simulation, the number of time points stored gets very big.
    ScalarType save_exponent = 4;
    // The IDE model will be simulated using a fixed step size dt=10^{-ide_exponent} for ide_exponent in ide_exponents.
    std::vector<ScalarType> ide_exponents = {1, 2, 3, 4};

    for (ScalarType ide_exponent : ide_exponents) {
        mio::IOResult<void> result =
            simulate_ode_and_ide(t0, tmax, ode_exponent, ide_exponent, save_exponent, result_dir);

        if (!result) {
            printf("%s\n", result.error().formatted_message().c_str());
            return -1;
        }
    }

    return 0;
}
