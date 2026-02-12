/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "ide_secir/model.h"
#include "ide_secir/simulation.h"

#include "boost/filesystem.hpp"
#include <iostream>
#include <string>

/** This file can be used to create simulation results in order to compare ODE and IDE models.
* This means that the parameters of the IDE model are set in such a way that the continuous version of the IDE model
* is reduced to the ODE model used. 
* The IDE model is initialized with flows from the ODE model so that both models are comparable. 
* If these results are generated for different accuracies, the convergence rate of the solution scheme of the IDE model
* can be determined numerically with the ODE model with a high accuracy as ground truth.
*/

// Parameters for the simulation.
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
* @brief Takes the cumulative flows of an ODE simulation and computes the respective flows per time interval, i.e. 
* $\tilde{sigma}.
*
* With t_max and t_window, it can be determined for which time window the flows will be computed. 
* By default, we compute the flows with the same time step size.
* It is also possible to compute the corresponding flows for a bigger time step which can be given by dt_target.
*
* @param[in] model_ode ODE-SECIR model used.
* @param[in] cumulative_flows TimeSeries containing cumulative flows from an ODE simulation. 
*           The time steps must be equidistant.
* @param[out] flows_per_interval TimeSeries where the computed flows per time interval will be stored. 
* @param[in] t_max Maximal time for which the flows are computed.
* @param[in] t_window Time window before t_max for which flows will be computed.
* @param[in] dt_target Time step size used for the resulting TimeSeries flows. 
*       Default is the time step size of the ODE simulation result compartments. 
*       dt_target should be a multiple of the step size used for the simulation result in compartments.
*/
void get_flows_per_time_interval_from_cumulative_flows(mio::osecir::Model<ScalarType>& model_ode,
                                                       mio::TimeSeries<ScalarType> cumulative_flows,
                                                       mio::TimeSeries<ScalarType>& flows_per_interval,
                                                       ScalarType t_max, ScalarType t_window, ScalarType dt_target = 0.)
{
    ScalarType dt_ode = cumulative_flows.get_time(1) - cumulative_flows.get_time(0);

    // Assert that dt_target is a multiple of dt_ode.
    assert(mio::floating_point_equal(std::fmod(dt_target, dt_ode), 0., mio::Limits<ScalarType>::zero_tolerance()));

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;
    // Check that the TimeSeries flows_per_interval is empty as expected.
    if (flows_per_interval.get_num_time_points() > 0) {
        flows_per_interval = mio::TimeSeries<ScalarType>(num_transitions);
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

    // Define flow indices to access (cumulative) flows from ODE simulation.
    std::vector<size_t> flow_indices_ode = {
        model_ode.get_flat_flow_index<mio::osecir::InfectionState::Susceptible, mio::osecir::InfectionState::Exposed>(
            {mio::AgeGroup(0)}),
        model_ode.get_flat_flow_index<mio::osecir::InfectionState::Exposed,
                                      mio::osecir::InfectionState::InfectedNoSymptoms>({mio::AgeGroup(0)}),
        model_ode.get_flat_flow_index<mio::osecir::InfectionState::InfectedNoSymptoms,
                                      mio::osecir::InfectionState::InfectedSymptoms>({mio::AgeGroup(0)}),
        model_ode.get_flat_flow_index<mio::osecir::InfectionState::InfectedNoSymptoms,
                                      mio::osecir::InfectionState::Recovered>({mio::AgeGroup(0)}),
        model_ode.get_flat_flow_index<mio::osecir::InfectionState::InfectedSymptoms,
                                      mio::osecir::InfectionState::InfectedSevere>({mio::AgeGroup(0)}),
        model_ode.get_flat_flow_index<mio::osecir::InfectionState::InfectedSymptoms,
                                      mio::osecir::InfectionState::Recovered>({mio::AgeGroup(0)}),
        model_ode.get_flat_flow_index<mio::osecir::InfectionState::InfectedSevere,
                                      mio::osecir::InfectionState::InfectedCritical>({mio::AgeGroup(0)}),
        model_ode.get_flat_flow_index<mio::osecir::InfectionState::InfectedSevere,
                                      mio::osecir::InfectionState::Recovered>({mio::AgeGroup(0)}),
        model_ode.get_flat_flow_index<mio::osecir::InfectionState::InfectedCritical, mio::osecir::InfectionState::Dead>(
            {mio::AgeGroup(0)}),
        model_ode.get_flat_flow_index<mio::osecir::InfectionState::InfectedCritical,
                                      mio::osecir::InfectionState::Recovered>({mio::AgeGroup(0)})};

    Eigen::Index flows_start_index = t_max_index - t_window_index + 1;

    for (Eigen::Index i = flows_start_index; i <= t_max_index; i++) {
        // Add time point.
        flows_per_interval.add_time_point(i * dt_target,
                                          mio::TimeSeries<ScalarType>::Vector::Constant(num_transitions, 0));

        // Compute flows per time interval from cumulative flows for every transition.
        for (Eigen::Index transition = 0; transition < (Eigen::Index)mio::isecir::InfectionTransition::Count;
             transition++) {
            flows_per_interval.get_last_value()[transition] =
                cumulative_flows[(Eigen::Index)(scale_timesteps * i)][flow_indices_ode[transition]] -
                cumulative_flows[(Eigen::Index)(scale_timesteps * (i - 1))][flow_indices_ode[transition]];
        }
    }
}

/**
* @brief Computes the initial flows (defined per time interval) that are needed for an IDE simulation given simulation 
* results obtained with an ODE model (containing results for both compartments and cumulative flows) for an adequate 
* time window before t0_ide.
* 
* The results of the ODE model can be obtained uing the simulate_flows() function returning both compartments and 
* cumulative flows. 
*
* Here, we assume that the ODE and the IDE model are matching, i.e. that the parameters of the IDE model are chosen 
* such that the continous version reduces to the ODE model. This is achieved by choosing exponentially distributed 
* transitions with the corresponding mean stay times etc. 
* The time step size of ODE and IDE simulation can be chosen independently.
* However, we assume that the time step size of the IDE model is a multiple of the one of the ODE model.
*
* @param[in] model_ode ODE model that is used.
* @param[in] model_ide IDE model that is used.
* @param[in] results_ode Vector of TimeSeries containing results obtained from a simulation with model_ode; first entry 
* contains compartments and second entry contains cumulative flows. 
* @param[in] t0_ide Start time of IDE simulation that we want to compute initial flows for. 
* @param[in] dt_ide Time step size of IDE simulation. 
*/
void compute_initial_flows_for_ide_from_ode(mio::osecir::Model<ScalarType>& model_ode, mio::isecir::Model& model_ide,
                                            std::vector<mio::TimeSeries<ScalarType>> results_ode, ScalarType t0_ide,
                                            ScalarType dt_ide)
{
    mio::TimeSeries<ScalarType> compartments_ode     = results_ode[0];
    mio::TimeSeries<ScalarType> cumulative_flows_ode = results_ode[1];

    // Use t_window=t0_ide to get flows from t0 onwards.
    get_flows_per_time_interval_from_cumulative_flows(model_ode, cumulative_flows_ode, model_ide.transitions, t0_ide,
                                                      t0_ide, dt_ide);
    ScalarType dt_ode = compartments_ode.get_time(1) - compartments_ode.get_time(0);

    // Remove time series from previous run and set initial values in populations.
    if (model_ide.populations.get_num_time_points() > 0) {
        model_ide.populations = mio::TimeSeries<ScalarType>((int)mio::isecir::InfectionState::Count);
    }
    model_ide.populations.add_time_point<Eigen::VectorXd>(
        model_ide.transitions.get_last_time(),
        mio::TimeSeries<ScalarType>::Vector::Constant((int)mio::isecir::InfectionState::Count, 0));
    model_ide.populations[0][Eigen::Index(mio::isecir::InfectionState::Dead)] =
        compartments_ode[(Eigen::Index)compartments_ode.get_num_time_points() -
                         (Eigen::Index)((compartments_ode.get_last_time() - t0_ide) / dt_ode) - 1]
                        [(Eigen::Index)mio::osecir::InfectionState::Dead];
}

/**
* @brief Function to postprocess simulation results before saving. This allows to remove some time points from the 
* simulation results so that not every point has to be saved afterwards. Furthermore, in the case of flows, we can 
* scale the results so that we store the values at time points (\hat{sigma}) instead of flows over a time interval 
* (\tilde{sigma}) as used in the simulation.
*
* @param[in] simulation_result TimeSeries containing simulation results. Can contain compartments or flows.
* @param[in] saving_dt Step size in between the time points of the TimeSeries with less time points.
*   This should be a multiple of the time step size used in simulation_results.
* @param[in] scale_flows Factor by which the TimeSeries values of flows are scaled to obtain flows at time points. 
* @returns TimeSeries with simulation results where some time points have been removed and/or scaled. 
*/
mio::TimeSeries<ScalarType> postprocess_timeseries(const mio::TimeSeries<ScalarType>& simulation_result,
                                                   ScalarType saving_dt, ScalarType scale_flows = 1.)
{
    mio::TimeSeries<ScalarType> removed(simulation_result.get_num_elements());
    ScalarType time = simulation_result.get_time(0);
    for (int i = 0; i < simulation_result.get_num_time_points(); i++) {
        if (mio::floating_point_greater_equal(simulation_result.get_time(i), time, 1e-8)) {
            removed.add_time_point(simulation_result.get_time(i), scale_flows * simulation_result[i]);
            time += saving_dt;
        }
    }
    return removed;
}

mio::UncertainContactMatrix<ScalarType> scale_contact_matrix(ScalarType damping, ScalarType scaling_time)
{

    mio::ContactMatrixGroup contact_matrix = mio::ContactMatrixGroup<ScalarType>(1, 1);
    if (damping <= 1.) {
        // Perform simulation with a decrease in contacts.
        contact_matrix[0] =
            mio::ContactMatrix<ScalarType>(Eigen::MatrixXd::Constant(1, 1, simulation_parameter["cont_freq"]));
        contact_matrix[0].add_damping(0., mio::SimulationTime<ScalarType>(scaling_time));
        contact_matrix[0].add_damping(damping, mio::SimulationTime<ScalarType>(scaling_time + 0.00001));
    }
    // else {
    //     // Perform simulation with an increase in contacts.
    //     contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, damping * cont_freq));
    //     contact_matrix[0].add_damping(1 - 1. / damping, mio::SimulationTime(-100.));
    //     contact_matrix[0].add_damping(1 - 1. / damping, mio::SimulationTime(scaling_time));
    //     contact_matrix[0].add_damping(0., mio::SimulationTime(scaling_time + 0.1));
    // }

    return mio::UncertainContactMatrix(contact_matrix);
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
* run, e.g., to create a ground truth.
* However, we assume that the time step size of the IDE model is a multiple of the one of the ODE model.
*
* @param[in] t0 Start time of the ODE simulation. 
* @param[in] tmax Time up to which we simulate. 
* @param[in] ode_exponent The ODE model is simulated using a fixed step size dt=10^{-ode_exponent}.
* @param[in] ide_exponents The IDE model is simulated using fixed step sizes dt=10^{-ide_exponent} for ide_exponent in 
* ide_exponents.
* @param[in] save_exponent The results of the ODE model will be saved using the step size 10^{-save_exponent}, should 
* not be larger than the maximum ide_exponent.
* @param[in] save_dir Directory where simulation results will be saved. Default is an empty string leading to the 
* results not being saved. 
* @returns Any io errors that happen. 
*/
mio::IOResult<void> simulate_ode_and_ide(ScalarType t0, ScalarType t0_ide, ScalarType tmax, ScalarType ode_exponent,
                                         std::vector<ScalarType> ide_exponents, ScalarType save_exponent,
                                         ScalarType damping, ScalarType scaling_time, std::string save_dir = "")
{
    /**********************************
    *         ODE simulation          *
    **********************************/

    // The ODE model is simulated using a fixed step size dt=10^{-ode_exponent}.
    ScalarType dt_ode = pow(10, -ode_exponent);

    mio::osecir::Model<ScalarType> model_ode(1);

    // Set initial values for compartments. These values are not set realistically as we are considering a synthetic
    //scenario here.

    model_ode.populations.set_total(simulation_parameter["total_population"]);
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Exposed}]                     = 20;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptoms}]          = 20;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedNoSymptomsConfirmed}] = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}]            = 3;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptomsConfirmed}]   = 0;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSevere}]              = 1;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedCritical}]            = 1;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Recovered}]                   = 10;
    model_ode.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::Dead}]                        = 0;
    model_ode.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                    simulation_parameter["total_population"]);

    // Set parameters.
    // If Seasonality=0, then cont_freq_eff is equal to the contact frequency as defined in contact_matrix.
    model_ode.parameters.set<mio::osecir::Seasonality<ScalarType>>(simulation_parameter["Seasonality"]);
    mio::ContactMatrixGroup<ScalarType>& contact_matrix =
        model_ode.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix[0] =
        mio::ContactMatrix<ScalarType>(Eigen::MatrixXd::Constant(1, 1, simulation_parameter["cont_freq"]));
    model_ode.parameters.get<mio::osecir::ContactPatterns<ScalarType>>() = scale_contact_matrix(damping, scaling_time);
    // mio::UncertainContactMatrix(contact_matrix);

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
        std::make_unique<mio::ControlledStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta_cash_karp54>>();
    // Choose dt_min = dt_max so that we have a fixed time step and can compare to IDE.
    integrator->set_dt_min(dt_ode);
    integrator->set_dt_max(dt_ode);
    // Set tolerance as follows so that every time step is only computed once (found out by trying).
    integrator->set_rel_tolerance(1e-1);
    integrator->set_abs_tolerance(1e-1);

    std::cout << "Starting simulation with ODE model. \n";

    // Vector that contains ODE simulation results. First entry contains values for compartments, second entry contains
    // values for cumulative flows.
    std::vector<mio::TimeSeries<ScalarType>> results_ode =
        mio::osecir::simulate_flows<ScalarType>(t0, tmax, dt_ode, model_ode, std::move(integrator));

    if (!save_dir.empty() && save_exponent > 0) {
        // Create result directory if not existent yet.
        boost::filesystem::path res_dir(save_dir);
        boost::filesystem::create_directory(res_dir);

        // Save compartments.
        auto save_result_status_ode = mio::save_result(
            {postprocess_timeseries(results_ode[0], pow(10, -save_exponent))}, {0}, 1, save_dir + "result_ode.h5");

        // Compute flows per time interval from cumulative flows of ODE simulation, i.e. we are computing \tilde{\sigma}
        // here.
        mio::TimeSeries<ScalarType> flows_per_timestep_ode((int)mio::isecir::InfectionTransition::Count);
        get_flows_per_time_interval_from_cumulative_flows(model_ode, results_ode[1], flows_per_timestep_ode, tmax,
                                                          tmax - t0);
        // Save flows.
        // To be able to compare flows between different time step sizes (of ODE and IDE) we need to divide by dt to get
        // \hat{\sigma}. This is done while saving the results.
        auto save_result_status_ode_flows =
            mio::save_result({postprocess_timeseries(flows_per_timestep_ode, pow(10, -save_exponent), 1. / dt_ode)},
                             {0}, 1, save_dir + "result_ode_flows.h5");
    }

    /**********************************
    *         IDE simulation          *
    **********************************/
    if (!ide_exponents.empty()) {
        // Set up IDE model.

        // Number of deaths will be set according to the ODE model later in the function where also the transitions are calculated.
        ScalarType deaths_init_value = 0.;

        // Initialize model.
        mio::TimeSeries<ScalarType> init_transitions((int)mio::isecir::InfectionTransition::Count);
        size_t num_agegroups = 1;
        mio::CustomIndexArray<ScalarType, mio::AgeGroup> total_population =
            mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups),
                                                             simulation_parameter["total_population"]);
        mio::CustomIndexArray<ScalarType, mio::AgeGroup> deaths =
            mio::CustomIndexArray<ScalarType, mio::AgeGroup>(mio::AgeGroup(num_agegroups), deaths_init_value);
        mio::isecir::Model model_ide(std::move(init_transitions), total_population, deaths, num_agegroups);

        // Set parameters.
        // contact_matrix was already defined for ODE.
        model_ide.parameters.get<mio::isecir::ContactPatterns>() = scale_contact_matrix(damping, scaling_time);
        //mio::UncertainContactMatrix(contact_matrix);

        // To compare with the ODE model we use ExponentialSurvivalFunctions functions as TransitionDistributions.
        // We set the parameters so that they correspond to the above ODE model.
        mio::ExponentialSurvivalFunction exponential(10.0);
        mio::StateAgeFunctionWrapper<ScalarType> delaydistribution(exponential);
        std::vector<mio::StateAgeFunctionWrapper<ScalarType>> vec_delaydistrib(
            (int)mio::isecir::InfectionTransition::Count, delaydistribution);
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

        // Set probabilities.
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
        mio::ConstantFunction<ScalarType> constfunc_proboncontact(
            model_ode.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>()[(mio::AgeGroup)0]);
        mio::StateAgeFunctionWrapper<ScalarType> proboncontact(constfunc_proboncontact);
        model_ide.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(proboncontact);

        mio::ConstantFunction<ScalarType> constfunc_reltransnosympt(
            model_ode.parameters.get<mio::osecir::RelativeTransmissionNoSymptoms<ScalarType>>()[(mio::AgeGroup)0]);
        mio::StateAgeFunctionWrapper<ScalarType> reltransnosympt(constfunc_reltransnosympt);
        model_ide.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(reltransnosympt);

        mio::ConstantFunction<ScalarType> constfunc_riskofinf(
            model_ode.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()[(mio::AgeGroup)0]);
        mio::StateAgeFunctionWrapper<ScalarType> riskofinf(constfunc_riskofinf);
        model_ide.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(riskofinf);

        // Compute initial flows and simulate for each ide_exponent.
        for (ScalarType ide_exponent : ide_exponents) {

            // The IDE model is simulated using a fixed step size dt=10^{-ide_exponent}.
            ScalarType dt_ide = pow(10, -ide_exponent);

            // Compute initial flows from results of ODE simulation and set initial values for populations.
            compute_initial_flows_for_ide_from_ode(model_ode, model_ide, results_ode, t0_ide, dt_ide);

            model_ide.check_constraints(dt_ide);

            // Carry out simulation.
            std::cout << "Starting simulation with IDE model. \n";
            mio::isecir::Simulation sim(model_ide, dt_ide);
            sim.advance(tmax);

            std::cout << "Initialization method of the IDE model: "
                      << sim.get_model().get_initialization_method_compartments() << "\n";

            if (!save_dir.empty()) {
                // Save compartments.
                mio::TimeSeries<ScalarType> secihurd_ide = sim.get_result();
                auto save_result_status_ide = mio::save_result({secihurd_ide}, {0}, 1, save_dir + "result_ide.h5");
                // Save flows.
                // To be able to compare flows between different timesteps (of ODE and IDE) we need to divide by dt to get
                // \hat{\sigma}. This is done while saving the results.
                mio::TimeSeries<ScalarType> secihurd_ide_flows = sim.get_transitions();
                auto save_result_status_ide_flows =
                    mio::save_result({postprocess_timeseries(secihurd_ide_flows, dt_ide, 1. / dt_ide)}, {0}, 1,
                                     save_dir + "result_ide_flows.h5");
            }
        }
    }
    return mio::success();
}

/** 
* Usage: ide_convergence_rate <result_dir>
* The command line argument is optional. Default values are provided if not specified.
*/
int main(int argc, char** argv)
{
    // Default path is valid if script is executed in e.g. memilio-simulations/2024_Wendler_et_al_Nonstandard_numerical_scheme_IDE.
    std::string result_dir = "../../simulation_results/2026-02-03/secir_changepoints/";

    // Set result_dir via command line.
    if (argc == 2) {
        result_dir = argv[1];
    }

    // Make folder if not existent yet.
    boost::filesystem::path dir(result_dir);
    boost::filesystem::create_directories(dir);

    // General set up.
    ScalarType t0     = 0.;
    ScalarType t0_ide = 35.;
    ScalarType tmax   = 60.;

    ScalarType damping      = 0.7;
    ScalarType scaling_time = 40.;

    // The ODE model will be simulated using a fixed step size dt=10^{-ode_exponent}.
    ScalarType ode_exponent = 2;
    // The results of the ODE model will be saved using the step size 10^{-save_exponent}
    // as for very small step sizes used for the simulation, the number of time points stored gets very big.
    ScalarType save_exponent = 2;
    // The IDE model will be simulated using a fixed step size dt=10^{-ide_exponent} for ide_exponent in ide_exponents.
    std::vector<ScalarType> ide_exponents = {2};

    mio::IOResult<void> result = simulate_ode_and_ide(t0, t0_ide, tmax, ode_exponent, ide_exponents, save_exponent,
                                                      damping, scaling_time, result_dir);

    if (!result) {
        printf("%s\n", result.error().formatted_message().c_str());
        return -1;
    }

    return 0;
}