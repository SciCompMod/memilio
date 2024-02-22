/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Anna Wendler, Lena Ploetzke, Martin J. Kuehn
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
#include "ide_secir/model_ide.h"
#include "ide_secir/parameters.h"
#include "ode_secir/model.h"
#include "infection_state.h"
#include "memilio/config.h"
#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"

namespace mio
{
namespace isecir
{

Model::Model(TimeSeries<ScalarType>&& init, ScalarType N_init, ScalarType Dead_before,
             const ParameterSet& Parameterset_init)
    : parameters{Parameterset_init}
    , m_transitions{std::move(init)}
    , m_populations{TimeSeries<ScalarType>(Eigen::Index(InfectionState::Count))}
    , m_N{N_init}
    , m_deaths_before{Dead_before}
{
    // // add first timepoint to m_populations at last time from m_transitions
    // m_populations.add_time_point<Eigen::VectorXd>(
    //     m_transitions.get_last_time(), TimeSeries<ScalarType>::Vector::Constant((int)InfectionState::Count, 0));
}

/********************************************************
* Compute initial flows based on given compartment data *
********************************************************/

void Model::compute_initial_flows_from_ode_compartments(mio::osecir::Model model_ode,
                                                        mio::TimeSeries<ScalarType> secihurd_ode, ScalarType t0_ide,
                                                        ScalarType dt)
{
    std::cout << "Computing initial flows. \n";
    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // get (global) support_max to determine how many flows in the past we have to compute
    ScalarType global_support_max         = this->get_global_support_max(dt);
    Eigen::Index global_support_max_index = std::ceil(global_support_max / dt);

    // remove time point
    this->m_transitions.remove_last_time_point();

    ScalarType t0_ide_index = std::ceil(t0_ide / dt);
    unused(secihurd_ode);
    int init_start_index = t0_ide_index - global_support_max_index + 1;
    // flow from S to E for -6*global_support_max, ..., 0 (directly from compartments)
    // add time points to init_transitions here
    for (int i = init_start_index; i <= t0_ide_index; i++) {
        this->m_transitions.add_time_point(i * dt, mio::TimeSeries<ScalarType>::Vector::Constant(num_transitions, 0));
        this->m_transitions.get_last_value()[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)] =
            secihurd_ode[i - 1][Eigen::Index(mio::isecir::InfectionState::Susceptible)] -
            secihurd_ode[i][Eigen::Index(mio::isecir::InfectionState::Susceptible)];
    }

    // compute resulting flows as combination of change in compartments and previously computed flows
    // Eigen::Index start_shift = t0_ide_index - 6 * global_support_max_index;

    // flow from E to C for -global_support_max, ..., 0
    for (int i = init_start_index; i <= t0_ide_index; i++) {
        this->m_transitions[i - init_start_index]
                           [Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)] =
            secihurd_ode[i - 1][Eigen::Index(mio::isecir::InfectionState::Exposed)] -
            secihurd_ode[i][Eigen::Index(mio::isecir::InfectionState::Exposed)] +
            this->m_transitions[i - init_start_index]
                               [Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)];
    }

    // flow from C to I and from C to R for -global_support_max, ..., 0
    for (int i = init_start_index; i <= t0_ide_index; i++) {
        this->m_transitions[i - init_start_index]
                           [Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] =
            (1 - model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]) *
            (secihurd_ode[i - 1][Eigen::Index(mio::isecir::InfectionState::InfectedNoSymptoms)] -
             secihurd_ode[i][Eigen::Index(mio::isecir::InfectionState::InfectedNoSymptoms)] +
             this->m_transitions[i - init_start_index]
                                [Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)]);
        this->m_transitions[i - init_start_index]
                           [Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered)] =
            (model_ode.parameters.get<mio::osecir::RecoveredPerInfectedNoSymptoms>()[(mio::AgeGroup)0]) *
            (secihurd_ode[i - 1][Eigen::Index(mio::isecir::InfectionState::InfectedNoSymptoms)] -
             secihurd_ode[i][Eigen::Index(mio::isecir::InfectionState::InfectedNoSymptoms)] +
             this->m_transitions[i - init_start_index]
                                [Eigen::Index(mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms)]);
        ;
    }

    // flow from I to H and from I to R for -global_support_max, ..., 0
    for (int i = init_start_index; i <= t0_ide_index; i++) {
        this->m_transitions[i - init_start_index]
                           [Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)] =
            model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0] *
            (secihurd_ode[i - 1][Eigen::Index(mio::isecir::InfectionState::InfectedSymptoms)] -
             secihurd_ode[i][Eigen::Index(mio::isecir::InfectionState::InfectedSymptoms)] +
             this->m_transitions[i - init_start_index]
                                [Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]);
        this->m_transitions[i - init_start_index]
                           [Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToRecovered)] =
            (1 - model_ode.parameters.get<mio::osecir::SeverePerInfectedSymptoms>()[(mio::AgeGroup)0]) *
            (secihurd_ode[i - 1][Eigen::Index(mio::isecir::InfectionState::InfectedSymptoms)] -
             secihurd_ode[i][Eigen::Index(mio::isecir::InfectionState::InfectedSymptoms)] +
             this->m_transitions[i - init_start_index]
                                [Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)]);
    }

    // flow from H to U and from H to R for -global_support_max, ..., 0
    for (int i = init_start_index; i <= t0_ide_index; i++) {
        this->m_transitions[i - init_start_index]
                           [Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)] =
            model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0] *
            (secihurd_ode[i - 1][Eigen::Index(mio::isecir::InfectionState::InfectedSevere)] -
             secihurd_ode[i][Eigen::Index(mio::isecir::InfectionState::InfectedSevere)] +
             this->m_transitions[i - init_start_index]
                                [Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)]);
        this->m_transitions[i - init_start_index]
                           [Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToRecovered)] =
            (1 - model_ode.parameters.get<mio::osecir::CriticalPerSevere>()[(mio::AgeGroup)0]) *
            (secihurd_ode[i - 1][Eigen::Index(mio::isecir::InfectionState::InfectedSevere)] -
             secihurd_ode[i][Eigen::Index(mio::isecir::InfectionState::InfectedSevere)] +
             this->m_transitions[i - init_start_index]
                                [Eigen::Index(mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere)]);
    }

    // flow from U to D and from U to R for -global_support_max, ..., 0
    for (int i = init_start_index; i <= t0_ide_index; i++) {
        this->m_transitions[i - init_start_index]
                           [Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToDead)] =
            model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0] *
            (secihurd_ode[i - 1][Eigen::Index(mio::isecir::InfectionState::InfectedCritical)] -
             secihurd_ode[i][Eigen::Index(mio::isecir::InfectionState::InfectedCritical)] +
             this->m_transitions[i - init_start_index]
                                [Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)]);
        this->m_transitions[i - init_start_index]
                           [Eigen::Index(mio::isecir::InfectionTransition::InfectedCriticalToRecovered)] =
            (1 - model_ode.parameters.get<mio::osecir::DeathsPerCritical>()[(mio::AgeGroup)0]) *
            (secihurd_ode[i - 1][Eigen::Index(mio::isecir::InfectionState::InfectedCritical)] -
             secihurd_ode[i][Eigen::Index(mio::isecir::InfectionState::InfectedCritical)] +
             this->m_transitions[i - init_start_index]
                                [Eigen::Index(mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical)]);
    }
}

// define function that computes flows needed for initalization of IDE for a given result/compartments of the ODE model
// we assume that the ODE simulation starts at t0=0
void Model::compute_initial_flows_from_flows(mio::TimeSeries<ScalarType> secihurd_ode, ScalarType t0_ide, ScalarType dt)
{
    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // get (global) support_max to determine how many flows in the past we have to compute
    ScalarType global_support_max         = this->get_global_support_max(dt);
    Eigen::Index global_support_max_index = std::ceil(global_support_max / dt);

    // remove time point
    this->m_transitions.remove_last_time_point();

    ScalarType t0_ide_index = std::ceil(t0_ide / dt);
    unused(secihurd_ode);
    // flow from S to E for -6*global_support_max, ..., 0 (directly from compartments)
    // add time points to init_transitions here
    for (int i = t0_ide_index - 6 * global_support_max_index + 1; i <= t0_ide_index; i++) {
        this->m_transitions.add_time_point(i * dt, mio::TimeSeries<ScalarType>::Vector::Constant(num_transitions, 0));
        this->m_transitions.get_last_value()[Eigen::Index(mio::isecir::InfectionTransition::SusceptibleToExposed)] =
            secihurd_ode[i - 1][Eigen::Index(mio::isecir::InfectionState::Susceptible)] -
            secihurd_ode[i][Eigen::Index(mio::isecir::InfectionState::Susceptible)];
    }

    std::cout << "Computing flows for initialization. \n";

    // then use compute_flow function to compute following flows

    Eigen::Index start_shift = t0_ide_index - 6 * global_support_max_index;

    // flow from E to C for -5*global_support_max, ..., 0
    for (int i = t0_ide_index - 5 * global_support_max_index + 1; i <= t0_ide_index; i++) {
        this->compute_flow(1, 0, dt, true, i - start_shift);
    }

    // flow from C to I for -4*global_support_max, ..., 0
    for (int i = t0_ide_index - 4 * global_support_max_index + 1; i <= t0_ide_index; i++) {
        // C to I
        this->compute_flow(2, 1, dt, true, i - start_shift);
    }

    // flow from I to H for -3*global_support_max, ..., 0
    for (int i = t0_ide_index - 3 * global_support_max_index + 1; i <= t0_ide_index; i++) {
        // I to H
        this->compute_flow(4, 2, dt, true, i - start_shift);
    }

    // flow from H to U for -2*global_support_max, ..., 0
    for (int i = t0_ide_index - 2 * global_support_max_index + 1; i <= t0_ide_index; i++) {
        // H to U
        this->compute_flow(6, 4, dt, true, i - start_shift);
    }

    // flow from U to D and C, I, H, U to R for -1*global_support_max, ..., 0
    for (int i = t0_ide_index - 1 * global_support_max_index + 1; i <= t0_ide_index; i++) {
        // U to D
        this->compute_flow(8, 6, dt, true, i - start_shift);
        // C to R
        this->compute_flow(3, 1, dt, true, i - start_shift);
        // I to R
        this->compute_flow(5, 2, dt, true, i - start_shift);
        // H to R
        this->compute_flow(7, 4, dt, true, i - start_shift);
        // U to R
        this->compute_flow(9, 6, dt, true, i - start_shift);
    }
}

/***********************
* Solver for IDE model *
***********************/

void Model::initialize_solver(ScalarType dt)
{
    // TODO: think of where its best to add time point to m_populations,
    // especially wrt to check_constraints function
    // add first timepoint to m_populations at last time from m_transitions
    // fo now do it here to get the right time from m_transitions after initialization of flows
    // but this way we can't write in m_populations
    // another way would be to leave it in the constructor and overwrite it after initializing the flows
    m_populations.add_time_point<Eigen::VectorXd>(
        m_transitions.get_last_time(), TimeSeries<ScalarType>::Vector::Constant((int)InfectionState::Count, 0));

    // compute deaths at time t0
    m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Dead)] =
        m_deaths_before + m_transitions.get_last_value()[Eigen::Index(InfectionTransition::InfectedCriticalToDead)];

    // compute Susceptibles at time 0  and m_forceofinfection at time -m_dt as initial values for discretization scheme
    // use m_forceofinfection at -m_dt to be consistent with further calculations of S (see compute_susceptibles()),
    // where also the value of m_forceofinfection for the previous timestep is used
    update_forceofinfection(dt, true);

    if (m_forceofinfection > 0) {
        std::cout << "Computing compartments at time 0. \n";
        m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] =
            m_transitions.get_last_value()[Eigen::Index(InfectionTransition::SusceptibleToExposed)] /
            (dt * m_forceofinfection);

        //calculate other compartment sizes for t=0
        other_compartments_current_timestep(dt);

        //R; need an initial value for R, therefore do not calculate via compute_recovered()
        m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Recovered)] =
            m_N - m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Exposed)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::InfectedNoSymptoms)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSymptoms)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSevere)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::InfectedCritical)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Dead)];
    }
    else if (m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] > 1e-12) {
        //take initialized value for Susceptibles if value can't be calculated via the standard formula
        //calculate other compartment sizes for t=0
        other_compartments_current_timestep(dt);

        //R; need an initial value for R, therefore do not calculate via compute_recovered()
        m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Recovered)] =
            m_N - m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Exposed)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::InfectedNoSymptoms)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSymptoms)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSevere)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::InfectedCritical)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Dead)];
    }
    else if (m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Recovered)] > 1e-12) {
        //if value for Recovered is initialized and standard method is not applicable, calculate Susceptibles via other compartments
        //determining other compartment sizes is not dependent of Susceptibles(0), just of the transitions of the past.
        //calculate other compartment sizes for t=0
        other_compartments_current_timestep(dt);

        m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Susceptible)] =
            m_N - m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Exposed)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::InfectedNoSymptoms)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSymptoms)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::InfectedSevere)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::InfectedCritical)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Recovered)] -
            m_populations[Eigen::Index(0)][Eigen::Index(InfectionState::Dead)];
    }
    else {
        log_error("Error occured while initializing compartments: Force of infection is evaluated to 0 and neither "
                  "Susceptibles nor Recovered for time 0 were set. One of them should be larger 0.");
    }

    // compute m_forceofinfection at time 0 needed for further simulation
    update_forceofinfection(dt);
}

void Model::compute_susceptibles(ScalarType dt)
{
    Eigen::Index num_time_points = m_populations.get_num_time_points();
    // using number of susceptibles from previous time step and force of infection from previous time step:
    // compute current number of susceptibles and store susceptibles in m_populations
    m_populations.get_last_value()[Eigen::Index(InfectionState::Susceptible)] =
        m_populations[num_time_points - 2][Eigen::Index(InfectionState::Susceptible)] / (1 + dt * m_forceofinfection);
}

void Model::compute_flow(int idx_InfectionTransitions, Eigen::Index idx_IncomingFlow, ScalarType dt, bool initial_flow,
                         Eigen::Index current_time_index)
{
    ScalarType sum = 0;
    /* In order to satisfy TransitionDistribution(m_dt*i) = 0 for all i >= k, k is determined by the maximum support of the distribution.
    Since we are using a backwards difference scheme to compute the derivative, we have that the
    derivative of TransitionDistribution(m_dt*i) = 0 for all i >= k+1.

    Hence calc_time_index goes until std::ceil(support_max/m_dt) since for std::ceil(support_max/m_dt)+1 all terms are already zero. 
    This needs to be adjusted if we are changing the finite difference scheme */

    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(
        parameters.get<TransitionDistributions>()[idx_InfectionTransitions].get_support_max(dt) / dt);

    Eigen::Index num_time_points;

    if (initial_flow) {
        num_time_points = current_time_index;
    }

    else {
        num_time_points = m_transitions.get_num_time_points();
    }

    for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points - 1; i++) {

        // (num_time_points - 1 - i)* m_dt is the time, the individuals has already spent in this state.
        ScalarType state_age = (num_time_points - 1 - i) * dt;

        // backward difference scheme to approximate first derivative
        sum += (parameters.get<TransitionDistributions>()[idx_InfectionTransitions].eval(state_age) -
                parameters.get<TransitionDistributions>()[idx_InfectionTransitions].eval(state_age - dt)) /
               dt * m_transitions[i + 1][idx_IncomingFlow];
    }

    m_transitions.get_value(num_time_points - 1)[Eigen::Index(idx_InfectionTransitions)] =
        (-dt) * parameters.get<TransitionProbabilities>()[idx_InfectionTransitions] * sum;
}

void Model::flows_current_timestep(ScalarType dt)
{
    // calculate flow from S to E with force of infection from previous time step und susceptibles from current time step
    m_transitions.get_last_value()[Eigen::Index(InfectionTransition::SusceptibleToExposed)] =
        dt * m_forceofinfection * m_populations.get_last_value()[Eigen::Index(InfectionState::Susceptible)];
    // calculate all other flows with compute_flow
    // flow from E to C
    compute_flow((int)InfectionTransition::ExposedToInfectedNoSymptoms,
                 Eigen::Index(InfectionTransition::SusceptibleToExposed), dt);
    // flow from C to I
    compute_flow((int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms,
                 Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt);
    // flow from C to R
    compute_flow((int)InfectionTransition::InfectedNoSymptomsToRecovered,
                 Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms), dt);
    // flow from I to H
    compute_flow((int)InfectionTransition::InfectedSymptomsToInfectedSevere,
                 Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt);
    // flow from I to R
    compute_flow((int)InfectionTransition::InfectedSymptomsToRecovered,
                 Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms), dt);
    // flow from H to U
    compute_flow((int)InfectionTransition::InfectedSevereToInfectedCritical,
                 Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt);
    // flow from to H to R
    compute_flow((int)InfectionTransition::InfectedSevereToRecovered,
                 Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere), dt);
    // flow from U to D
    compute_flow((int)InfectionTransition::InfectedCriticalToDead,
                 Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt);
    // flow from U to R
    compute_flow((int)InfectionTransition::InfectedCriticalToRecovered,
                 Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical), dt);
}

void Model::compute_deaths()
{
    Eigen::Index num_time_points = m_populations.get_num_time_points();

    m_populations.get_last_value()[Eigen::Index(InfectionState::Dead)] =
        m_populations[num_time_points - 2][Eigen::Index(InfectionState::Dead)] +
        m_transitions.get_last_value()[Eigen::Index(InfectionTransition::InfectedCriticalToDead)];
}

void Model::update_forceofinfection(ScalarType dt, bool initialization)
{
    m_forceofinfection = 0;

    // determine the relevant calculation area = union of the supports of the relevant transition distributions
    ScalarType calc_time = std::max(
        {parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
             .get_support_max(dt),
         parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedNoSymptomsToRecovered]
             .get_support_max(dt),
         parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSymptomsToInfectedSevere]
             .get_support_max(dt),
         parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSymptomsToRecovered]
             .get_support_max(dt)});

    // corresponding index
    /* need calc_time_index timesteps in sum,
     subtract 1 because in the last summand all TransitionDistributions evaluate to 0 (by definition of support_max)*/
    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(calc_time / dt) - 1;

    Eigen::Index num_time_points;
    ScalarType current_time;
    ScalarType deaths;

    if (initialization) {
        // determine m_forceofinfection at time -m_dt which is the penultimate timepoint in m_transitions
        num_time_points = m_transitions.get_num_time_points() - 1;
        current_time    = -dt;
        deaths          = m_deaths_before;
    }
    else {
        // determine m_forceofinfection for current last time in m_transitions.
        num_time_points = m_transitions.get_num_time_points();
        current_time    = m_transitions.get_last_time();
        deaths          = m_populations.get_last_value()[Eigen::Index(InfectionState::Dead)];
    }

    for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points - 1; i++) {

        ScalarType state_age = (num_time_points - 1 - i) * dt;

        m_forceofinfection +=
            parameters.get<TransmissionProbabilityOnContact>().eval(state_age) *
            parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(current_time)(0, 0) *
            ((parameters
                      .get<TransitionProbabilities>()[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] *
                  parameters
                      .get<TransitionDistributions>()[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
                      .eval(state_age) +
              parameters.get<TransitionProbabilities>()[(int)InfectionTransition::InfectedNoSymptomsToRecovered] *
                  parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedNoSymptomsToRecovered]
                      .eval(state_age)) *
                 m_transitions[i + 1][Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms)] *
                 parameters.get<RelativeTransmissionNoSymptoms>().eval(state_age) +
             (parameters.get<TransitionProbabilities>()[(int)InfectionTransition::InfectedSymptomsToInfectedSevere] *
                  parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSymptomsToInfectedSevere]
                      .eval(state_age) +
              parameters.get<TransitionProbabilities>()[(int)InfectionTransition::InfectedSymptomsToRecovered] *
                  parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSymptomsToRecovered].eval(
                      state_age)) *
                 m_transitions[i + 1][Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] *
                 parameters.get<RiskOfInfectionFromSymptomatic>().eval(state_age));
    }
    m_forceofinfection = 1 / (m_N - deaths) * m_forceofinfection;
}

void Model::compute_compartment(Eigen::Index idx_InfectionState, Eigen::Index idx_IncomingFlow,
                                int idx_TransitionDistribution1, int idx_TransitionDistribution2, ScalarType dt)
{
    ScalarType sum = 0;

    // determine relevant calculation area and corresponding index
    ScalarType calc_time =
        std::max(parameters.get<TransitionDistributions>()[idx_TransitionDistribution1].get_support_max(dt),
                 parameters.get<TransitionDistributions>()[idx_TransitionDistribution2].get_support_max(dt));

    Eigen::Index calc_time_index = (Eigen::Index)std::ceil(calc_time / dt) - 1;

    Eigen::Index num_time_points = m_transitions.get_num_time_points();

    for (Eigen::Index i = num_time_points - 1 - calc_time_index; i < num_time_points - 1; i++) {

        ScalarType state_age = (num_time_points - 1 - i) * dt;

        sum += (parameters.get<TransitionProbabilities>()[idx_TransitionDistribution1] *
                    parameters.get<TransitionDistributions>()[idx_TransitionDistribution1].eval(state_age) +
                (1 - parameters.get<TransitionProbabilities>()[idx_TransitionDistribution1]) *
                    parameters.get<TransitionDistributions>()[idx_TransitionDistribution2].eval(state_age)) *
               m_transitions[i + 1][idx_IncomingFlow];
    }

    m_populations.get_last_value()[idx_InfectionState] = sum;
}

void Model::other_compartments_current_timestep(ScalarType dt)
{
    // E
    compute_compartment(Eigen::Index(InfectionState::Exposed), Eigen::Index(InfectionTransition::SusceptibleToExposed),
                        (int)InfectionTransition::ExposedToInfectedNoSymptoms, 0,
                        dt); // this is a dummy index as there is no transition from E to R in our model,
    // write any transition here as probability from E to R is 0
    // C
    compute_compartment(Eigen::Index(InfectionState::InfectedNoSymptoms),
                        Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                        (int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms,
                        (int)InfectionTransition::InfectedNoSymptomsToRecovered, dt);
    // I
    compute_compartment(Eigen::Index(InfectionState::InfectedSymptoms),
                        Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
                        (int)InfectionTransition::InfectedSymptomsToInfectedSevere,
                        (int)InfectionTransition::InfectedSymptomsToRecovered, dt);
    // H
    compute_compartment(Eigen::Index(InfectionState::InfectedSevere),
                        Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere),
                        (int)InfectionTransition::InfectedSevereToInfectedCritical,
                        (int)InfectionTransition::InfectedSevereToRecovered, dt);
    // U
    compute_compartment(Eigen::Index(InfectionState::InfectedCritical),
                        Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical),
                        (int)InfectionTransition::InfectedCriticalToDead,
                        (int)InfectionTransition::InfectedCriticalToRecovered, dt);
}

void Model::compute_recovered()
{
    Eigen::Index num_time_points = m_populations.get_num_time_points();

    m_populations.get_last_value()[Eigen::Index(InfectionState::Recovered)] =
        m_populations[num_time_points - 2][Eigen::Index(InfectionState::Recovered)] +
        m_transitions.get_last_value()[Eigen::Index(InfectionTransition::InfectedNoSymptomsToRecovered)] +
        m_transitions.get_last_value()[Eigen::Index(InfectionTransition::InfectedSymptomsToRecovered)] +
        m_transitions.get_last_value()[Eigen::Index(InfectionTransition::InfectedSevereToRecovered)] +
        m_transitions.get_last_value()[Eigen::Index(InfectionTransition::InfectedCriticalToRecovered)];
}

/*******************
* Helper functions *
*******************/

ScalarType Model::get_global_support_max(ScalarType dt) const
{
    ScalarType global_support_max = std::max(
        {parameters.get<TransitionDistributions>()[(int)InfectionTransition::ExposedToInfectedNoSymptoms]
             .get_support_max(dt),
         parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedNoSymptomsToInfectedSymptoms]
             .get_support_max(dt),
         parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedNoSymptomsToRecovered]
             .get_support_max(dt),
         parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSymptomsToInfectedSevere]
             .get_support_max(dt),
         parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSymptomsToRecovered]
             .get_support_max(dt),
         parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSevereToInfectedCritical]
             .get_support_max(dt),
         parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedSevereToRecovered].get_support_max(
             dt),
         parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedCriticalToDead].get_support_max(
             dt),
         parameters.get<TransitionDistributions>()[(int)InfectionTransition::InfectedCriticalToRecovered]
             .get_support_max(dt)});

    return global_support_max;
}

} // namespace isecir
} // namespace mio
