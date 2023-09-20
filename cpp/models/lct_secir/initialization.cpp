/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke
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
#include "lct_secir/initialization.h"
#include "lct_secir/infection_state.h"
#include "lct_secir/parameters.h"
#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include "memilio/math/floating_point.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/logging.h"
#include "memilio/epidemiology/state_age_function.h"

namespace mio
{
namespace lsecir
{

void Initializer::check_constraints() const
{
    if (!((int)InfectionTransition::Count == (int)m_flows.get_num_elements())) {
        log_error("Initial condition size does not match Subcompartments.");
    }

    parameters.check_constraints();

    for (int i = 1; i < m_flows.get_num_time_points(); i++) {
        if (!(floating_point_equal(m_dt, m_flows.get_time(i) - m_flows.get_time(i - 1), 1.0, 1e-14))) {
            log_error("Time points in flows have to be equidistant.");
        }
    }

    if (!(m_dt < 1)) {
        log_warning("Step size was set very large. The result could be distorted.");
    }
}

Eigen::VectorXd Initializer::compute_compartment(InfectionStateBase base, Eigen::Index idx_incoming_flow,
                                                 ScalarType transition_rate) const
{
    int num_infectionstates = infectionStates.get_number(base);
    Eigen::VectorXd subcompartments(num_infectionstates);
    // Initialize relevant density for the compartment base.
    // For the first subcompartment a shape parameter of one is needed.
    ErlangDensity erlang(1, 1. / (num_infectionstates * transition_rate));

    // Initialize other relevant parameters.
    ScalarType calc_time{0};
    Eigen::Index calc_time_index{0};
    ScalarType state_age{0};
    ScalarType sum{0};
    Eigen::Index num_time_points = m_flows.get_num_time_points();

    // Calculate number of people in each subcomaprtment.
    for (int j = 0; j < num_infectionstates; j++) {
        // For subcompartment number j+1, shape parameter j+1 is needed.
        erlang.set_parameter(j + 1);

        // Determine relevant calculation area and corresponding index.
        calc_time       = erlang.get_support_max(m_dt, 1e-6);
        calc_time_index = (Eigen::Index)std::ceil(calc_time / m_dt) - 1;

        if (num_time_points < calc_time_index) {
            log_error("Initialization failed. Not enough time points for the transitions are given.  {} are needed but "
                      "just {} are given.",
                      calc_time_index, num_time_points);
        }

        // Approximate integral with non-standard scheme.
        for (Eigen::Index i = num_time_points - calc_time_index; i < num_time_points; i++) {
            state_age = (num_time_points - i) * m_dt;
            sum += erlang.eval(state_age) * m_flows[i][idx_incoming_flow];
        }
        subcompartments[j] = 1 / (num_infectionstates * transition_rate) * sum;
        sum                = 0;
    }
    return subcompartments;
}

Eigen::VectorXd Initializer::compute_initializationvector(ScalarType total_population, ScalarType deaths,
                                                          ScalarType total_confirmed_cases) const
{
    check_constraints();

    int infectionStates_count = infectionStates.get_count();
    Eigen::VectorXd init(infectionStates_count);

    //E
    init.segment(infectionStates.get_firstindex(InfectionStateBase::Exposed),
                 infectionStates.get_number(InfectionStateBase::Exposed)) =
        compute_compartment(InfectionStateBase::Exposed, Eigen::Index(InfectionTransition::SusceptibleToExposed),
                            1 / parameters.get<TimeExposed>());
    //C
    init.segment(infectionStates.get_firstindex(InfectionStateBase::InfectedNoSymptoms),
                 infectionStates.get_number(InfectionStateBase::InfectedNoSymptoms)) =
        compute_compartment(InfectionStateBase::InfectedNoSymptoms,
                            Eigen::Index(InfectionTransition::ExposedToInfectedNoSymptoms),
                            1 / parameters.get<TimeInfectedNoSymptoms>());
    //I
    init.segment(infectionStates.get_firstindex(InfectionStateBase::InfectedSymptoms),
                 infectionStates.get_number(InfectionStateBase::InfectedSymptoms)) =
        compute_compartment(InfectionStateBase::InfectedSymptoms,
                            Eigen::Index(InfectionTransition::InfectedNoSymptomsToInfectedSymptoms),
                            1 / parameters.get<TimeInfectedSymptoms>());
    //H
    init.segment(infectionStates.get_firstindex(InfectionStateBase::InfectedSevere),
                 infectionStates.get_number(InfectionStateBase::InfectedSevere)) =
        compute_compartment(InfectionStateBase::InfectedSevere,
                            Eigen::Index(InfectionTransition::InfectedSymptomsToInfectedSevere),
                            1 / parameters.get<TimeInfectedSevere>());
    //U
    init.segment(infectionStates.get_firstindex(InfectionStateBase::InfectedCritical),
                 infectionStates.get_number(InfectionStateBase::InfectedCritical)) =
        compute_compartment(InfectionStateBase::InfectedCritical,
                            Eigen::Index(InfectionTransition::InfectedSevereToInfectedCritical),
                            1 / parameters.get<TimeInfectedCritical>());
    //R
    // Number of recovered is equal to the cumulative number of confirmed cases minus the number of people who are infected at the moment.
    init[infectionStates_count - 2] = total_confirmed_cases -
                                      init.segment(infectionStates.get_firstindex(InfectionStateBase::InfectedSymptoms),
                                                   infectionStates.get_number(InfectionStateBase::InfectedSymptoms) +
                                                       infectionStates.get_number(InfectionStateBase::InfectedSevere) +
                                                       infectionStates.get_number(InfectionStateBase::InfectedCritical))
                                          .sum() -
                                      deaths;

    //S
    init[0] = total_population - init.segment(1, infectionStates_count - 2).sum() - deaths;

    //D
    init[infectionStates_count - 1] = deaths;

    return init;
}

} // namespace lsecir
} // namespace mio