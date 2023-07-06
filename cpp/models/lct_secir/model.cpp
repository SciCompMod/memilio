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
#include "lct_secir/model.h"
#include "Eigen/src/Core/Matrix.h"
#include "infection_state.h"
#include "memilio/config.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include <string>

namespace mio
{
namespace lsecir
{

Model::Model(Eigen::VectorXd init, const InfectionState InfectionState_init, const ParameterSet& Parameterset_init)
    : parameters{Parameterset_init}
    , infectionStates{InfectionState_init}
    , m_initial_values{std::move(init)}
{
    m_N0 = m_initial_values.sum();
}

void Model::check_constraints() const
{
    if (!(infectionStates.get_count() == m_initial_values.size())) {
        log_error("Initial condition size does not match Subcompartments.");
    }
    parameters.check_constraints();
}

void Model::eval_right_hand_side(Eigen::Ref<const Eigen::VectorXd> y, double t, Eigen::Ref<Eigen::VectorXd> dydt) const
{
    dydt.setZero();

    ScalarType C     = 0;
    ScalarType I     = 0;
    ScalarType dummy = 0;

    //calculate sum of all subcompartments for InfectedNoSymptoms
    C = y.segment(infectionStates.get_firstindex(InfectionStateBase::InfectedNoSymptoms),
                  infectionStates.get_number(InfectionStateBase::InfectedNoSymptoms))
            .sum();
    //calculate sum of all subcompartments for InfectedSymptoms
    I = y.segment(infectionStates.get_firstindex(InfectionStateBase::InfectedSymptoms),
                  infectionStates.get_number(InfectionStateBase::InfectedSymptoms))
            .sum();

    //S'
    dydt[0] =
        -y[0] / (m_N0 - y[infectionStates.get_firstindex(InfectionStateBase::Dead)]) *
        parameters.get<TransmissionProbabilityOnContact>() *
        parameters.get<ContactPatterns>().get_cont_freq_mat().get_matrix_at(t)(0, 0) *
        (parameters.get<RelativeTransmissionNoSymptoms>() * C + parameters.get<RiskOfInfectionFromSymptomatic>() * I);

    //E'
    dydt[1] = -dydt[0];
    for (int i = 0; i < infectionStates.get_number(InfectionStateBase::Exposed); i++) {
        // dummy stores the value of the flow from dydt[1 + i] to dydt[2 + i].
        // 1+i is always the index of a (sub-)compartment of E and 2+i can also be the index of the first (sub-)compartment of C
        dummy =
            infectionStates.get_number(InfectionStateBase::Exposed) * (1 / parameters.get<TimeExposed>()) * y[1 + i];
        // subtract flow from dydt[1 + i] and add to dydt[2 + i]
        dydt[1 + i] = dydt[1 + i] - dummy;
        dydt[2 + i] = dummy;
    }

    //C'
    for (int i = 0; i < infectionStates.get_number(InfectionStateBase::InfectedNoSymptoms); i++) {
        dummy = infectionStates.get_number(InfectionStateBase::InfectedNoSymptoms) *
                (1 / parameters.get<TimeInfectedNoSymptoms>()) *
                y[infectionStates.get_firstindex(InfectionStateBase::InfectedNoSymptoms) + i];
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedNoSymptoms) + i] =
            dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedNoSymptoms) + i] - dummy;
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedNoSymptoms) + i + 1] = dummy;
    }

    //I'
    // flow from last (sub-) compartment of C must be split between I_1 and R
    dydt[infectionStates.get_firstindex(InfectionStateBase::Recovered)] =
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedSymptoms)] *
        parameters.get<RecoveredPerInfectedNoSymptoms>();
    dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedSymptoms)] =
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedSymptoms)] *
        (1 - parameters.get<RecoveredPerInfectedNoSymptoms>());

    for (int i = 0; i < infectionStates.get_number(InfectionStateBase::InfectedSymptoms); i++) {
        dummy = infectionStates.get_number(InfectionStateBase::InfectedSymptoms) *
                (1 / parameters.get<TimeInfectedSymptoms>()) *
                y[infectionStates.get_firstindex(InfectionStateBase::InfectedSymptoms) + i];
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedSymptoms) + i] =
            dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedSymptoms) + i] - dummy;
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedSymptoms) + i + 1] = dummy;
    }

    // H'
    dydt[infectionStates.get_firstindex(InfectionStateBase::Recovered)] =
        dydt[infectionStates.get_firstindex(InfectionStateBase::Recovered)] +
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedSevere)] *
            (1 - parameters.get<SeverePerInfectedSymptoms>());
    dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedSevere)] =
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedSevere)] *
        parameters.get<SeverePerInfectedSymptoms>();
    for (int i = 0; i < infectionStates.get_number(InfectionStateBase::InfectedSevere); i++) {
        dummy = infectionStates.get_number(InfectionStateBase::InfectedSevere) *
                (1 / parameters.get<TimeInfectedSevere>()) *
                y[infectionStates.get_firstindex(InfectionStateBase::InfectedSevere) + i];
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedSevere) + i] =
            dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedSevere) + i] - dummy;
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedSevere) + i + 1] = dummy;
    }

    // U'
    dydt[infectionStates.get_firstindex(InfectionStateBase::Recovered)] =
        dydt[infectionStates.get_firstindex(InfectionStateBase::Recovered)] +
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedCritical)] *
            (1 - parameters.get<CriticalPerSevere>());
    dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedCritical)] =
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedCritical)] *
        parameters.get<CriticalPerSevere>();
    for (int i = 0; i < infectionStates.get_number(InfectionStateBase::InfectedCritical) - 1; i++) {
        dummy = infectionStates.get_number(InfectionStateBase::InfectedCritical) *
                (1 / parameters.get<TimeInfectedCritical>()) *
                y[infectionStates.get_firstindex(InfectionStateBase::InfectedCritical) + i];
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedCritical) + i] =
            dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedCritical) + i] - dummy;
        dydt[infectionStates.get_firstindex(InfectionStateBase::InfectedCritical) + i + 1] = dummy;
    }
    // last flow from U has to be divided between R and D.
    // Must be calculated separately in order not to overwrite the already calculated values ​​for R.
    dummy = infectionStates.get_number(InfectionStateBase::InfectedCritical) *
            (1 / parameters.get<TimeInfectedCritical>()) *
            y[infectionStates.get_firstindex(InfectionStateBase::Recovered) - 1];
    dydt[infectionStates.get_firstindex(InfectionStateBase::Recovered) - 1] =
        dydt[infectionStates.get_firstindex(InfectionStateBase::Recovered) - 1] - dummy;
    dydt[infectionStates.get_firstindex(InfectionStateBase::Recovered)] =
        dydt[infectionStates.get_firstindex(InfectionStateBase::Recovered)] +
        (1 - parameters.get<DeathsPerCritical>()) * dummy;
    dydt[infectionStates.get_firstindex(InfectionStateBase::Dead)] = parameters.get<DeathsPerCritical>() * dummy;
}

TimeSeries<ScalarType> Model::calculate_populations(const TimeSeries<ScalarType>& result) const
{
    if (!(infectionStates.get_count() == result.get_num_elements())) {
        log_error("Result does not match InfectionStates of the Model.");
    }
    TimeSeries<ScalarType> populations((int)InfectionStateBase::Count);
    Eigen::VectorXd dummy((int)InfectionStateBase::Count);
    for (Eigen::Index i = 0; i < result.get_num_time_points(); ++i) {
        for (int j = 0; j < dummy.size(); ++j) {
            // use segment of vector of the rsult with subcompartments of InfectionStateBase with index j and sum up values of subcompartments.
            dummy[j] = result[i]
                           .segment(Eigen::Index(infectionStates.get_firstindex(j)),
                                    Eigen::Index(infectionStates.get_number(j)))
                           .sum();
        }
        populations.add_time_point(result.get_time(i), dummy);
    }

    return populations;
}

std::string Model::get_heading_CompartmentsBase() const
{
    return "S | E | C | I | H | U | R | D";
}

std::string Model::get_heading_Subcompartments() const
{

    std::string heading = "S";
    std::string filler  = " | ";
    // if E has subcompartments in the model, append E1 to En to the heading
    if (infectionStates.get_number(InfectionStateBase::Exposed) > 1) {
        for (int i = 0; i < infectionStates.get_number(InfectionStateBase::Exposed); i++) {
            heading = heading + filler + "E" + std::to_string(i + 1);
        }
    }
    else { // if E does not have subcompartments in the model, just append E
        heading = heading + filler + "E";
    }
    if (infectionStates.get_number(InfectionStateBase::InfectedNoSymptoms) > 1) {
        for (int i = 0; i < infectionStates.get_number(InfectionStateBase::InfectedNoSymptoms); i++) {
            heading = heading + filler + "C" + std::to_string(i + 1);
        }
    }
    else {
        heading = heading + filler + "C";
    }
    if (infectionStates.get_number(InfectionStateBase::InfectedSymptoms) > 1) {
        for (int i = 0; i < infectionStates.get_number(InfectionStateBase::InfectedSymptoms); i++) {
            heading = heading + filler + "I" + std::to_string(i + 1);
        }
    }
    else {
        heading = heading + filler + "I";
    }
    if (infectionStates.get_number(InfectionStateBase::InfectedSevere) > 1) {
        for (int i = 0; i < infectionStates.get_number(InfectionStateBase::InfectedSevere); i++) {
            heading = heading + filler + "H" + std::to_string(i + 1);
        }
    }
    else {
        heading = heading + filler + "H";
    }
    if (infectionStates.get_number(InfectionStateBase::InfectedCritical) > 1) {
        for (int i = 0; i < infectionStates.get_number(InfectionStateBase::InfectedCritical); i++) {
            heading = heading + filler + "U" + std::to_string(i + 1);
        }
    }
    else {
        heading = heading + filler + "U";
    }
    heading = heading + filler + "R";
    heading = heading + filler + "D";
    return heading;
}

} // namespace lsecir
} // namespace mio
