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
#include "lct_secir/infection_state.h"
#include "lct_secir/simulation.h"
#include "memilio/config.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/math/eigen.h"
#include "memilio/io/result_io.h"

int main()
{
    // Set vector that specifies the number of subcompartments
    std::vector<int> SubcompartmentNumbers((int)mio::lsecir::InfectionStateBase::Count, 1);
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::Exposed]            = 2;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedNoSymptoms] = 3;
    SubcompartmentNumbers[(int)mio::lsecir::InfectionStateBase::InfectedCritical]   = 5;
    mio::lsecir::InfectionState InfState(SubcompartmentNumbers);

    ScalarType tmax = 20;
    // define initial population distribution in infection states, one entry per Subcompartment
    Eigen::VectorXd init(InfState.get_count());
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Susceptible)]            = 750;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Exposed)]                = 30;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Exposed) + 1]            = 20;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms)]     = 20;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms) + 1] = 10;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedNoSymptoms) + 2] = 10;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedSymptoms)]       = 50;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedSevere)]         = 50;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical)]       = 10;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 1]   = 10;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 2]   = 5;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 3]   = 3;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::InfectedCritical) + 4]   = 2;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Recovered)]              = 20;
    init[InfState.get_firstindex(mio::lsecir::InfectionStateBase::Dead)]                   = 10;

    // initialize model
    mio::lsecir::Model model(std::move(init), InfState);

    // Set Parameters of the model
    model.parameters.get<mio::lsecir::TimeExposed>()            = 2 * 4.2 - 5.2;
    model.parameters.get<mio::lsecir::TimeInfectedNoSymptoms>() = 2 * (5.2 - 4.2);
    model.parameters.get<mio::lsecir::TimeInfectedSymptoms>()   = 5.8;
    model.parameters.get<mio::lsecir::TimeInfectedSevere>()     = 9.5;
    // also possible to change values with setter
    model.parameters.set<mio::lsecir::TimeInfectedCritical>(7.1);

    model.parameters.get<mio::lsecir::TransmissionProbabilityOnContact>() = 0.05;

    mio::ContactMatrixGroup& contact_matrix = model.parameters.get<mio::lsecir::ContactPatterns>();
    contact_matrix[0]                       = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 10));
    contact_matrix[0].add_damping(0.7, mio::SimulationTime(5.));

    model.parameters.get<mio::lsecir::RelativeTransmissionNoSymptoms>() = 0.7;
    model.parameters.get<mio::lsecir::RiskOfInfectionFromSymptomatic>() = 0.25;
    model.parameters.get<mio::lsecir::RecoveredPerInfectedNoSymptoms>() = 0.09;
    model.parameters.get<mio::lsecir::SeverePerInfectedSymptoms>()      = 0.2;
    model.parameters.get<mio::lsecir::CriticalPerSevere>()              = 0.25;
    model.parameters.set<mio::lsecir::DeathsPerCritical>(0.3);

    // perform simulation
    mio::TimeSeries<ScalarType> result = mio::lsecir::simulate(0, tmax, 0.5, model);
    // print result with Subcompartments
    // mio::lsecir::print_TimeSeries(result, model.get_heading_Subcompartments());
    // calculate the distribution in Infectionstates without subcompartments of the result
    mio::TimeSeries<ScalarType> populations = model.calculate_populations(result);
    // print it
    mio::lsecir::print_TimeSeries(populations, model.get_heading_CompartmentsBase());

    bool save_result = false;
    if (save_result) {
        auto save_result_status_subcompartments = mio::save_result({result}, {0}, 1, "result_lct_subcompartments.h5");
        auto save_result_status                 = mio::save_result({populations}, {0}, 1, "result_lct.h5");
    }
}