/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Anna Wendler, Lena Ploetzke
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

#include "ide_secir/model.h"
#include "ide_secir/infection_state.h"
#include "ide_secir/parameters.h"
#include "ide_secir/simulation.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/time_series.h"
#include "memilio/epidemiology/uncertain_matrix.h"

int main()
{
    using Vec = mio::TimeSeries<ScalarType>::Vector;

    ScalarType tmax   = 10;
    ScalarType N      = 10000;
    ScalarType deaths = 13.10462213;
    ScalarType dt     = 1;

    int num_transitions = (int)mio::isecir::InfectionTransition::Count;

    // Create TimeSeries with num_transitions elements where transitions needed for simulation will be stored.
    mio::TimeSeries<ScalarType> init(num_transitions);

    //Add initial time point to TimeSeries.
    Vec vec_init(num_transitions);
    vec_init[(int)mio::isecir::InfectionTransition::SusceptibleToExposed]                 = 25.0;
    vec_init[(int)mio::isecir::InfectionTransition::ExposedToInfectedNoSymptoms]          = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms] = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToRecovered]        = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToInfectedSevere]     = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSymptomsToRecovered]          = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToInfectedCritical]     = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedSevereToRecovered]            = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToDead]               = 0.0;
    vec_init[(int)mio::isecir::InfectionTransition::InfectedCriticalToRecovered]          = 0.0;
    vec_init                                                                              = vec_init * dt;

    // Add initial time point to time series.
    init.add_time_point(-10, vec_init);
    // Add further time points until time 0.
    while (init.get_last_time() < 0) {
        init.add_time_point(init.get_last_time() + dt, vec_init);
    }

    // Initialize model.
    mio::isecir::Model model(std::move(init), N, deaths, 20);

    // Uncomment one of these lines to use a different method to initialize the model using the TimeSeries init.
    // model.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Susceptible] = 1000;
    //model.m_populations.get_last_value()[(Eigen::Index)mio::isecir::InfectionState::Recovered] = 1;

    // Set working parameters.
    mio::SmootherCosine expdecay2(4.0);
    mio::StateAgeFunctionWrapper delaydistribution(expdecay2);
    std::vector<mio::StateAgeFunctionWrapper> vec_delaydistrib(num_transitions, delaydistribution);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::SusceptibleToExposed].set_parameter(3.0);
    vec_delaydistrib[(int)mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms].set_parameter(4.0);
    model.parameters.set<mio::isecir::TransitionDistributions>(vec_delaydistrib);

    std::vector<ScalarType> vec_prob((int)mio::isecir::InfectionTransition::Count, 1.0);
    vec_prob[Eigen::Index(mio::isecir::InfectionTransition::InfectedNoSymptomsToInfectedSymptoms)] = 0.0;
    model.parameters.set<mio::isecir::TransitionProbabilities>(vec_prob);

    mio::ContactMatrixGroup contact_matrix               = mio::ContactMatrixGroup(1, 1);
    contact_matrix[0]                                    = mio::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, 1.));
    model.parameters.get<mio::isecir::ContactPatterns>() = mio::UncertainContactMatrix(contact_matrix);

    mio::ConstantFunction expdecay(0.0);
    mio::StateAgeFunctionWrapper prob(expdecay);
    model.parameters.set<mio::isecir::TransmissionProbabilityOnContact>(prob);
    model.parameters.set<mio::isecir::RelativeTransmissionNoSymptoms>(prob);
    model.parameters.set<mio::isecir::RiskOfInfectionFromSymptomatic>(prob);

    model.check_constraints(dt);

    // Carry out simulation.
    mio::isecir::Simulation sim(model, 0, dt);
    sim.advance(tmax);

    sim.get_result().print_table({"S", "E", "C", "I", "H", "U", "R", "D "}, 16, 8);
    sim.get_transitions().print_table({"S->E", "E->C", "C->I", "C->R", "I->H", "I->R", "H->U", "H->R", "U->D", "U->R"},
                                      16, 8);

    // mio::TimeSeries<ScalarType> result = sim.get_result();

    // std::cout << "Compartments at first time step of IDE:\n";
    // std::cout << "# time  |  S  |  E  |  C  |  I  |  H  |  U  |  R  |  D  |" << std::endl;
    // for (Eigen::Index j = 0; j < result.get_num_elements(); ++j) {
    //     std::cout << "  |  " << std::fixed << std::setprecision(8) << result[0][j];
    // }
    // std::cout << "\n";
    // std::cout << "Compartments at last time step of IDE:\n";
    // std::cout << "# time  |  S  |  E  |  C  |  I  |  H  |  U  |  R  |  D  |" << std::endl;
    // for (Eigen::Index j = 0; j < result.get_num_elements(); ++j) {
    //     std::cout << "  |  " << std::fixed << std::setprecision(8) << result.get_last_value()[j];
    // }
    // std::cout << "\n";
}
