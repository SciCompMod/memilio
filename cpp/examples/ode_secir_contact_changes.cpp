/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Martin J. Kuehn
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
#include "ode_secir/model.h"
#include "memilio/compartments/flow_simulation.h"
#include "memilio/utils/logging.h"
#include "memilio/math/euler.h"

/*
 * This example demonstrates how to realize contact behavior changes with any of our ordinary differential
 * equation-based models. This example can thus easily be adapted for other models like osecirvvs, oseir, osir etc.
 * We print out the flows (i.e., new transmissions, infections, hospitalizations etc. per time point.) As we use an
 * fixed-time step Explicit Euler, we can compare them.
*/
int main()
{
    mio::set_log_level(mio::LogLevel::warn);

    double t0   = 0;
    double dt   = 0.1;
    double tmax = 1;

    double cont_freq = 10;

    double nb_total_t0 = 1000, nb_inf_t0 = 10;

    // default model run to be compared against
    mio::osecir::Model<double> model_a(1);
    const auto indx_flow_SE =
        model_a.get_flat_flow_index<mio::osecir::InfectionState::Susceptible, mio::osecir::InfectionState::Exposed>(
            {mio::AgeGroup(0)});

    model_a.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] = nb_inf_t0;
    model_a.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                  nb_total_t0);
    mio::ContactMatrixGroup<double>& contact_matrix_a = model_a.parameters.get<mio::osecir::ContactPatterns<double>>();
    contact_matrix_a[0] = mio::ContactMatrix<double>(Eigen::MatrixX<ScalarType>::Constant(1, 1, cont_freq));
    // set probability of transmission and risk of infection to 1.
    model_a.parameters.get<mio::osecir::TransmissionProbabilityOnContact<ScalarType>>() = 1.0;
    model_a.parameters.get<mio::osecir::RiskOfInfectionFromSymptomatic<ScalarType>>()   = 1.0;

    mio::EulerIntegratorCore<ScalarType> integrator;
    auto result_a = mio::simulate_flows<ScalarType>(t0, tmax, dt, model_a, integrator.clone());
    result_a[1].print_table({"S->E", "E->I_NS", "I_NS->I_Sy", "I_NS->R", "I_NSC->I_SyC", "I_NSC->R", "I_Sy->I_Sev",
                             "I_Sy->R", "I_SyC->I_Sev", "I_SyC->R", "I_Sev->I_Crit", "I_Sev->R", "I_Sev->D",
                             "I_Crit->D", "I_Crit->R"},
                            4, 4);
    std::cout << "With default contacts, the number of new transmissions (flow from S->E) in first time step is: "
              << result_a[1].get_value(1)[indx_flow_SE] << ".\n";

    // The contacts are halfed: reduced transmission through damping with value 0.5
    mio::osecir::Model<double> model_b{model_a};
    model_b.populations.set_total(nb_total_t0);
    model_b.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] = nb_inf_t0;
    model_b.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                  nb_total_t0);
    mio::ContactMatrixGroup<ScalarType>& contact_matrix_b =
        model_b.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix_b[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, cont_freq));
    contact_matrix_b[0].add_damping(0.5, mio::SimulationTime<ScalarType>(0.)); // contact reduction happens here!

    auto result_b = mio::simulate_flows<ScalarType>(t0, tmax, dt, model_b, integrator.clone());
    result_b[1].print_table({"S->E", "E->I_NS", "I_NS->I_Sy", "I_NS->R", "I_NSC->I_SyC", "I_NSC->R", "I_Sy->I_Sev",
                             "I_Sy->R", "I_SyC->I_Sev", "I_SyC->R", "I_Sev->I_Crit", "I_Sev->R", "I_Sev->D",
                             "I_Crit->D", "I_Crit->R"},
                            4, 4);
    std::cout << "With contacts reduced to a half of the original example, the number of new transmissions (flow from "
                 "S->E) in first time step is: "
              << result_b[1].get_value(1)[indx_flow_SE] << ".\n";

    // No contacts at all: no transmission through damping with value 1.
    mio::osecir::Model<double> model_c{model_a};
    model_c.populations.set_total(nb_total_t0);
    model_c.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] = nb_inf_t0;
    model_c.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                  nb_total_t0);
    mio::ContactMatrixGroup<ScalarType>& contact_matrix_c =
        model_c.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix_c[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, cont_freq));
    contact_matrix_c[0].add_damping(1., mio::SimulationTime<ScalarType>(0.)); // contact reduction happens here!

    auto result_c = mio::simulate_flows<ScalarType>(t0, tmax, dt, model_c, integrator.clone());
    result_c[1].print_table({"S->E", "E->I_NS", "I_NS->I_Sy", "I_NS->R", "I_NSC->I_SyC", "I_NSC->R", "I_Sy->I_Sev",
                             "I_Sy->R", "I_SyC->I_Sev", "I_SyC->R", "I_Sev->I_Crit", "I_Sev->R", "I_Sev->D",
                             "I_Crit->D", "I_Crit->R"},
                            4, 4);
    std::cout
        << "With contacts reduced to zero, the number of new transmissions (flow from S->E) in first time step is: "
        << result_c[1].get_value(1)[indx_flow_SE] << ".\n";

    // The contacts are doubled: increased transmission through damping with value -1.
    mio::osecir::Model model_d{model_a};
    model_d.populations.set_total(nb_total_t0);
    model_d.populations[{mio::AgeGroup(0), mio::osecir::InfectionState::InfectedSymptoms}] = nb_inf_t0;
    model_d.populations.set_difference_from_total({mio::AgeGroup(0), mio::osecir::InfectionState::Susceptible},
                                                  nb_total_t0);
    mio::ContactMatrixGroup<ScalarType>& contact_matrix_d =
        model_d.parameters.get<mio::osecir::ContactPatterns<ScalarType>>();
    contact_matrix_d[0] = mio::ContactMatrix<ScalarType>(Eigen::MatrixX<ScalarType>::Constant(1, 1, cont_freq));
    contact_matrix_d[0].add_damping(-1., mio::SimulationTime<ScalarType>(0.)); // contact increase happens here!

    auto result_d = mio::simulate_flows<ScalarType>(t0, tmax, dt, model_d, integrator.clone());
    result_d[1].print_table({"S->E", "E->I_NS", "I_NS->I_Sy", "I_NS->R", "I_NSC->I_SyC", "I_NSC->R", "I_Sy->I_Sev",
                             "I_Sy->R", "I_SyC->I_Sev", "I_SyC->R", "I_Sev->I_Crit", "I_Sev->R", "I_Sev->D",
                             "I_Crit->D", "I_Crit->R"},
                            4, 4);
    std::cout << "With contacts doubled, the number of new transmissions (flow from S->E) in first time step is: "
              << result_d[1].get_value(1)[indx_flow_SE] << "\n";
}
