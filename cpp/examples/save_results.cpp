/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Wadim Koslow
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
#include "secir/secir.h"
#include "secir/secir_result_io.h"

#include <iostream>

int main(int argc, char** argv)
{

    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 0.5, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    int nb_groups = 1;
    double fact   = 1.0 / (double)nb_groups;

    mio::SecirParams params(nb_groups);

    params.set_start_day(0);
    params.set_seasonality(0);

    for (size_t i = 0; i < nb_groups; i++) {
        params.times[i].set_incubation(tinc);
        params.times[i].set_infectious_mild(tinfmild);
        params.times[i].set_serialinterval(tserint);
        params.times[i].set_hospitalized_to_home(thosp2home);
        params.times[i].set_home_to_hospitalized(thome2hosp);
        params.times[i].set_hospitalized_to_icu(thosp2icu);
        params.times[i].set_icu_to_home(ticu2home);
        params.times[i].set_infectious_asymp(tinfasy);
        params.times[i].set_icu_to_death(ticu2death);

        params.populations.set({i, mio::SecirCompartments::E}, fact * nb_exp_t0);
        params.populations.set({i, mio::SecirCompartments::C}, fact * nb_car_t0);
        params.populations.set({i, mio::SecirCompartments::I}, fact * nb_inf_t0);
        params.populations.set({i, mio::SecirCompartments::H}, fact * nb_hosp_t0);
        params.populations.set({i, mio::SecirCompartments::U}, fact * nb_icu_t0);
        params.populations.set({i, mio::SecirCompartments::R}, fact * nb_rec_t0);
        params.populations.set({i, mio::SecirCompartments::D}, fact * nb_dead_t0);
        params.populations.set_difference_from_group_total({i, mio::SecirCompartments::S}, mio::SecirCategory::AgeGroup,
                                                           i, fact * nb_total_t0);

        params.probabilities[i].set_infection_from_contact(0.98);
        params.probabilities[i].set_carrier_infectability(0.67);
        params.probabilities[i].set_asymp_per_infectious(alpha);
        params.probabilities[i].set_risk_from_symptomatic(beta);
        params.probabilities[i].set_hospitalized_per_infectious(rho);
        params.probabilities[i].set_icu_per_hospitalized(theta);
        params.probabilities[i].set_dead_per_icu(delta);
    }

    mio::ContactMatrixGroup& contact_matrix = params.get_contact_patterns();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(nb_groups, nb_groups, fact * cont_freq));
    contact_matrix.add_damping(0.3, mio::SimulationTime(30.));

    auto result_from_sim = simulate(t0, tmax, dt, params);

    mio::save_result(result_from_sim, "test_result.h5");
}