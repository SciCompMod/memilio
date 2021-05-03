#include "load_test_data.h"
#include "epidemiology/model/simulation.h"
#include "epidemiology/secir/secir.h"
#include "epidemiology/utils/time_series.h"
#include <epidemiology_io/secir_result_io.h>
#include <gtest/gtest.h>

TEST(TestSaveResult, compareResultWithH5)
{
    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 10, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    epi::SecirModel model(1);
    auto& params     = model.parameters;
    size_t nb_groups = params.get_num_groups();;

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

        model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}] = nb_exp_t0;
        model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}] = nb_car_t0;
        model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}] = nb_inf_t0;
        model.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}] = nb_hosp_t0;
        model.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}] = nb_icu_t0;
        model.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}] = nb_rec_t0;
        model.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}] = nb_dead_t0;
        model.populations.set_difference_from_total({epi::AgeGroup(0), epi::InfectionState::Susceptible}, nb_total_t0);

        params.probabilities[i].set_infection_from_contact(0.06);
        params.probabilities[i].set_carrier_infectability(0.67);
        params.probabilities[i].set_asymp_per_infectious(alpha);
        params.probabilities[i].set_risk_from_symptomatic(beta);
        params.probabilities[i].set_hospitalized_per_infectious(rho);
        params.probabilities[i].set_icu_per_hospitalized(theta);
        params.probabilities[i].set_dead_per_icu(delta);
    }

    epi::ContactMatrixGroup& contact_matrix = params.get_contact_patterns();
    contact_matrix[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant(nb_groups, nb_groups, cont_freq));
    contact_matrix[0].add_damping(0.7, epi::SimulationTime(30.));

    auto result_from_sim                                  = simulate(t0, tmax, dt, model);
    std::vector<epi::TimeSeries<double>> results_from_sim = {result_from_sim, result_from_sim};
    std::vector<int> ids                                  = {1, 2};
    epi::save_result(results_from_sim, ids, "test_result.h5");

    std::vector<epi::SecirSimulationResult> results_from_file{
        epi::read_result("test_result.h5", static_cast<int>(nb_groups))};
    auto result_from_file = results_from_file[0];

    ASSERT_EQ(result_from_file.get_groups().get_num_time_points(), result_from_sim.get_num_time_points());
    ASSERT_EQ(result_from_file.get_totals().get_num_time_points(), result_from_sim.get_num_time_points());
    for (Eigen::Index i = 0; i < result_from_sim.get_num_time_points(); i++) {
        ASSERT_EQ(result_from_file.get_groups().get_num_elements(), result_from_sim.get_num_elements())
            << "at row " << i;
        ASSERT_EQ(result_from_file.get_totals().get_num_elements(),
                  result_from_sim.get_num_elements() / static_cast<Eigen::Index>(nb_groups))
            << "at row " << i;
        ASSERT_NEAR(result_from_sim.get_time(i), result_from_file.get_groups().get_time(i), 1e-10) << "at row " << i;
        ASSERT_NEAR(result_from_sim.get_time(i), result_from_file.get_totals().get_time(i), 1e-10) << "at row " << i;
        for (Eigen::Index l = 0; l < result_from_file.get_totals().get_num_elements(); l++) {
            double total = 0.0;
            for (Eigen::Index j = 0; j < Eigen::Index(nb_groups); j++) {
                total += result_from_sim[i][j * (size_t)epi::InfectionState::Count + l];
                EXPECT_NEAR(result_from_file.get_groups()[i][j * (size_t)epi::InfectionState::Count + l],
                            result_from_sim[i][j * (size_t)epi::InfectionState::Count + l], 1e-10)
                    << " at row " << i << " at row " << l << " at Group " << j;
            }
            EXPECT_NEAR(result_from_file.get_totals()[i][l], total, 1e-10) << " at row " << i << " at row " << l;
        }
    }
}
