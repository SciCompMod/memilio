#include "load_test_data.h"
#include "epidemiology/secir.h"
#include <epidemiology/save_result.h>
#include <epidemiology/save_parameters.h>
#include <gtest/gtest.h>

TEST(TestSaveResult, compareResultWithH5)
{
	double t0 = 0;
	double tmax = 50;
	double dt = 0.1;

	double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
		tinfasy = 6.2, ticu2death = 5;

	double cont_freq = 0.5, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

	double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
		nb_rec_t0 = 10, nb_dead_t0 = 0;

	int nb_groups = 1;
	double fact = 1.0 / (double)nb_groups;

	std::vector<epi::SecirParams> params{ epi::SecirParams{} };
	epi::ContactFrequencyMatrix contact_freq_matrix{ (size_t)nb_groups };
	for (size_t i = 1; i < nb_groups; i++) {
		params.push_back(epi::SecirParams{});
	}

	for (size_t i = 0; i < nb_groups; i++) {
		params[i].times.set_incubation(tinc);
		params[i].times.set_infectious_mild(tinfmild);
		params[i].times.set_serialinterval(tserint);
		params[i].times.set_hospitalized_to_home(thosp2home);
		params[i].times.set_home_to_hospitalized(thome2hosp);
		params[i].times.set_hospitalized_to_icu(thosp2icu);
		params[i].times.set_icu_to_home(ticu2home);
		params[i].times.set_infectious_asymp(tinfasy);
		params[i].times.set_icu_to_death(ticu2death);

		params[i].populations.set_total_t0(fact * nb_total_t0);
		params[i].populations.set_exposed_t0(fact * nb_exp_t0);
		params[i].populations.set_carrier_t0(fact * nb_car_t0);
		params[i].populations.set_infectious_t0(fact * nb_inf_t0);
		params[i].populations.set_hospital_t0(fact * nb_hosp_t0);
		params[i].populations.set_icu_t0(fact * nb_icu_t0);
		params[i].populations.set_recovered_t0(fact * nb_rec_t0);
		params[i].populations.set_dead_t0(fact * nb_dead_t0);

		params[i].probabilities.set_asymp_per_infectious(alpha);
		params[i].probabilities.set_risk_from_symptomatic(beta);
		params[i].probabilities.set_hospitalized_per_infectious(rho);
		params[i].probabilities.set_icu_per_hospitalized(theta);
		params[i].probabilities.set_dead_per_icu(delta);
	}

	epi::Damping dummy(30., 0.3);
	for (int i = 0; i < nb_groups; i++) {
		for (int j = i; j < nb_groups; j++) {
			contact_freq_matrix.set_cont_freq(fact * cont_freq, i, j);
			contact_freq_matrix.add_damping(dummy, i, j);
		}
	}


	std::vector<Eigen::VectorXd> secihurd(0);
	auto t = simulate(t0, tmax, dt, contact_freq_matrix, params, secihurd);

	save_result(t, secihurd, "test_result.h5");

	result test_result{ read_result("test_result.h5", nb_groups) };




	ASSERT_EQ(test_result.time.size(), t.size());
	ASSERT_EQ(test_result.groups[0].size(), secihurd.size());
	for (size_t i = 0; i < test_result.time.size(); i++) {
		ASSERT_EQ(test_result.groups[0][i].size(), secihurd[i].size() / nb_groups) << "at row " << i;
		ASSERT_NEAR(t[i], test_result.time[i], 1e-10) << "at row " << i;
		for (size_t l = 0; l < test_result.groups[0][i].size(); l++)
		{
			double dummy = 0.0;
			for (size_t j = 0; j < nb_groups; j++)
			{
				dummy += secihurd[i][j * 8 + l];
				EXPECT_NEAR(test_result.groups[j][i][l], secihurd[i][j*8 + l], 1e-10) << " at row " << i << " at row " << l << " at Group " << j;
			}
			EXPECT_NEAR(test_result.total[i][l], dummy, 1e-10) << " at row " << i << " at row " << l;

		}
	}
}
