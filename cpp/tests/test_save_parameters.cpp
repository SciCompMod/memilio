#include "load_test_data.h"
#include "epidemiology/secir.h"
#include <epidemiology_io/secir_result_io.h>
#include <epidemiology_io/secir_parameters_io.h>
//#include <epidemiology/parameter_studies/parameter_space.h>
//#include <epidemiology/parameter_studies/parameter_studies.h>
#include <gtest/gtest.h>

void check_dist(const epi::ParameterDistribution& dist, const epi::ParameterDistribution& dist_read)
{

    struct CheckDistEqVisitor : public epi::ConstParameterDistributionVisitor {
        CheckDistEqVisitor(const epi::ParameterDistribution& other_dist)
            : other(other_dist)
        {
        }

        void visit(const epi::ParameterDistributionNormal& self) override
        {
            auto p_other_normal_dist = dynamic_cast<const epi::ParameterDistributionNormal*>(&other);
            ASSERT_TRUE(p_other_normal_dist != nullptr);

            EXPECT_EQ(self.get_mean(), p_other_normal_dist->get_mean());
            EXPECT_EQ(self.get_standard_dev(), p_other_normal_dist->get_standard_dev());
            EXPECT_EQ(self.get_lower_bound(), p_other_normal_dist->get_lower_bound());
            EXPECT_EQ(self.get_upper_bound(), p_other_normal_dist->get_upper_bound());
        }
        void visit(const epi::ParameterDistributionUniform& self) override
        {
            auto p_other_uniform_dist = dynamic_cast<const epi::ParameterDistributionUniform*>(&other);
            ASSERT_TRUE(p_other_uniform_dist != nullptr);

            EXPECT_EQ(self.get_lower_bound(), p_other_uniform_dist->get_lower_bound());
            EXPECT_EQ(self.get_upper_bound(), p_other_uniform_dist->get_upper_bound());
        }
        const epi::ParameterDistribution& other;
    };

    CheckDistEqVisitor visitor(dist_read);
    dist.accept(visitor);
}

TEST(TestSaveParameters, compareSingleRun)
{
    double t0   = 0.0;
    double tmax = 50.5;
    double dt   = 0.1;

    double tinc = 5.2, tinfmild = 6, tserint = 4.2, thosp2home = 12, thome2hosp = 5, thosp2icu = 2, ticu2home = 8,
           tinfasy = 6.2, ticu2death = 5;

    double cont_freq = 0.5, alpha = 0.09, beta = 0.25, delta = 0.3, rho = 0.2, theta = 0.25;

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    int nb_groups = 2;
    double fact   = 1.0 / (double)nb_groups;

    std::vector<epi::SecirParams> params{epi::SecirParams{}};
    epi::ContactFrequencyMatrix contact_freq_matrix{(size_t)nb_groups};
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
    int nb_runs      = 5;
    std::string path = "/Parameters";
    TixiDocumentHandle handle;

    tixiCreateDocument("Parameters", &handle);
    epi::ParameterStudy study(epi::simulate, contact_freq_matrix, params, t0, tmax, 0.0, nb_runs);

    epi::write_parameter_study(handle, path, study);
    tixiSaveDocument(handle, "TestParameters.xml");
    tixiCloseDocument(handle);

    tixiOpenDocument("TestParameters.xml", &handle);
    epi::ParameterStudy read_study = epi::read_parameter_study(handle, path);
    tixiCloseDocument(handle);

    ASSERT_EQ(study.get_nb_runs(), read_study.get_nb_runs());
    ASSERT_EQ(study.get_t0(), read_study.get_t0());
    ASSERT_EQ(study.get_tmax(), read_study.get_tmax());

    const epi::ParameterSpace& space      = study.get_parameter_space();
    const epi::ParameterSpace& read_space = read_study.get_parameter_space();

    const epi::ContactFrequencyVariableElement& contact      = space.get_cont_freq_matrix_variable();
    const epi::ContactFrequencyVariableElement& read_contact = read_space.get_cont_freq_matrix_variable();

    int num_groups      = space.get_dead().size();
    int num_groups_read = read_space.get_dead().size();
    ASSERT_EQ(num_groups, num_groups_read);

    for (int i = 0; i < num_groups; i++) {
        ASSERT_EQ(space.get_dead()[i], read_space.get_dead()[i]);
        ASSERT_EQ(space.get_hospitalized()[i], read_space.get_hospitalized()[i]);
        ASSERT_EQ(space.get_icu()[i], read_space.get_icu()[i]);
        ASSERT_EQ(space.get_total()[i], read_space.get_total()[i]);

        check_dist(space.get_dist_exposed(i), read_space.get_dist_exposed(i));
        check_dist(space.get_dist_carrier(i), read_space.get_dist_carrier(i));
        check_dist(space.get_dist_recovered(i), read_space.get_dist_recovered(i));
        check_dist(space.get_dist_infectious(i), read_space.get_dist_infectious(i));

        check_dist(space.get_dist_incubation(i), read_space.get_dist_incubation(i));
        check_dist(space.get_dist_inf_mild(i), read_space.get_dist_inf_mild(i));
        check_dist(space.get_dist_serial_int_incub_diff(i), read_space.get_dist_serial_int_incub_diff(i));
        check_dist(space.get_dist_hosp_to_rec(i), read_space.get_dist_hosp_to_rec(i));
        check_dist(space.get_dist_inf_to_hosp(i), read_space.get_dist_inf_to_hosp(i));
        check_dist(space.get_dist_inf_asymp(i), read_space.get_dist_inf_asymp(i));
        check_dist(space.get_dist_hosp_to_icu(i), read_space.get_dist_hosp_to_icu(i));
        check_dist(space.get_dist_icu_to_rec(i), read_space.get_dist_icu_to_rec(i));
        check_dist(space.get_dist_icu_to_death(i), read_space.get_dist_icu_to_death(i));

        check_dist(space.get_dist_inf_from_cont(i), read_space.get_dist_inf_from_cont(i));
        check_dist(space.get_dist_risk_from_symp(i), read_space.get_dist_risk_from_symp(i));
        check_dist(space.get_dist_asymp_per_inf(i), read_space.get_dist_asymp_per_inf(i));
        check_dist(space.get_dist_death_per_icu(i), read_space.get_dist_death_per_icu(i));
        check_dist(space.get_dist_hosp_per_inf(i), read_space.get_dist_hosp_per_inf(i));
        check_dist(space.get_dist_icu_per_hosp(i), read_space.get_dist_icu_per_hosp(i));

        for (int j = 0; j < num_groups; j++) {
            ASSERT_EQ(contact.get_cont_freq().get_cont_freq(i, j), read_contact.get_cont_freq().get_cont_freq(i, j));
        }

        check_dist(contact.get_dist_num_dampings(), read_contact.get_dist_num_dampings());
        check_dist(contact.get_dist_day(), read_contact.get_dist_day());
        check_dist(contact.get_dist_damp_diag_base(), read_contact.get_dist_damp_diag_base());
        check_dist(contact.get_dist_damp_diag_rel(), read_contact.get_dist_damp_diag_rel());
        check_dist(contact.get_dist_damp_offdiag_rel(), read_contact.get_dist_damp_offdiag_rel());
    }
}
