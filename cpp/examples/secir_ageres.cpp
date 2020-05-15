#include <epidemiology/secir.h>
#include <epidemiology/logging.h>

int main()
{
    epi::set_log_level(epi::LogLevel::debug);

    double t0   = 0;
    double tmax = 100;
    double dt   = 0.1;

    epi::log_info("Simulating SECIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    double tinc    = 5.2, // R_2^(-1)+R_3^(-1)
        tinfmild   = 6, // 4-14  (=R4^(-1))
        tserint    = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        thosp2home = 12, // 7-16 (=R5^(-1))
        thome2hosp = 5, // 2.5-7 (=R6^(-1))
        thosp2icu  = 2, // 1-3.5 (=R7^(-1))
        ticu2home  = 8, // 5-16 (=R8^(-1))
        tinfasy    = 6.2, // (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
        ticu2death = 5; // 3.5-7 (=R5^(-1))

    double tinfasy2 = 1.0 / (0.5 / (tinfmild - tserint) + 0.5 / tinfmild);
    if (fabs(tinfasy2 - tinfasy) > 0) {
        epi::log_warning("----> TODO / To consider: In the HZI paper, tinfasy (the asymptomatic infectious time) or "
                         "R9^(-1)=R_3^(-1)+0.5*R_4^(-1) is directly given by R_3 and R_4 and maybe should not be an "
                         "'additional parameter'");
    }

    double cont_freq = 0.5, // 0.2-0.75
        alpha        = 0.09, // 0.01-0.16
        beta         = 0.25, // 0.05-0.5
        delta        = 0.3, // 0.15-0.77
        rho          = 0.2, // 0.1-0.35
        theta        = 0.25; // 0.15-0.4

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    // alpha = alpha_in; // percentage of asymptomatic cases
    // beta  = beta_in; // risk of infection from the infected symptomatic patients
    // rho   = rho_in; // hospitalized per infected
    // theta = theta_in; // icu per hospitalized
    // delta = delta_in; // deaths per ICUs

    int nb_groups = 3;
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
    for (size_t i = 0; i < nb_groups; i++) {
        for (size_t j = i; j < nb_groups; j++) {
            contact_freq_matrix.set_cont_freq(fact * cont_freq, i, j);
            contact_freq_matrix.add_damping(dummy, i, j);
        }
    }

    print_secir_params(params);

    std::vector<Eigen::VectorXd> secir(0);

    simulate(t0, tmax, dt, contact_freq_matrix, params, secir);
    char vars[] = {'S', 'E', 'C', 'I', 'H', 'U', 'R', 'D'};
    printf("People in\n");
    for (size_t k = 0; k < 8; k++) {
        double dummy = 0;

        for (size_t i = 0; i < params.size(); i++) {
            printf("\t %c[%d]: %.0f", vars[k], (int)i, secir[secir.size() - 1][k]);
            dummy += secir[secir.size() - 1][k];
        }

        printf("\t %c_otal: %.0f\n", vars[k], dummy);
    }

    // printf("People in\n S[0]:\t %.0f\tS[1]:\t %.0f, S_total:\t %.0f", secir[secir.size() - 1][0],
    //        secir[secir.size() - 1][8], secir[secir.size() - 1][0] + secir[secir.size() - 1][8]);
}
