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

    epi::SecirParams params{};

    params.times.set_incubation(tinc);
    params.times.set_infectious_mild(tinfmild);
    params.times.set_serialinterval(tserint);
    params.times.set_hospitalized_to_home(thosp2home);
    params.times.set_home_to_hospitalized(thome2hosp);
    params.times.set_hospitalized_to_icu(thosp2icu);
    params.times.set_icu_to_home(ticu2home);
    params.times.set_infectious_asymp(tinfasy);
    params.times.set_icu_to_death(ticu2death);
    params.times.set_cont_freq(cont_freq);

    params.populations.set_total_t0(nb_total_t0);
    params.populations.set_exposed_t0(nb_exp_t0);
    params.populations.set_carrier_t0(nb_car_t0);
    params.populations.set_infectious_t0(nb_inf_t0);
    params.populations.set_hospital_t0(nb_hosp_t0);
    params.populations.set_icu_t0(nb_icu_t0);
    params.populations.set_recovered_t0(nb_rec_t0);
    params.populations.set_dead_t0(nb_dead_t0);

    params.probabilities.set_asymp_per_infectious(alpha);
    params.probabilities.set_risk_from_symptomatic(beta);
    params.probabilities.set_hospitalized_per_infectious(rho);
    params.probabilities.set_icu_per_hospitalized(theta);
    params.probabilities.set_dead_per_icu(delta);

    params.dampings.add(epi::Damping(30., 0.3));

    print_secir_params(params);

    std::vector<std::vector<double>> secir(0);

    simulate(t0, tmax, dt, params, secir);

    printf("number total: %f", secir[secir.size() - 1][0] + secir[secir.size() - 1][1] + secir[secir.size() - 1][2] +
                                   secir[secir.size() - 1][3] + secir[secir.size() - 1][4] +
                                   secir[secir.size() - 1][5] + secir[secir.size() - 1][6] +
                                   secir[secir.size() - 1][7]);
}
