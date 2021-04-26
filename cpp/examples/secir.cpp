#include <epidemiology/secir/secir.h>
#include <epidemiology/model/simulation.h>
#include <epidemiology/utils/logging.h>

int main()
{
    epi::set_log_level(epi::LogLevel::debug);

    double t0   = 0;
    double tmax = 50;
    double dt   = 0.1;

    epi::log_info("Simulating SECIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    // working_params
    double tinc    = 5.2, // R_2^(-1)+R_3^(-1)
        tinfmild   = 6, // 4-14  (=R4^(-1))
        tserint    = 4.2, // 4-4.4 // R_2^(-1)+0.5*R_3^(-1)
        thosp2home = 12, // 7-16 (=R5^(-1))
        thome2hosp = 5, // 2.5-7 (=R6^(-1))
        thosp2icu  = 2, // 1-3.5 (=R7^(-1))
        ticu2home  = 8, // 5-16 (=R8^(-1))
        // tinfasy    = 6.2, // (=R9^(-1)=R_3^(-1)+0.5*R_4^(-1))
        ticu2death = 5; // 3.5-7 (=R5^(-1))

    double cont_freq = 10, // see Polymod study
        inf_prob = 0.05, carr_infec = 1,
           alpha = 0.09, // 0.01-0.16
        beta     = 0.25, // 0.05-0.5
        delta    = 0.3, // 0.15-0.77
        rho      = 0.2, // 0.1-0.35
        theta    = 0.25; // 0.15-0.4

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    epi::SecirModel model(1);

    // params.set_icu_capacity(20);
    model.parameters.set_start_day(0);
    model.parameters.set_seasonality(0);

    model.parameters.times[0].set_incubation(tinc);
    model.parameters.times[0].set_infectious_mild(tinfmild);
    model.parameters.times[0].set_serialinterval(tserint);
    model.parameters.times[0].set_hospitalized_to_home(thosp2home);
    model.parameters.times[0].set_home_to_hospitalized(thome2hosp);
    model.parameters.times[0].set_hospitalized_to_icu(thosp2icu);
    model.parameters.times[0].set_icu_to_home(ticu2home);
    model.parameters.times[0].set_icu_to_death(ticu2death);

    epi::ContactMatrixGroup& contact_matrix = model.parameters.get_contact_patterns();
    contact_matrix[0]                       = epi::ContactMatrix(Eigen::MatrixXd::Constant(1, 1, cont_freq));
    contact_matrix[0].add_damping(0.7, epi::SimulationTime(30.));

    model.populations.set_total(nb_total_t0);
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Exposed}] = nb_exp_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Carrier}] = nb_car_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Infected}] = nb_inf_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Hospitalized}] = nb_hosp_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::ICU}] = nb_icu_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Recovered}] = nb_rec_t0;
    model.populations[{epi::AgeGroup(0), epi::InfectionState::Dead}] = nb_dead_t0;
    model.populations.set_difference_from_total({epi::AgeGroup(0), epi::InfectionState::Susceptible}, nb_total_t0);

    model.parameters.probabilities[0].set_infection_from_contact(inf_prob);
    model.parameters.probabilities[0].set_carrier_infectability(carr_infec);
    model.parameters.probabilities[0].set_asymp_per_infectious(alpha);
    model.parameters.probabilities[0].set_risk_from_symptomatic(beta);
    model.parameters.probabilities[0].set_hospitalized_per_infectious(rho);
    model.parameters.probabilities[0].set_icu_per_hospitalized(theta);
    model.parameters.probabilities[0].set_dead_per_icu(delta);

    model.apply_constraints();

    epi::TimeSeries<double> secir = simulate(t0, tmax, dt, model);

    bool print_to_terminal = true;

    if (print_to_terminal) {
        char vars[] = {'S', 'E', 'C', 'I', 'H', 'U', 'R', 'D'};
        printf("\n # t");
        for (size_t k = 0; k < (size_t)epi::InfectionState::Count; k++) {
            printf(" %c", vars[k]);
        }
        auto num_points = static_cast<size_t>(secir.get_num_time_points());
        for (size_t i = 0; i < num_points; i++) {
            printf("\n%.14f ", secir.get_time(i));
            Eigen::VectorXd res_j = secir.get_value(i);
            for (size_t j = 0; j < (size_t)epi::InfectionState::Count; j++) {
                printf(" %.14f", res_j[j]);
            }
        }

        Eigen::VectorXd res_j = secir.get_last_value();
        printf("number total: %f",
               res_j[0] + res_j[1] + res_j[2] + res_j[3] + res_j[4] + res_j[5] + res_j[6] + res_j[7]);
    }
}
