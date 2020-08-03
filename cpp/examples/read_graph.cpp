#include <epidemiology/migration.h>
#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology/seir.h>
#include <epidemiology/secir.h>
#include <iostream>

int main(int argc, char** argv)
{
    const auto t0   = 0.;
    const auto tmax = 10.;
    const auto dt   = 1.; //time step of migration, not integration

    /* std::vector<epi::SecirParams> params(1);
    epi::ContactFrequencyMatrix cont_freq(1);
    params[0].populations.set_exposed_t0(0);
    params[0].populations.set_infectious_t0(0);
    params[0].populations.set_total_t0(10000);
    params[0].populations.set_recovered_t0(0);
    params[0].times.set_incubation(1);
    //params.times.set_cont_freq(2.7);
    //params.times.set_infectious(1);

    //two mostly identical groups
    auto params_group1 = params;
    auto params_group2 = params;
    //some contact restrictions in group 1
    //params_group1.dampings.add({5, 0.5});
    //infection starts in group 1
    params_group1[0].populations.set_exposed_t0(10);*/

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

    int nb_groups = 1;
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

        params[i].probabilities.set_infection_from_contact(1.0);
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
        }
    }

    epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge> graph;
    graph.add_node(contact_freq_matrix, params, t0);
    graph.add_node(contact_freq_matrix, params, t0);
    graph.add_edge(0, 1, Eigen::VectorXd::Constant(8 * nb_groups, 0.01));
    graph.add_edge(1, 0, Eigen::VectorXd::Constant(8 * nb_groups, 0.01));

    epi::write_graph(graph, t0, tmax);

    epi::Graph<epi::ModelNode<epi::SecirSimulation>, epi::MigrationEdge> g = epi::read_graph();

    epi::print_graph(std::cout, g);
    auto node1 = g.nodes()[0].model.get_params();

    auto edges = g.edges()[0].property.coefficients;

    for (int i = 0; i < edges.size(); i++) {
        std::cout << g.nodes().size() << std::endl;
        std::cout << g.edges().size() << std::endl;
        std::cout << edges[i] << std::endl;
        std::cout << g.edges()[0].start_node_idx << std::endl;
        std::cout << g.edges()[0].end_node_idx << std::endl;
    }

    auto sim = epi::make_migration_sim(t0, dt, g);

    sim.advance(10);

    return 0;
}
