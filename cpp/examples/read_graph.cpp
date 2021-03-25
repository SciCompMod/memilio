#include <epidemiology/migration/migration.h>
#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology_io/mobility_io.h>

#include <iostream>

void print_usage()
{
    std::cout << "Usage: read_graph MIGRATION_FILE"
              << "\n\n";
    std::cout << "This example performs a simulation based on twitter "
                 "migration data."
              << std::endl;
}

int main(int argc, char** argv)
{
    if (argc != 2) {
        print_usage();
        return -1;
    }

    auto filename = argv[1];

    const auto t0   = 0.;
    const auto tmax = 10.;
    const auto dt   = 1.; //time step of migration, not integration

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
        inf_prob = 0.05, carr_infec = 0.67,
           alpha = 0.09, // 0.01-0.16
        beta     = 0.25, // 0.05-0.5
        delta    = 0.3, // 0.15-0.77
        rho      = 0.2, // 0.1-0.35
        theta    = 0.25; // 0.15-0.4

    double nb_total_t0 = 10000, nb_exp_t0 = 100, nb_inf_t0 = 50, nb_car_t0 = 50, nb_hosp_t0 = 20, nb_icu_t0 = 10,
           nb_rec_t0 = 10, nb_dead_t0 = 0;

    epi::SecirModel<epi::AgeGroup1> model;
    size_t nb_groups = model.parameters.get_num_groups();
    double fact   = 1.0 / (double)nb_groups;

    auto& params = model.parameters;

    params.set_icu_capacity(std::numeric_limits<double>::max());
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
        params.times[i].set_icu_to_death(ticu2death);

        model.populations[{(epi::AgeGroup1)i, epi::InfectionType::E}] = fact * nb_exp_t0;
        model.populations[{(epi::AgeGroup1)i, epi::InfectionType::C}] = fact * nb_car_t0;
        model.populations[{(epi::AgeGroup1)i, epi::InfectionType::I}] = fact * nb_inf_t0;
        model.populations[{(epi::AgeGroup1)i, epi::InfectionType::H}] = fact * nb_hosp_t0;
        model.populations[{(epi::AgeGroup1)i, epi::InfectionType::U}] = fact * nb_icu_t0;
        model.populations[{(epi::AgeGroup1)i, epi::InfectionType::R}] = fact * nb_rec_t0;
        model.populations[{(epi::AgeGroup1)i, epi::InfectionType::D}] = fact * nb_dead_t0;
        model.populations.set_difference_from_group_total(fact * nb_total_t0, (epi::AgeGroup1)i, (epi::AgeGroup1)i,
                                                          epi::InfectionType::S);

        params.probabilities[i].set_infection_from_contact(inf_prob);
        params.probabilities[i].set_carrier_infectability(carr_infec);
        params.probabilities[i].set_asymp_per_infectious(alpha);
        params.probabilities[i].set_risk_from_symptomatic(beta);
        params.probabilities[i].set_hospitalized_per_infectious(rho);
        params.probabilities[i].set_icu_per_hospitalized(theta);
        params.probabilities[i].set_dead_per_icu(delta);
    }

    params.apply_constraints();

    epi::ContactMatrixGroup& contact_matrix = params.get_contact_patterns();
    contact_matrix[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant(nb_groups, nb_groups, fact * cont_freq));

    epi::set_params_distributions_normal(model, t0, tmax, 0.2);

    std::cout << "Reading Migration File..." << std::flush;
    Eigen::MatrixXd twitter_migration_2018 = epi::read_mobility_plain(filename);
    std::cout << "Done" << std::endl;

    std::cout << "Intializing Graph..." << std::flush;
    epi::Graph<epi::SecirModel<epi::AgeGroup1>, epi::MigrationParameters> graph;
    for (int node = 0; node < twitter_migration_2018.rows(); node++) {
        graph.add_node(node, model);
    }
    for (int row = 0; row < twitter_migration_2018.rows(); row++) {
        for (int col = 0; col < twitter_migration_2018.cols(); col++) {
            graph.add_edge(row, col,
                           Eigen::VectorXd::Constant(8 * nb_groups, twitter_migration_2018(row, col) /
                                                                        graph.nodes()[row].property.populations.get_total()));
        }
    }
    std::cout << "Done" << std::endl;

    std::cout << "Writing XML Files..." << std::flush;
    epi::write_graph(graph, "graph_parameters");
    std::cout << "Done" << std::endl;

    std::cout << "Reading XML Files..." << std::flush;
    epi::Graph<epi::SecirModel<epi::AgeGroup1>, epi::MigrationParameters> graph_read =
        epi::read_graph<epi::SecirModel<epi::AgeGroup1>>("graph_parameters");

    std::cout << "Done" << std::endl;

    std::cout << "Running Simulations..." << std::flush;
    auto study = epi::ParameterStudy<epi::SecirModel<epi::AgeGroup1>>(graph_read, t0, tmax, 0.5, 2);
    std::cout << "Done" << std::endl;

    return 0;
}
