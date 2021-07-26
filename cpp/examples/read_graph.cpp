#include <epidemiology/migration/migration.h>
#include <epidemiology_io/secir_parameters_io.h>
#include <epidemiology_io/mobility_io.h>
#include <data_dir.h>
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

    std::string dir = DATA_DIR;

    auto filename   = epi::path_join(dir, "migration", (std::string)argv[1]);
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

    epi::SecirModel model(1);
    epi::AgeGroup nb_groups = model.parameters.get_num_groups();
    double fact   = 1.0 / (double)(size_t)nb_groups;

    auto& params = model.parameters;

    params.set<epi::ICUCapacity>(std::numeric_limits<double>::max());
    params.set<epi::StartDay>(0);
    params.set<epi::Seasonality>(0);

    for (auto i = epi::AgeGroup(0); i < nb_groups; i++) {
        params.get<epi::IncubationTime>()[i] = tinc;
        params.get<epi::InfectiousTimeMild>()[i] = tinfmild;
        params.get<epi::SerialInterval>()[i] = tserint;
        params.get<epi::HospitalizedToHomeTime>()[i] = thosp2home;
        params.get<epi::HomeToHospitalizedTime>()[i] = thome2hosp;
        params.get<epi::HospitalizedToICUTime>()[i] = thosp2icu;
        params.get<epi::ICUToHomeTime>()[i] = ticu2home;
        params.get<epi::ICUToDeathTime>()[i] = ticu2death;

        model.populations[{i, epi::InfectionState::Exposed}] = fact * nb_exp_t0;
        model.populations[{i, epi::InfectionState::Carrier}] = fact * nb_car_t0;
        model.populations[{i, epi::InfectionState::Infected}] = fact * nb_inf_t0;
        model.populations[{i, epi::InfectionState::Hospitalized}] = fact * nb_hosp_t0;
        model.populations[{i, epi::InfectionState::ICU}] = fact * nb_icu_t0;
        model.populations[{i, epi::InfectionState::Recovered}] = fact * nb_rec_t0;
        model.populations[{i, epi::InfectionState::Dead}] = fact * nb_dead_t0;
        model.populations.set_difference_from_group_total<epi::AgeGroup>({i, epi::InfectionState::Susceptible},
                                                                         fact * nb_total_t0);

        params.get<epi::InfectionProbabilityFromContact>()[i] = inf_prob;
        params.get<epi::RelativeCarrierInfectability>()[i] = carr_infec;
        params.get<epi::AsymptoticCasesPerInfectious>()[i] = alpha;
        params.get<epi::RiskOfInfectionFromSympomatic>()[i] = beta;
        params.get<epi::HospitalizedCasesPerInfectious>()[i] = rho;
        params.get<epi::ICUCasesPerHospitalized>()[i] = theta;
        params.get<epi::DeathsPerHospitalized>()[i] = delta;
    }

    params.apply_constraints();

    epi::ContactMatrixGroup& contact_matrix = params.get<epi::ContactPatterns>();
    contact_matrix[0] = epi::ContactMatrix(Eigen::MatrixXd::Constant((size_t)nb_groups, (size_t)nb_groups, fact * cont_freq));

    epi::set_params_distributions_normal(model, t0, tmax, 0.2);

    std::cout << "Reading Migration File..." << std::flush;
    auto read_mobility_result = epi::read_mobility_plain(filename);
    if (!read_mobility_result)
    {
        std::cout << read_mobility_result.error().formatted_message() << '\n';
        return -1;
    }
    auto& twitter_migration_2018 = read_mobility_result.value();
    std::cout << "Done" << std::endl;

    std::cout << "Intializing Graph..." << std::flush;
    epi::Graph<epi::SecirModel, epi::MigrationParameters> graph;
    for (int node = 0; node < twitter_migration_2018.rows(); node++) {
        graph.add_node(node, model);
    }
    for (int row = 0; row < twitter_migration_2018.rows(); row++) {
        for (int col = 0; col < twitter_migration_2018.cols(); col++) {
            graph.add_edge(
                row, col,
                Eigen::VectorXd::Constant(8 * (size_t)nb_groups, twitter_migration_2018(row, col) /
                                                             graph.nodes()[row].property.populations.get_total()));
        }
    }
    std::cout << "Done" << std::endl;

    std::cout << "Writing Json Files..." << std::flush;
    auto write_status = epi::write_graph(graph, "graph_parameters");
    if (!write_status) {
        std::cout << "Error: " << write_status.error().formatted_message();
    }
    std::cout << "Done" << std::endl;

    std::cout << "Reading Json Files..." << std::flush;
    auto graph_read_result = epi::read_graph<epi::SecirModel>("graph_parameters");
        
    if (!graph_read_result) {
        std::cout << "Error: " << graph_read_result.error().formatted_message();
    }
    std::cout << "Done" << std::endl;
    auto& graph_read = graph_read_result.value();

    std::cout << "Running Simulations..." << std::flush;
    auto study = epi::ParameterStudy<epi::SecirSimulation<>>(graph_read, t0, tmax, 0.5, 2);
    std::cout << "Done" << std::endl;

    return 0;
}
