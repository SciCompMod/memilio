#include "ode_ui/model.h"
#include "ode_ui/infection_state.h"
#include "ode_ui/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 400;
    double dt   = 0.1;

    mio::log_info("Simulating UI; t={} ... {} with dt = {}.", t0, tmax, dt);

    

    //double nb_total_t0 = 10000;//  nb_inf_t0 = 100;
    const size_t num_groups = 10;

    mio::oui::Model<double> model(num_groups);
    double fact = 1.0 / num_groups;
    double cont_freq = num_groups;

    auto& params = model.parameters;

    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(num_groups); i++) {
        model.populations[{i, mio::oui::InfectionState::InfectedV1}]  = 0;//fact * nb_inf_t0;
        model.populations[{i, mio::oui::InfectionState::InfectedV2}]  = 0;//fact * nb_inf_t0;
        model.populations[{i, mio::oui::InfectionState::InfectedV3}]  = 0;//fact * nb_inf_t0;
        
        //model.populations.set_difference_from_group_total<mio::AgeGroup>({i, mio::oui::InfectionState::Susceptible},
        //                                                                 fact * (nb_total_t0 - 1000));
    }

    model.parameters.get<mio::oui::TransmissionProbabilityOnContactV1<double>>()        = 0.1;
    model.parameters.get<mio::oui::TransmissionProbabilityOnContactV2<double>>()        = 0.1;
    model.parameters.get<mio::oui::TransmissionProbabilityOnContactV3<double>>()        = 0.1;
    model.parameters.get<mio::oui::TimeInfectedV1<double>>()                            = 12;
    model.parameters.get<mio::oui::TimeInfectedV2<double>>()                            = 12 * 1.35;
    model.parameters.get<mio::oui::TimeInfectedV3<double>>()                            = 12 * 1.35 * 1.35;

    model.populations[{mio::AgeGroup(0), mio::oui::InfectionState::Susceptible}] = 9000;
    model.populations[{mio::AgeGroup(0), mio::oui::InfectionState::InfectedV1}]  = 1000;

    mio::ContactMatrixGroup& contact_matrix = params.get<mio::oui::ContactPatterns<double>>();
    contact_matrix[0] = mio::ContactMatrix(Eigen::MatrixXd::Constant(num_groups, num_groups, fact * cont_freq));
    //contact_matrix.add_damping(Eigen::MatrixXd::Constant(num_groups, num_groups, 0.7), mio::SimulationTime(30.));

    Eigen::MatrixXd rec_1 = Eigen::MatrixXd::Constant(num_groups, num_groups, 0.0);
    rec_1 << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 2, 0, 0, 100, 0, 0, 0, 0, 0,
             0, 0, 2, 0, 0, 100, 0, 0, 0, 0,
             0, 0, 0, 2, 0, 0, 100, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 100, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 100, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 100;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& recovery_matrix_v1 = params.get<mio::oui::RecoveryPatternsV1<double>>();
    //recovery_matrix_v1[0] = mio::ContactMatrix(Eigen::MatrixXd::Identity(num_groups, num_groups));
    recovery_matrix_v1 = rec_1;

    Eigen::MatrixXd rec_2 = Eigen::MatrixXd::Constant(num_groups, num_groups, 0.0);
    rec_2 << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 2, 0, 0, 2, 0, 100, 0, 0,
             0, 0, 0, 2, 0, 0, 2, 0, 100, 0, 
             0, 0, 0, 0, 0, 0, 0, 0, 0, 100;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& recovery_matrix_v2 = params.get<mio::oui::RecoveryPatternsV2<double>>();
    //recovery_matrix_v2[0] = mio::ContactMatrix(Eigen::MatrixXd::Identity(num_groups, num_groups));
    recovery_matrix_v2 = rec_2;

    Eigen::MatrixXd rec_3 = Eigen::MatrixXd::Constant(num_groups, num_groups, 0.0);
    rec_3 << 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 1, 0, 0, 1, 0, 0, 0, 0, 0,
             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
             0, 0, 1, 0, 0, 1, 0, 1, 0, 0,
             0, 0, 0, 2, 0, 0, 2, 0, 2, 100;

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>& recovery_matrix_v3 = params.get<mio::oui::RecoveryPatternsV3<double>>();
    //recovery_matrix_v3[0] = mio::ContactMatrix(Eigen::MatrixXd::Identity(num_groups, num_groups));
    recovery_matrix_v3 = rec_3;

    //model.apply_constraints();

    auto ui = simulate(t0, tmax, dt, model);

    std::vector<std::string> vars = {"U", "I1", "I2", "I3"};
    printf("Number of time points :%d\n", static_cast<int>(ui.get_num_time_points()));
    printf("People in\n");

    for (size_t k = 0; k < (size_t)mio::oui::InfectionState::Count; k++) {
        //double dummy = 0;

        for (size_t i = 0; i < (size_t)params.get_num_groups(); i++) {
            printf("\t %s[%d]: %.0f", vars[k].c_str(), (int)i,
                   ui.get_last_value()[k + (size_t)mio::oui::InfectionState::Count * (int)i]);
            //dummy += ui.get_last_value()[k + (size_t)mio::oui::InfectionState::Count * (int)i];
        }

        //printf("\t %s_Total: %.0f\n", vars[k].c_str(), dummy);
    }

    std::ofstream fileout;
    fileout.open("filename.txt");

    ui.print_table({"Susceptible_0", "InfectedV1_0", "InfectedV2_0", "InfectedV3_0", 
                    "Susceptible_1", "InfectedV1_1", "InfectedV2_1", "InfectedV3_1", 
                    "Susceptible_2", "InfectedV1_2", "InfectedV2_2", "InfectedV3_2", 
                    "Susceptible_3", "InfectedV1_3", "InfectedV2_3", "InfectedV3_3",
                    "Susceptible_4", "InfectedV1_4", "InfectedV2_4", "InfectedV3_4",
                    "Susceptible_5", "InfectedV1_5", "InfectedV2_5", "InfectedV3_5",
                    "Susceptible_6", "InfectedV1_6", "InfectedV2_6", "InfectedV3_6",
                    "Susceptible_7", "InfectedV1_7", "InfectedV2_7", "InfectedV3_7",
                    "Susceptible_8", "InfectedV1_8", "InfectedV2_8", "InfectedV3_8",
                    "Susceptible_9", "InfectedV1_9", "InfectedV2_9", "InfectedV3_9",
                    "Susceptible_10", "InfectedV1_10", "InfectedV2_10", "InfectedV3_10"}, 16, 5, fileout);
}
