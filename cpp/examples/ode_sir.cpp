
#include "ode_sir/model.h"
#include "ode_sir/infection_state.h"
#include "ode_sir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0.;
    double tmax = 50.;
    double dt   = 0.1002004008016032;

    double total_population = 1061000;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::osir::Model model;

    model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Infected)}]  = 1000;
    model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Removed)}] = 1000;
    model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Removed)}];
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    model.parameters.set<mio::osir::TimeInfected>(2);
    model.parameters.set<mio::osir::TransmissionProbabilityOnContact>(1);
    model.parameters.get<mio::osir::ContactPatterns>().get_baseline()(0, 0) = 2.7;
    model.parameters.get<mio::osir::ContactPatterns>().add_damping(0.6, mio::SimulationTime(12.5));

    auto integrator                          = std::make_shared<mio::EulerIntegratorCore>();

    model.check_constraints();
    // print_seir_params(model);    

    auto sir = simulate(t0, tmax, dt, model, integrator);

    //printf("\n number total: %f\n",
    //       sir.get_last_value()[0] + sir.get_last_value()[1] + sir.get_last_value()[2]);

     bool print_to_terminal = true;


    if (print_to_terminal) {
        std::vector<std::string> vars = {"S", "I", "R"};
        printf("\n # t");
        for (size_t k = 0; k < (size_t)mio::osir::InfectionState::Count; k++) {
            printf(" %s", vars[k].c_str());
        }

        auto num_points = static_cast<size_t>(sir.get_num_time_points());
        for (size_t i = 0; i < num_points; i++) {
            printf("\n%.14f ", sir.get_time(i));
            Eigen::VectorXd res_j = sir.get_value(i);
            for (size_t j = 0; j < (size_t)mio::osir::InfectionState::Count; j++) {
                printf(" %.14f", res_j[j]);
            }
        }

        Eigen::VectorXd res_j = sir.get_last_value();
        printf("number total: %f",
               res_j[0] + res_j[1] + res_j[2]);
    }       
}
