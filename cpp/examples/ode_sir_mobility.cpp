
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "memilio/epidemiology/age_group.h"
#include "ode_sir_mobility/infection_state.h"
#include "ode_sir_mobility/model.h"
#include "ode_sir_mobility/parameters.h"
#include "ode_sir_mobility/regions.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0.;
    double tmax = 50.;
    double dt   = 0.1;

    double total_population = 1061000;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::osirmobility::Model model(5);

    model.populations[{mio::Index<mio::osirmobility::InfectionState, mio::osirmobility::Region>(mio::osirmobility::InfectionState::Infected, mio::osirmobility::Region(1))}]  = 1000;
    model.populations[{mio::Index<mio::osirmobility::InfectionState, mio::osirmobility::Region>(mio::osirmobility::InfectionState::Recovered, mio::osirmobility::Region(1))}] = 1000;
    model.populations[{mio::Index<mio::osirmobility::InfectionState, mio::osirmobility::Region>(mio::osirmobility::InfectionState::Susceptible, mio::osirmobility::Region(1))}] =
        total_population -
        model.populations[{mio::Index<mio::osirmobility::InfectionState, mio::osirmobility::Region>(mio::osirmobility::InfectionState::Infected, mio::osirmobility::Region(1))}] -
        model.populations[{mio::Index<mio::osirmobility::InfectionState, mio::osirmobility::Region>(mio::osirmobility::InfectionState::Recovered, mio::osirmobility::Region(1))}];
    model.parameters.set<mio::osirmobility::TimeInfected>(2);
    model.parameters.set<mio::osirmobility::TransmissionProbabilityOnContact>(0.04);
    model.parameters.get<mio::osirmobility::ContactPatterns>().get_baseline()(0, 0) = 2.7;
    model.parameters.get<mio::osirmobility::ContactPatterns>().add_damping(0.6, mio::SimulationTime(12.5));

    auto integrator = std::make_shared<mio::EulerIntegratorCore>();

    model.check_constraints();

    auto sir = simulate(t0, tmax, dt, model, integrator);

    bool print_to_terminal = true;

    if (print_to_terminal) {
        std::vector<std::string> vars = {"S", "I", "R"};
        printf("\n # t");
        for (size_t k = 0; k < (size_t)mio::osirmobility::InfectionState::Count; k++) {
            printf(" %s", vars[k].c_str());
        }

        auto num_points = static_cast<size_t>(sir.get_num_time_points());
        for (size_t i = 0; i < num_points; i++) {
            printf("\n%.14f ", sir.get_time(i));
            Eigen::VectorXd res_j = sir.get_value(i);
            for (size_t j = 0; j < (size_t)mio::osirmobility::InfectionState::Count; j++) {
                printf(" %.14f", res_j[j]);
            }
        }

        Eigen::VectorXd res_j = sir.get_last_value();
        printf("\nnumber total: %f \n", res_j[0] + res_j[1] + res_j[2]);
    }
}
