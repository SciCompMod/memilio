
#include "ode_sir/model.h"
#include "ode_sir/infection_state.h"
#include "ode_sir/parameters.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0;
    double tmax = 1;
    double dt   = 0.001;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::osir::Model model;

    double total_population                                                                            = 10000;
    model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Infected)}]  = 100;
    model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Removed)}] = 100;
    model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Susceptible)}] =
        total_population -
        model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Infected)}] -
        model.populations[{mio::Index<mio::osir::InfectionState>(mio::osir::InfectionState::Removed)}];
    // suscetible now set with every other update
    // params.nb_sus_t0   = params.nb_total_t0 - params.nb_exp_t0 - params.nb_inf_t0 - params.nb_rec_t0;
    model.parameters.set<mio::osir::TimeInfected>(6);
    model.parameters.set<mio::osir::TransmissionProbabilityOnContact>(0.04);
    model.parameters.get<mio::osir::ContactPatterns>().get_baseline()(0, 0) = 10;

    model.check_constraints();
    // print_seir_params(model);    

    auto sir = simulate(t0, tmax, dt, model);

    printf("\n number total: %f\n",
           sir.get_last_value()[0] + sir.get_last_value()[1] + sir.get_last_value()[2]);
}
