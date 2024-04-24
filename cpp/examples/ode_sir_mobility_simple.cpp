
#include "memilio/compartments/simulation.h"
#include "memilio/math/euler.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/io/mobility_io.h"
#include "ode_sir_mobility/infection_state.h"
#include "ode_sir_mobility/model.h"
#include "ode_sir_mobility/parameters.h"
#include "ode_sir_mobility/regions.h"
#include "ode_sir_mobility/contact_location.h"

int main()
{
    mio::set_log_level(mio::LogLevel::debug);

    double t0   = 0.;
    double tmax = 50.;
    double dt   = 1;

    size_t number_regions              = 2;
    size_t total_population_per_region = 10;

    mio::log_info("Simulating SIR; t={} ... {} with dt = {}.", t0, tmax, dt);

    mio::osirmobility::Model model(number_regions);

    for (size_t i = 0; i < number_regions; i++) {
        // model.populations[{mio::Index<mio::osirmobility::Region, mio::osirmobility::InfectionState>(
        //     mio::osirmobility::Region(i), mio::osirmobility::InfectionState::Infected)}]  = 1;
        model.populations[{mio::Index<mio::osirmobility::Region, mio::osirmobility::InfectionState>(
            mio::osirmobility::Region(1), mio::osirmobility::InfectionState::Infected)}]  = 5;
        model.populations[{mio::Index<mio::osirmobility::Region, mio::osirmobility::InfectionState>(
            mio::osirmobility::Region(0), mio::osirmobility::InfectionState::Infected)}]  = 0;
        model.populations[{mio::Index<mio::osirmobility::Region, mio::osirmobility::InfectionState>(
            mio::osirmobility::Region(i), mio::osirmobility::InfectionState::Recovered)}] = 0;
        model.populations[{mio::Index<mio::osirmobility::Region, mio::osirmobility::InfectionState>(
            mio::osirmobility::Region(i), mio::osirmobility::InfectionState::Susceptible)}] =
            total_population_per_region -
            model.populations[{mio::Index<mio::osirmobility::Region, mio::osirmobility::InfectionState>(
                mio::osirmobility::Region(i), mio::osirmobility::InfectionState::Infected)}] -
            model.populations[{mio::Index<mio::osirmobility::Region, mio::osirmobility::InfectionState>(
                mio::osirmobility::Region(i), mio::osirmobility::InfectionState::Recovered)}];
    }

    model.parameters.set<mio::osirmobility::TimeInfected>(2);
    model.parameters.set<mio::osirmobility::TransmissionProbabilityOnContact>(0.04);
    model.parameters.get<mio::osirmobility::ContactPatterns>().get_baseline()(0, 0) = 1.;
    // model.parameters.get<mio::osirmobility::ContactPatterns>().add_damping(0.6, mio::SimulationTime(12.5));
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(0), mio::osirmobility::Region(1), 0.});
    model.parameters.get<mio::osirmobility::CommutingRatio>().push_back(
        {mio::osirmobility::Region(1), mio::osirmobility::Region(0), 0.5});

    auto integrator = std::make_shared<mio::EulerIntegratorCore>();

    model.check_constraints();

    auto sir = simulate(t0, tmax, dt, model, integrator);

    bool print_to_terminal = true;

    sir.print_table();

    if (print_to_terminal) {

        std::vector<std::string> vars = {"S", "I", "R"};
        printf("Number of time points :%d\n", static_cast<int>(sir.get_num_time_points()));
        printf("People in\n");

        for (size_t k = 0; k < (size_t)mio::osirmobility::InfectionState::Count; k++) {
            double dummy = 0;

            for (size_t i = 0; i < (size_t)model.parameters.get_num_regions(); i++) {
                printf("\t %s[%d]: %.2f", vars[k].c_str(), (int)i,
                       sir.get_last_value()[k + (size_t)mio::osirmobility::InfectionState::Count * (int)i]);
                dummy += sir.get_last_value()[k + (size_t)mio::osirmobility::InfectionState::Count * (int)i];
            }

            printf("\t %s_total: %.2f\n", vars[k].c_str(), dummy);
        }
    }
}
