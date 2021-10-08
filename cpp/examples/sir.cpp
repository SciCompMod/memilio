#include "epidemiology/secir/sir.h"
#include "epidemiology/model/simulation.h"

#include <iostream>
#include <vector>
#include <string>

template<class SC>
void print_to_terminal(const epi::TimeSeries<SC>& results, const std::vector<std::string>& state_names) {
    printf("| %-16s |", "Time");
    for (size_t k = 0; k < state_names.size(); k++) {
        printf(" %-16s |", state_names[k].data()); // print underlying char*
    }
    auto num_points = static_cast<size_t>(results.get_num_time_points());
    for (size_t i = 0; i < num_points; i++) {
        printf("\n| %16.6f |", results.get_time(i));
        auto res_i = results.get_value(i);
        for (size_t j = 0; j < state_names.size(); j++) {
            printf(" %16.6f |", res_i[j]);
        }
    }
    printf("\n");
}

int main() {
    epi::set_log_level(epi::LogLevel::debug);

    // set start time, stop time and time step length
    double t0   = 0;
    double tmax = 1;
    double dt   = 0.001;

    epi::log_info("Simulating SIR_EXAMPLE; t={} ... {} with dt = {}.", t0, tmax, dt);

    gen::sir_example model; // create sir_example model instance

    // set sir_example model parameters        model.parameters.set<gen::ContactsPerDay>(15);
    model.parameters.set<gen::RecoveryRate>(0.5);

    // set starting population for each compartment        double pop_S=9900, pop_I=100, pop_R=0;
    model.populations[{epi::Index<gen::Compartments>(gen::Compartments::S)}] = pop_S;
    model.populations[{epi::Index<gen::Compartments>(gen::Compartments::I)}] = pop_I;
    model.populations[{epi::Index<gen::Compartments>(gen::Compartments::R)}] = pop_R;

    printf("Starting population: S = %f, I = %f, R = %f\n", pop_S, pop_I, pop_R);
    printf("Total: %f\n", pop_S + pop_I + pop_R);

    // run sir_example model simulation
    epi::TimeSeries<double> sir_example = epi::simulate(t0, tmax, dt, model);

    print_to_terminal<double>(sir_example, std::vector<std::string>{"S", "I", "R"});
}
