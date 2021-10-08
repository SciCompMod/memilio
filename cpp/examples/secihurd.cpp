#include "epidemiology/secir/secihurd.h"
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
    double tmax = 50;
    double dt   = 0.1;

    epi::log_info("Simulating SECIHURD_EXAMPLE; t={} ... {} with dt = {}.", t0, tmax, dt);

    gen::secihurd_example model; // create secihurd_example model instance

    // set secihurd_example model parameters        model.parameters.set<gen::StartDay>(0);
    model.parameters.set<gen::Seasonality>(0);
    model.parameters.set<gen::IncubationTime>(5.2);
    model.parameters.set<gen::InfectiousTimeMild>(6);
    model.parameters.set<gen::SerialInterval>(4.2);
    model.parameters.set<gen::HospitalizedToHomeTime>(1);
    model.parameters.set<gen::HomeToHospitalizedTime>(5);
    model.parameters.set<gen::HospitalizedToICUTime>(2);
    model.parameters.set<gen::ICUToHomeTime>(8);
    model.parameters.set<gen::ICUToDeathTime>(5);
    model.parameters.set<gen::InfectionProbabilityFromContact>(0.05);
    model.parameters.set<gen::RelativeCarrierInfectability>(1);
    model.parameters.set<gen::AsymptoticCasesPerInfectious>(0.09);
    model.parameters.set<gen::RiskOfInfectionFromSympomatic>(0.25);
    model.parameters.set<gen::HospitalizedCasesPerInfectious>(0.2);
    model.parameters.set<gen::ICUCasesPerHospitalized>(0.25);
    model.parameters.set<gen::DeathsPerHospitalized>(0.3);
    model.parameters.set<gen::TestAndTraceCapacity>(double(std::numeric_limits<double>::max()));
    model.parameters.set<gen::ICUCapacity>(0);
    model.parameters.set<gen::InfectiousTimeAsymptomatic>(5.0);

    // set starting population for each compartment        double pop_Susceptible=99760, pop_Exposed=100, pop_Carrier=50, pop_Infected=50, pop_Hospitalized=20, pop_ICU=10, pop_Recovered=10, pop_Dead=0;
    model.populations[{epi::Index<gen::Compartments>(gen::Compartments::Susceptible)}] = pop_Susceptible;
    model.populations[{epi::Index<gen::Compartments>(gen::Compartments::Exposed)}] = pop_Exposed;
    model.populations[{epi::Index<gen::Compartments>(gen::Compartments::Carrier)}] = pop_Carrier;
    model.populations[{epi::Index<gen::Compartments>(gen::Compartments::Infected)}] = pop_Infected;
    model.populations[{epi::Index<gen::Compartments>(gen::Compartments::Hospitalized)}] = pop_Hospitalized;
    model.populations[{epi::Index<gen::Compartments>(gen::Compartments::ICU)}] = pop_ICU;
    model.populations[{epi::Index<gen::Compartments>(gen::Compartments::Recovered)}] = pop_Recovered;
    model.populations[{epi::Index<gen::Compartments>(gen::Compartments::Dead)}] = pop_Dead;

    printf("Starting population: Susceptible = %f, Exposed = %f, Carrier = %f, Infected = %f, Hospitalized = %f, ICU = %f, Recovered = %f, Dead = %f\n", pop_Susceptible, pop_Exposed, pop_Carrier, pop_Infected, pop_Hospitalized, pop_ICU, pop_Recovered, pop_Dead);
    printf("Total: %f\n", pop_Susceptible + pop_Exposed + pop_Carrier + pop_Infected + pop_Hospitalized + pop_ICU + pop_Recovered + pop_Dead);

    // run secihurd_example model simulation
    epi::TimeSeries<double> secihurd_example = epi::simulate(t0, tmax, dt, model);

    print_to_terminal<double>(secihurd_example, std::vector<std::string>{"Susceptible", "Exposed", "Carrier", "Infected", "Hospitalized", "ICU", "Recovered", "Dead"});
}
