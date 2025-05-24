#include "examples/state_estimators.h"
#include "ode_seir/model.h"
#include "memilio/compartments/flow_simulation.h"
#include <iostream>
#include <string>

using namespace mio::examples;

int main(int argc, char* argv[])
{

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <method>" << std::endl;
        std::cerr << "Available methods: euler, highorder, flow, prob" << std::endl;
        return 1;
    }
    std::string method = argv[1];

    // std::string method = "prob";
    // mio::unused(argc, argv);

    // Setup model
    mio::oseir::Model<ScalarType> model(1);
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Susceptible}] = 9700;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Exposed}]     = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Infected}]    = 100;
    model.populations[{mio::AgeGroup(0), mio::oseir::InfectionState::Recovered}]   = 100;
    model.parameters.set<mio::oseir::TimeExposed<ScalarType>>(5.2);
    model.parameters.set<mio::oseir::TimeInfected<ScalarType>>(6);
    model.parameters.set<mio::oseir::TransmissionProbabilityOnContact<ScalarType>>(0.1);
    model.check_constraints();

    // Create simulation
    ScalarType t0   = 0.0;
    ScalarType dt   = 0.5;
    ScalarType tmax = 100.0;
    auto sim        = mio::FlowSimulation<ScalarType, mio::oseir::Model<ScalarType>>(model, t0, dt);
    auto integrator = std::make_shared<mio::ExplicitStepperWrapper<ScalarType, boost::numeric::odeint::runge_kutta4>>();
    sim.set_integrator(integrator);
    sim.advance(tmax);
    const auto& seir_res = sim.get_result();

    // Create mobile population
    Eigen::VectorXd mobile_pop(4);
    mobile_pop << 2910, 30, 30, 30; // 30% of total population

    std::cout << "Running memory benchmark with method: " << method << std::endl;

    // Run the specified method
    for (ScalarType t = t0; t < tmax; t += dt) {
        const auto closest_idx = static_cast<size_t>(t / dt);
        const auto& total_pop  = seir_res.get_value(closest_idx);

        if (method == "euler") {
            integrate_mobile_population_euler(mobile_pop, sim, total_pop, t, dt);
        }
        else if (method == "highorder") {
            integrate_mobile_population_high_order(mobile_pop, sim, total_pop, t, dt);
        }
        else if (method == "flow") {
            flow_based_mobility_returns(mobile_pop, sim, total_pop, t, dt);
        }
        else if (method == "prob") {
            probabilistic_mobility_returns(mobile_pop, sim, total_pop, t, dt);
        }
        else {
            std::cerr << "Unknown method: " << method << std::endl;
            return 1;
        }
    }

    std::cout << "Memory benchmark completed successfully." << std::endl;
    return 0;
}
