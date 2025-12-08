#include "boost/filesystem.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "memilio/utils/logging.h"
#include "memilio/math/eigen.h"
#include "memilio/epidemiology/damping.h"

#include "optimization_model/optimization_model.h"

#include "helpers/integrator_selector.h"
#include "helpers/make_time_grid.h"
#include "constraints/infection_state_utils.h"

#include <iostream>
#include <vector>
#include <string>
#include <set>

int main(int argc, char* argv[])
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " [data directory] [optimal_control_settings directory]\n";
        return 1;
    }

    mio::set_log_level(mio::LogLevel::critical);

    boost::filesystem::path data_directory                     = argv[1];
    boost::filesystem::path optimal_control_settings_directory = argv[2];

    boost::filesystem::path constraints_file = optimal_control_settings_directory / "constraints.json";

    size_t num_age_groups = 1;

    size_t simulation_weeks = 8;
    size_t days_per_week    = 7;
    size_t simulation_days  = simulation_weeks * days_per_week;

    // General dynamic NPIs settings
    double npi_duration = 14.0;
    double npi_interval = 4.0;
    double npi_base     = 100'000; // Measured per 100'000 people

    OptimizationModel optimization_model(data_directory, optimal_control_settings_directory, simulation_days,
                                         num_age_groups, npi_duration, npi_interval, npi_base);

    auto graph_model = optimization_model.get_graph_model<double>();

    size_t num_control_intervals = simulation_weeks; // Number of control intervals (e.g., one for each week)
    size_t pc_resolution         = days_per_week; // Path constraint resolution (e.g., 7 days per control interval)

    IntegratorType integrator_type = IntegratorType::ControlledCashKarp54;
    size_t integrator_resolution   = 10;

    double t0(0.0);
    double tmax(simulation_days);

    // Number of intervals to check for path constraints (each day).
    size_t num_intervals = num_control_intervals * pc_resolution;
    // Time step size used in the integrator.
    double dt = (tmax - t0) / (num_intervals * integrator_resolution);

    try {
        // Set the integrator for each node simulation
        for (auto& node : graph_model.nodes()) {
            node.property.get_simulation().set_integrator_core(std::move(make_integrator<double>(integrator_type, dt)));
        }

        // Create graph simulation (move graph_model since parameter was passed by value)
        auto graph_simulation = mio::make_mobility_sim<double>(t0, dt, std::move(graph_model));

        // Time grid
        std::vector<double> time_steps = make_time_grid<double>(t0, tmax, num_intervals);

        // Prepare the infection-state queries we want to sum across all nodes
        const std::vector<std::string> state_names = {"Infected", "Severe", "Critical", "Dead"};
        std::vector<std::vector<mio::osecir::InfectionState>> queried_states;
        queried_states.reserve(state_names.size());
        for (const auto& name : state_names) {
            queried_states.push_back(query_infection_states(name));
        }

        // Open output CSV
        std::ofstream out("solution_states.csv");
        if (!out.is_open()) {
            std::cerr << "Failed to open 'solution_states.csv' for writing.\n";
        }
        else {
            // Header
            out << "time";
            for (const auto& name : state_names) {
                out << "," << name;
            }
            out << "\n";

            // Helper lambda to compute sums at current simulation state
            auto compute_sums = [&](std::vector<double>& sums) {
                std::fill(sums.begin(), sums.end(), static_cast<double>(0));
                // iterate nodes
                const auto& graph = graph_simulation.get_graph();
                for (size_t node_idx = 0; node_idx < graph.nodes().size(); ++node_idx) {
                    const auto& node            = graph.nodes()[node_idx].property;
                    const auto& node_simulation = node.get_simulation();
                    const auto& node_model      = node_simulation.get_model();
                    size_t num_age_groups       = static_cast<size_t>(node_model.parameters.get_num_groups());

                    const mio::TimeSeries<double>& simulation_result      = node_simulation.get_result();
                    Eigen::Ref<const Eigen::VectorX<double>> latest_state = simulation_result.get_last_value();

                    // For each queried state, sum over age groups and matching infection states
                    for (size_t q = 0; q < queried_states.size(); ++q) {
                        for (size_t age_group = 0; age_group < num_age_groups; ++age_group) {
                            for (auto st : queried_states[q]) {
                                size_t index = age_group * num_infection_states() + static_cast<size_t>(st);
                                sums[q] += latest_state[index];
                            }
                        }
                    }
                }
            };

            // Write initial time (t0)
            std::vector<double> sums(queried_states.size(), static_cast<double>(0));
            compute_sums(sums);
            out << std::fixed << std::setprecision(8) << time_steps[0];
            for (auto v : sums) {
                out << "," << v;
            }
            out << "\n";

            // Advance through time grid and write states after each advance
            for (size_t interval = 0; interval < num_intervals; ++interval) {
                // advance to next time point
                graph_simulation.advance(time_steps[interval + 1]);

                compute_sums(sums);

                out << std::fixed << std::setprecision(8) << time_steps[interval + 1];
                for (auto v : sums) {
                    out << "," << v;
                }
                out << "\n";
            }

            out.close();
            std::cout << "Wrote summed infection states to 'solution_states.csv'.\n";
        }
    }
    catch (const std::exception& e) {
        std::cerr << "Exception while saving solution: " << e.what() << "\n";
    }

    return 0;
}
