#pragma once

#include <vector>

#include "../optimization_settings/secir_optimization.h"
#include "models/ode_secir/model.h"

#include "../control_parameters/damping_controls.h"
#include "../constraints/update_constraints.h"

#include "../helpers/integrator_selector.h"
#include "../helpers/make_time_grid.h"

#include "../constraints/infection_state_utils.h"

//     // assert(n == settings.num_dynamic_NPI_parameters() * graph_model.nodes().size());

//     // const size_t num_graph_nodes = graph_model.nodes().size();

//     FP objective = 0.0;

//     // std::vector<std::pair<FP, FP>> dynamic_NPI_values(settings.num_dynamic_NPI_parameters());
//     // std::vector<FP> commuter_testing_values(num_graph_nodes);

//     // for (const auto& dynamic_NPI : settings.dynamic_NPI_parameters()) {
//     //     size_t index                     = static_cast<size_t>(string_to_control(dynamic_NPI.name()));
//     //     dynamic_NPI_values[index].first  = ptr_parameters[2 * index + 0]; // threshold
//     //     dynamic_NPI_values[index].second = ptr_parameters[2 * index + 1]; // strength
//     // }
//     // for (size_t graph_node = 0; graph_node < num_graph_nodes; graph_node++) {
//     //     commuter_testing_values[graph_node] = ptr_parameters[2 * settings.num_dynamic_NPI_parameters() + graph_node];
//     // }

//     // set_dynamic_NPIs<FP>(settings, graph_model, dynamic_NPI_values);
//     // set_commuter_testing<FP>(settings, graph_model, commuter_testing_values);

//     // for (auto& node : graph_model.nodes()) {
//     //     node.property.get_simulation().set_integrator_core(
//     //         std::move(make_integrator<FP>(settings.integrator_type(), settings.dt())));
//     // }

//     // std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.simulation_days());
//     // auto graph_sim_mobility    = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), std::move(graph_model));

//     // for (size_t day = 0; day < settings.simulation_days(); day++) {

//     //     // for (const auto& node_entry : graph_model.nodes()) {
//     //     //     const auto& node_simulation = node_entry.property.get_simulation();
//     //     //     const auto& node_model      = node_simulation.get_model();
//     //     //     const auto& result          = node_simulation.get_result();

//     //     //     const auto num_age_groups      = node_model.parameters.get_num_groups();
//     //     //     const FP commuter_nondetection = node_model.parameters.get_commuter_nondetection();
//     //     //     const FP start_detection       = node_model.parameters.get_start_commuter_detection();
//     //     //     const FP end_detection         = node_model.parameters.get_end_commuter_detection();

//     //     //     const size_t num_time_points = result.get_num_time_points();

//     //     graph_sim_mobility.advance(time_steps[day + 1]);
//     // }

//     return objective;
//     // // ------------------------------------------------------------------ //
//     // // Evaluate the objective function of the model.                      //
//     // // Step 1. Define dampings based on 'const FP* parameters'.           //
//     // // Step 2. Evaluate the objective function based on                   //
//     // //         the parameters and the infection states in the simulation. //
//     // // ------------------------------------------------------------------ //
//     // assert(n == settings.num_control_parameters() * settings.num_control_intervals());
//     // const size_t num_graph_nodes = graph_model.nodes().size();

//     // FP objective = 0.0;

//     // std::vector<std::pair<FP, FP>> dynamic_NPI_values(settings.num_dynamic_NPI_parameters());
//     // std::vector<FP> commuter_testing_values(num_graph_nodes);

//     // for (const auto& dynamic_NPI : settings.dynamic_NPI_parameters()) {
//     //     size_t index                     = static_cast<size_t>(string_to_control(dynamic_NPI.name()));
//     //     dynamic_NPI_values[index].first  = ptr_parameters[2 * index + 0]; // threshold
//     //     dynamic_NPI_values[index].second = ptr_parameters[2 * index + 1]; // strength
//     // }
//     // for (size_t graph_node = 0; graph_node < num_graph_nodes; graph_node++) {
//     //     commuter_testing_values[graph_node] = ptr_parameters[2 * settings.num_dynamic_NPI_parameters() + graph_node];
//     // }

//     // set_dynamic_NPIs<FP>(settings, graph_model, dynamic_NPI_values);
//     // set_commuter_testing<FP>(settings, graph_model, commuter_testing_values);

//     // for (auto& node : graph_model.nodes()) {
//     //     node.property.get_simulation().set_integrator_core(
//     //         std::move(make_integrator<FP>(settings.integrator_type(), settings.dt())));
//     // }

//     // std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.simulation_days());
//     // auto graph_sim_mobility    = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), std::move(graph_model));

//     // for (size_t day = 0; day < settings.simulation_days(); day++) {
//     //     graph_sim_mobility.advance(time_steps[day + 1]);

//     //     auto& graph = graph_sim_mobility.get_graph();
//     //     for (size_t node_idx = 0; node_idx < graph.nodes().size(); node_idx++) {
//     //         const mio::SimulationNode<FP, mio::osecir::Simulation<FP>>& node = graph.nodes()[node_idx].property;
//     //         const mio::osecir::Simulation<FP>& node_simulation               = node.get_simulation();
//     //         const mio::osecir::Model<FP>& node_model                         = node_simulation.get_model();

//     //         const mio::TimeSeries<FP>& simulation_result      = node_simulation.get_result();
//     //         Eigen::Ref<const Eigen::VectorX<FP>> latest_state = simulation_result.get_last_value();

//     //         const size_t num_age_groups = node_model.parameters.get_num_groups();
//     //         const size_t num_states     = num_infection_states();

//     //         const auto infected_symptoms_states = query_infection_states("InfectedSymptoms");
//     //         const auto severe_states            = query_infection_states("Severe");
//     //         const auto critical_states          = query_infection_states("Critical");
//     //         const auto Dead                     = query_infection_states("Dead");

//     //         string infected_symptoms = "InfectedSymptoms"

//     //             {"Susceptible", InfectionState::Susceptible},
//     //                {"Exposed", InfectionState::Exposed}, {"InfectedNoSymptoms", InfectionState::InfectedNoSymptoms},
//     //                {"InfectedNoSymptomsConfirmed", InfectionState::InfectedNoSymptomsConfirmed},
//     //                {"InfectedSymptoms", InfectionState::InfectedSymptoms},
//     //                {"InfectedSymptomsConfirmed", InfectionState::InfectedSymptomsConfirmed},
//     //                {"InfectedSevere", InfectionState::InfectedSevere},
//     //                {"InfectedCritical", InfectionState::InfectedCritical}, {"Recovered", InfectionState::Recovered},
//     //         {
//     //             "Dead", InfectionState::Dead
//     //         }
//     //     };

//     //     // using std::max;

//     //     // for (size_t constraint_idx = 0; constraint_idx < optimization_settings.num_path_constraints(); constraint_idx++) {
//     //     //     const Constraint& constraint = optimization_settings.path_constraints()[constraint_idx];
//     //     //     const auto relevant_states   = query_infection_states(constraint.name());

//     //     //     FP total_value = 0.0;
//     //     //     for (size_t age_group = 0; age_group < num_age_groups; age_group++) {
//     //     //         for (const auto& state : relevant_states) {
//     //     //             size_t index = age_group * num_states + static_cast<size_t>(state);
//     //     //             total_value += latest_state[index];
//     //     //         }
//     //     //     }

//     //     //     path_constraints[constraint_idx] = max<FP>(path_constraints[constraint_idx], total_value);
//     //     // }

//     //     //         node_model.
//     // }

//     // Only add cost if exceeded threshold and was active
// }

// // interval_cost += cost("SchoolClosure") * param_at("SchoolClosure", node_index);
// // interval_cost += cost("HomeOffice") * param_at("HomeOffice", node_index);
// // interval_cost += cost("PhysicalDistancingSchool") * param_at("PhysicalDistancingSchool", node_index);
// // interval_cost += cost("PhysicalDistancingWork") * param_at("PhysicalDistancingWork", node_index);
// // interval_cost += cost("PhysicalDistancingOther") * param_at("PhysicalDistancingOther", node_index);

// // std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.num_intervals());
// // auto graph_sim_mobility    = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), graph_model);

// //         graph_sim_mobility.advance(time_steps[day + 1]);

// // for (size_t interval = 0; interval < settings.num_intervals(); interval++) {

// //     graph_sim_mobility.advance(time_steps[interval + 1]);

// //     size_t control_interval = interval / settings.pc_resolution();

// //     auto param_at = [&](const std::string& name, size_t nodex_index) {
// //         size_t control_index = static_cast<size_t>(string_to_control(name));
// //         return parameters[nodex_index][control_index + control_interval * settings.num_control_parameters()];
// //     };

// //     auto cost = [&](const std::string& name) {
// //         size_t control_index = static_cast<size_t>(string_to_control(name));
// //         return settings.control_parameters()[control_index].cost();
// //     };

// //     FP interval_cost = 0.0;
// //     for (size_t node_index = 0; node_index < num_graph_nodes; node_index++) {
// //         interval_cost += cost("SchoolClosure") * param_at("SchoolClosure", node_index);
// //         interval_cost += cost("HomeOffice") * param_at("HomeOffice", node_index);
// //         interval_cost += cost("PhysicalDistancingSchool") * param_at("PhysicalDistancingSchool", node_index);
// //         interval_cost += cost("PhysicalDistancingWork") * param_at("PhysicalDistancingWork", node_index);
// //         interval_cost += cost("PhysicalDistancingOther") * param_at("PhysicalDistancingOther", node_index);
// //     }

// //     objective += interval_cost / settings.num_intervals();
// // }

// // return objective;
// // }
// template <typename FP>
// FP objective_function(
//     mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>> graph_model,
//     const SecirOptimization& settings, const FP* ptr_parameters, size_t n)
// {

template <typename FP>
FP objective_function(
    mio::Graph<mio::SimulationNode<FP, mio::osecir::Simulation<FP>>, mio::MobilityEdge<FP>> graph_model,
    const SecirOptimization& settings, const FP* ptr_parameters, size_t n)
{
    const size_t num_states = 400;

    assert(n == 2 * settings.num_dynamic_NPI_parameters() + num_states);

    std::vector<std::pair<FP, FP>> dynamic_NPI_values(settings.num_dynamic_NPI_parameters());
    std::vector<FP> commuter_testing_values(num_states);

    for (const auto& dynamic_NPI : settings.dynamic_NPI_parameters()) {
        size_t index                     = static_cast<size_t>(string_to_control(dynamic_NPI.name()));
        dynamic_NPI_values[index].first  = ptr_parameters[2 * index + 0];
        dynamic_NPI_values[index].second = ptr_parameters[2 * index + 1];
    }
    for (size_t state_index = 0; state_index < num_states; ++state_index) {
        commuter_testing_values[state_index] = ptr_parameters[2 * settings.num_dynamic_NPI_parameters() + state_index];
    }

    // set_dynamic_NPIs<FP>(settings, graph_model, dynamic_NPI_values);
    set_commuter_testing<FP>(settings, graph_model, commuter_testing_values);

    for (auto& node : graph_model.nodes()) {
        node.property.get_simulation().set_integrator_core(
            std::move(make_integrator<FP>(settings.integrator_type(), settings.dt())));
    }

    std::vector<FP> time_steps = make_time_grid<FP>(settings.t0(), settings.tmax(), settings.simulation_days());
    auto graph_sim_mobility    = mio::make_mobility_sim<FP>(settings.t0(), settings.dt(), std::move(graph_model));

    FP objective_value = 0.0;

    for (size_t day = 0; day < settings.simulation_days(); ++day) {
        graph_sim_mobility.advance(time_steps[day + 1]);

        for (auto& node : graph_sim_mobility.get_graph().nodes()) {
            auto& sim          = node.property.get_simulation();
            auto& model        = sim.get_model();
            auto& dynamic_npis = model.parameters.template get<mio::osecir::DynamicNPIsInfectedSymptoms<FP>>();

            // --- Compute infections relative per 100k directly ---
            const auto& y = sim.get_result().get_last_value();
            FP sum_inf    = 0.0;
            for (mio::AgeGroup i = 0; i < model.parameters.get_num_groups(); ++i) {
                sum_inf += model.populations.get_from(y, {i, mio::osecir::InfectionState::InfectedSymptoms});
            }
            FP inf_rel = (sum_inf / model.populations.get_total()) * dynamic_npis.get_base_value();

            objective_value += sum_inf;

            // std::cout << inf_rel << " " << dynamic_NPI_values[0].first << " " << dynamic_NPI_values[0].second << "\n";

            // for (size_t dynamic_npi_index = 0; dynamic_npi_index < settings.num_dynamic_NPI_parameters();
            //      ++dynamic_npi_index) {
            //     const auto& [threshold, strength] = dynamic_NPI_values[dynamic_npi_index];
            //     if (inf_rel >= threshold) {
            //         // Strength multiplied by the number of people affected
            //         objective_value += strength * model.populations.get_total();
            //     }
            // }
        }
    }

    return objective_value;
}
