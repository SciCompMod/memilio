// clang-format off
#include "secirvvs_optimization.h"

#include <vector>
#include <string>
#include <tuple>
#include <utility>
#include <optional>
#include <cstddef>

#include "models/ode_secirvvs/model.h"
#include "memilio/compartments/simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/mobility/graph.h"

#include "integrator.h"
#include "ad_type.h"

#include "controls.h"

#include <iostream>
#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"

#include "initialize_simulation.h"
#include "secirvvs_ipopt.h"

#include "integrator.h"
#include "ad_type.h"

#include "time_grid.h"

#include "damping_helper.h"
#include "constraint_parsing.h"



void SecirvvsOptimization::update_path_constraint(auto graph_sim_mobility, std::vector<double>& path_constraint_values)
{
    assert(path_constraint_values.size() == numPathConstraints());

    for (size_t i = 0; i < numPathConstraints(); i++) {
        const ConstraintPolicy<double>& constraint = pathConstraints()[i];
        auto states = query_infection_states(constraint.get_name());
        std::optional<int> node_id = constraint.get_node_index();

        double value = 0.0;
        for (size_t node_idx = 0; node_idx < numGraphNodes(); node_idx++) {
            const auto& node  = graph_sim_mobility.get_graph().nodes()[node_idx];
            auto&       model = node.property.get_simulation().get_model();

            if (!node_id.has_value() || node_id.value() == node.id) {
                for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++){
                    for (auto state : states) {
                        value += model.populations[{age_group, state}];
                    }
                }
            }
        }
        path_constraint_values[i] = std::max(path_constraint_values[i], value);
    }
}


void SecirvvsOptimization::update_terminal_constraint(auto graph_sim_mobility, std::vector<double>& terminal_constraint_values)
{
    assert(terminal_constraint_values.size() == numTerminalConstraints());

    for (size_t i = 0; i < numTerminalConstraints(); i++) {
        const ConstraintPolicy<double>& constraint = terminalConstraints()[i];
        auto states = query_infection_states(constraint.get_name());
        std::optional<int> node_id = constraint.get_node_index();

        double value = 0.0;
        for (size_t node_idx = 0; node_idx < numGraphNodes(); node_idx++) {
            const auto& node  = graph_sim_mobility.get_graph().nodes()[node_idx];
            auto&       model = node.property.get_simulation().get_model();

            if (!node_id.has_value() || node_id.value() == node.id) {
                for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++){
                    for (auto state : states) {
                        value += model.populations[{age_group, state}];
                    }
                }
            }
        }
        terminal_constraint_values[i] += value;
    }
}




void SecirvvsOptimization::validateConstraints()
{
    std::vector<double> path_constraint_values(numPathConstraints(), 0.0);
    std::vector<double> terminal_constraint_values(numTerminalConstraints(), 0.0);

    auto graph_simulation = initialize_simulation<double>(t0(), tmax());
    set_restrictive_dampings(graph_simulation);

    for (auto& node : graph_simulation.nodes()){
        node.property.get_simulation().set_integrator(make_integrator<double>(integratorType(), dt()));
    }

    std::vector<double> timeSteps = make_time_grid<double>(t0(), tmax(), numIntervals());
    auto graph_sim_mobility = mio::make_mobility_sim<double>(t0(), dt(), graph_simulation);

    update_path_constraint(graph_sim_mobility, path_constraint_values);
    for (size_t interval = 0; interval < numIntervals(); interval++) {
        size_t controlIndex = interval / pcResolution();

        graph_sim_mobility.advance(timeSteps[interval+1]);

        for (auto& node : graph_sim_mobility.get_graph().nodes()) {
            auto& model = node.property.get_simulation().get_model();

            mio::TimeSeries<double> result = node.property.get_simulation().get_result();
            const auto& final_state = result.get_last_value();

            for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
                for (size_t state_index = 0; state_index < num_infection_states(); state_index++) {
                    size_t idx = age_group.get() * num_infection_states() + state_index;
                    model.populations[{age_group, InfectionState(state_index)}] = final_state[idx];
                }
            }
        }
        update_path_constraint(graph_sim_mobility, path_constraint_values);
    }
    update_terminal_constraint(graph_sim_mobility, terminal_constraint_values);


    // Validate and clamp path constraint values
    for (size_t i = 0; i < numPathConstraints(); ++i) {
        auto& constraint = pathConstraints()[i];
        double& value = path_constraint_values[i];
        double min_val = constraint.get_min();
        double max_val = constraint.get_max();
        if (value > max_val) {
            std::ostringstream msg;
            msg << "Path Constraint [" << i << "] \"" << constraint.get_name() << "\"";

            if (constraint.has_node_index()) {
                msg << " at node " << constraint.get_node_index().value();
            } else {
                msg << " (global)";
            }

            msg << ": Increasing "<< max_val << " to " << value << ".\n";

            constraint.range.second = value;

            std::cout << msg.str();
        }
    }

    // Validate and clamp terminal constraint values
    for (size_t i = 0; i < numTerminalConstraints(); ++i) {
        auto& constraint = terminalConstraints()[i];
        double& value = terminal_constraint_values[i];
        double min_val = constraint.get_min();
        double max_val = constraint.get_max();
        if (value > max_val) {
            std::ostringstream msg;
            msg << "Terminal Constraint [" << i << "] \"" << constraint.get_name() << "\"";

            if (constraint.has_node_index()) {
                msg << " at node " << constraint.get_node_index().value();
            } else {
                msg << " (global)";
            }

            msg << ": Increasing "<< max_val << " to " << value << ".\n";

            constraint.range.second = value;

            std::cout << msg.str();
        }
    }
}









void SecirvvsOptimization::set_restrictive_dampings(
    mio::Graph<mio::SimulationNode<double, mio::Simulation<double, mio::osecirvvs::Model<double>>>,
    mio::MobilityEdge<double>>& graph_simulation)
{
    std::vector<double> timeSteps = make_time_grid<double>(t0(), tmax(), numIntervals());

    for (auto& node : graph_simulation.nodes()) {
        auto& simulation = node.property.get_simulation();
        auto& model      = simulation.get_model();
        size_t num_age_groups = static_cast<size_t>(model.parameters.get_num_groups());

        mio::UncertainContactMatrix<double>& contacts =
            model.parameters.template get<mio::osecirvvs::ContactPatterns<double>>();
        std::vector<mio::DampingSampling<double>>& contact_dampings = contacts.get_dampings();
        contact_dampings.clear();

        for (size_t controlIndex = 0; controlIndex < numControlIntervals(); controlIndex++) {

            auto most_restrictive_control_value = [&](const std::string& name) {
                return m_control_policy[static_cast<size_t>(string_to_control(name))].get_max();
            };
            double school_closure             = most_restrictive_control_value("SchoolClosure");
            double home_office                = most_restrictive_control_value("HomeOffice");
            double physical_distancing_school = most_restrictive_control_value("PhysicalDistancingSchool");
            double physical_distancing_work   = most_restrictive_control_value("PhysicalDistancingWork");
            double physical_distancing_other  = most_restrictive_control_value("PhysicalDistancingOther");

            size_t timeStepIndex = controlIndex * pcResolution();
            mio::SimulationTime<double> time(timeSteps[timeStepIndex]);

            auto effectiveness = [&](const std::string& name) {
                return m_control_policy[static_cast<size_t>(string_to_control(name))].get_effectiveness();
            };

            contact_dampings.push_back(
                set_school_closure(time, effectiveness("SchoolClosure") * school_closure, num_age_groups));
            contact_dampings.push_back(
                set_home_office(time, effectiveness("HomeOffice") * home_office, num_age_groups));
            contact_dampings.push_back(
                set_physical_distancing_school(time, effectiveness("PhysicalDistancingSchool") * physical_distancing_school, num_age_groups));
            contact_dampings.push_back(
                set_physical_distancing_work(time, effectiveness("PhysicalDistancingWork") * physical_distancing_work, num_age_groups));
            contact_dampings.push_back(
                set_physical_distancing_other(time, effectiveness("PhysicalDistancingOther") * physical_distancing_other, num_age_groups));

        }
        contacts.make_matrix();
    }
}









            // for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
            //     for (size_t state_index = 0; state_index < 24; state_index++) {
            //         size_t idx = age_group.get() * 24 + state_index;

            //         std::cout<<model.populations[{age_group, InfectionState(state_index)}]<<std::endl;

            //         // model.populations[{age_group, InfectionState(state_index)}] = final_state[idx];
            //     }
            // }


            // node.property.get_simulation().set_integrator(make_integrator(integratorType(), dt());)
        // }



        // for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
        //     for (size_t state_index = 0; state_index < num_infection_states; state_index++) {
        //         size_t idx = age_group.get() * num_infection_states + state_index;
        //         model.populations[{age_group, InfectionState(state_index)}] = final_state[idx];
        //     }
        // }




        // mio::TimeSeries<FP> result = mio::simulate<FP, mio::osecirvvs::Model<FP>>(
        //     timeSteps[interval], timeSteps[interval+1], dt, model, integrator
        // );
        // const auto& final_state = result.get_last_value();

        // for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
        //     for (size_t state_index = 0; state_index < num_infection_states; state_index++) {
        //         size_t idx = age_group.get() * num_infection_states + state_index;
        //         model.populations[{age_group, InfectionState(state_index)}] = final_state[idx];
        //     }
        // }




    // // ---------------- //
    // // Start simulation //
    // // ---------------- //
    // for (size_t interval = 0; interval < numIntervals(); interval++) {

    //     const int controlIndex = interval / pcResolution();

    //     auto most_restrictive_control_value = [&](const std::string& name) {
    //         return m_control_policy[static_cast<size_t>(string_to_control(name))].range.second;
    //     };

    //     double school_closure             = most_restrictive_control_value("SchoolClosure");
    //     double home_office                = most_restrictive_control_value("HomeOffice");
    //     double physical_distancing_school = most_restrictive_control_value("PhysicalDistancingSchool");
    //     double physical_distancing_work   = most_restrictive_control_value("PhysicalDistancingWork");
    //     double physical_distancing_other  = most_restrictive_control_value("PhysicalDistancingOther");






    //     mio::TimeSeries<FP> result = mio::simulate<FP, mio::osecirvvs::Model<FP>>(
    //         timeSteps[interval], timeSteps[interval+1], dt, model, integrator
    //     );
    //     const auto& final_state = result.get_last_value();

    //     for (mio::AgeGroup age_group = 0; age_group < model.parameters.get_num_groups(); age_group++) {
    //         for (size_t state_index = 0; state_index < num_infection_states; state_index++) {
    //             size_t idx = age_group.get() * num_infection_states + state_index;
    //             model.populations[{age_group, InfectionState(state_index)}] = final_state[idx];
    //         }
    //     }








    // auto graph_sim_mobility = mio::make_mobility_sim<double>(t0(), dt(), graphSimulation());

//     // graph_sim_mobility.advance(tmax());
// }