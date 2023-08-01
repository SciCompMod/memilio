/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Julia Bicker, Daniel Abele
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef METAPOPULATION_MOBILITY_HYBRID_H
#define METAPOPULATION_MOBILITY_HYBRID_H

#include "memilio/mobility/metapopulation_mobility_instant.h"

namespace mio
{
template <class Agent, class AbmToOdeConversion, class OdeToAbmConversion, class OdeToAbmMapping>
class MigrationEdgeHybrid
{
public:
    //used if start node is abm node, because then we do not need MigrationParameters
    MigrationEdgeHybrid(const AbmToOdeConversion& abm_to_ode_fct, const OdeToAbmConversion& ode_to_abm_fct,
                        const OdeToAbmMapping& ode_to_abm_mapping, size_t num_compartments)
        : m_abm_to_ode_fct(abm_to_ode_fct)
        , m_ode_to_abm_fct(ode_to_abm_fct)
        , m_ode_to_abm_mapping(ode_to_abm_mapping)
        , m_migrated_compartments(num_compartments)
        , m_start_node_is_abm(true)
    {
    }

    //constructors with parameters are used if start node is ode node
    MigrationEdgeHybrid(const MigrationParameters& params, const AbmToOdeConversion& abm_to_ode_fct,
                        const OdeToAbmConversion& ode_to_abm_fct, const OdeToAbmMapping& ode_to_abm_mapping)
        : m_parameters(params)
        , m_migrated_compartments(params.get_coefficients().get_shape().rows())
        , m_abm_to_ode_fct(abm_to_ode_fct)
        , m_ode_to_abm_fct(ode_to_abm_fct)
        , m_ode_to_abm_mapping(ode_to_abm_mapping)
        , m_start_node_is_abm(false)
    {
    }

    MigrationEdgeHybrid(const Eigen::VectorXd& coeffs, const AbmToOdeConversion& abm_to_ode_fct,
                        const OdeToAbmConversion& ode_to_abm_fct, const OdeToAbmMapping& ode_to_abm_mapping)
        : m_parameters(coeffs)
        , m_migrated_compartments(coeffs.rows())
        , m_abm_to_ode_fct(abm_to_ode_fct)
        , m_ode_to_abm_fct(ode_to_abm_fct)
        , m_ode_to_abm_mapping(ode_to_abm_mapping)
        , m_start_node_is_abm(false)
    {
    }

    /**
     * @brief compute migration between an abm node and an ode node.
     * For the abm node, migration is based on agents' migration rules and for an ode node,
     * migration is based on coefficients. Migrations are deleted from the start node and added to the target node.
     * @param[in] t time of migration
     * @param[in, out] node_abm abm node. Can be start or target node.
     * @param[in, out] node_ode ode node. Can be start or target node.
    */
    template <class NodeAbm, class NodeOde>
    void apply_migration(double t, double dt, NodeAbm& node_abm, NodeOde& node_ode, NodeAbm& abm_node_to_ode_node)
    {
        if (m_start_node_is_abm) {
            //check whether agents should return
            if (m_migrated_agents.size() == 0) {
                size_t person_iter = 0;
                auto& persons_to_migrate =
                    node_abm.property.get_simulation().get_graph_world().get_persons_to_migrate();
                while (person_iter < persons_to_migrate.size()) {
                    //get the persons that want to go to the target node
                    if ((persons_to_migrate[person_iter]->get_location()).get_world_id() == node_ode.id) {
                        //agents a removed from start node world and added to migrated agents
                        m_migrated_agents.push_back(std::move(persons_to_migrate[person_iter]));
                        node_abm.get_simulation().get_graph_world().get_persons_to_migrate().erase(
                            node_abm.get_simulation().get_graph_world().get_persons_to_migrate().begin() + person_iter);
                        continue;
                    }
                    ++person_iter;
                }
                //convert migrating agents to ode compartments
                m_migrated_compartments.add_time_point(t, m_agents_to_compartments_fct(m_migrated_agents));
                //add migrating compartments to ode node
                node_ode.get_result().get_last_value() += m_migrated_compartments.get_last_value();
            }
            //agents return from ode node
            else {
                //returning compartments with euler step
                auto total =
                    find_value_reverse(node_ode.get_result(), m_migrated_compartments.get_last_time(), 1e-10, 1e-10);
                calculate_migration_returns(m_migrated_compartments.get_last_value(), node_ode.get_simulation(), *total,
                                            m_migrated_compartments.get_last_time(), dt);
                //the lower-order return calculation may in rare cases produce negative compartments,
                //especially at the beginning of the simulation.
                //fix by subtracting the supernumerous returns from the biggest compartment of the age group.
                Eigen::VectorXd remaining_after_return =
                    (node_ode.get_result().get_last_value() - m_migrated_compartments.get_last_value()).eval();
                for (Eigen::Index j = 0; j < node_ode.get_result().get_last_value().size(); ++j) {
                    if (remaining_after_return(j) < 0) {
                        auto num_comparts = (Eigen::Index)NodeOde::Sim::Model::Compartments::Count;
                        auto group        = Eigen::Index(j / num_comparts);
                        auto compart      = j % num_comparts;
                        log(remaining_after_return(j) < -1e-3 ? LogLevel::warn : LogLevel::info,
                            "Underflow during migration returns at time {}, compartment {}, age group {}: {}", t,
                            compart, group, remaining_after_return(j));
                        Eigen::Index max_index;
                        slice(remaining_after_return, {group * num_comparts, num_comparts}).maxCoeff(&max_index);
                        log_info("Transferring to compartment {}", max_index);
                        max_index += group * num_comparts;
                        m_migrated_compartments.get_last_value()(max_index) -= remaining_after_return(j);
                        m_migrated_compartments.get_last_value()(j) += remaining_after_return(j);
                    }
                }
                //substract returning compartments from ode nodes
                node_ode.get_result().get_last_value() -= m_migrated_compartments.get_last_value();
                //convert returning compartments to agents by setting migrated agents' status and location
                m_ode_to_abm_fct(m_migrated_compartments.get_last_value(), m_migrated_agents);
                //add returning agents to abm node and remove them from migrated agents
                while (m_migrated_agents.size() > 0) {
                    node_abm.get_simulation().get_graph_world().add_existing_person(std::move(m_migrated_agents[0]));
                    m_migrated_agents.erase(m_migrated_agents.begin());
                }
                //remove compartments that migrated to abm
                m_migrated_compartments.remove_time_point(m_migrated_compartments.get_last_time());
            }
        }
        //start node is ode
        else {
            //check whether agents should return
            if (m_migrated_compartments.get_num_time_points() == 0) {

                if ((m_parameters.get_coefficients().get_matrix_at(t).array() > 0.0).any()) {
                    //get migrating compartments
                    m_migrated_compartments.add_time_point(
                        t,
                        (node_ode.get_last_state().array() * m_parameters.get_coefficients().get_matrix_at(t).array() *
                         get_migration_factors(node_ode, t, node_ode.get_last_state()).array())
                            .matrix());
                    //get agents corresponding to migrating compartments
                    m_ode_to_abm_mapping(m_migrated_compartments.get_last_value(),
                                         abm_node_to_ode_node.property.get_simulation().get_graph_world().get_persons(),
                                         m_migrated_agents);
                    //add persons to abm node they migrate to
                    while (m_migrated_agents.size() > 0) {
                        auto& target_location = node_abm.get_simulation().get_graph_world().find_location(
                            (m_migrated_agents[0]->get_location()).get_type(), *m_migrated_agents[0]);
                        //let person migrate to location in target world
                        m_migrated_agents[0]->migrate_to_other_world(target_location, true);
                        //add person to abm node
                        node_abm.get_simulation().get_graph_world().add_existing_person(
                            std::move(m_migrated_agents[0]));
                        m_migrated_agents.erase(m_migrated_agents.begin());
                    }
                }
            }
            //agents should return
            else {
                //if agents return from abm node, they are in m_persons_to_migrate in the abm node they come from
                size_t person_iter = 0;
                auto& persons_to_migrate =
                    node_abm.property.get_simulation().get_graph_world().get_persons_to_migrate();
                while (person_iter < persons_to_migrate.size()) {
                    //get the persons that want to go to the target node
                    if ((persons_to_migrate[person_iter]->get_location()).get_world_id() == node_ode.id) {
                        //agents a removed from start node world and added to migrated agents
                        m_migrated_agents.push_back(std::move(persons_to_migrate[person_iter]));
                        node_abm.get_simulation().get_graph_world().get_persons_to_migrate().erase(
                            node_abm.get_simulation().get_graph_world().get_persons_to_migrate().begin() + person_iter);
                        continue;
                    }
                    ++person_iter;
                }
                //remove returned agents from migrated compartments
                m_migrated_compartments.remove_time_point(m_migrated_compartments.get_last_time());
                //get compartments that return from abm
                m_migrated_compartments.add_time_point(t, m_agents_to_compartments_fct(m_migrated_agents));
                //add compartments to ode node
                node_ode.get_result().get_last_value() += m_migrated_compartments.get_last_value();
                //remove compartments from m_migrated_compartments again after that have been added to ode node
                m_migrated_compartments.remove_time_point(m_migrated_compartments.get_last_time());
                //add agents to abm corresponding to ode node and remove them from m_migrated_agents
                while (m_migrated_agents.size() > 0) {
                    auto& target_location = abm_node_to_ode_node.get_simulation().get_graph_world().find_location(
                        (m_migrated_agents[0]->get_location()).get_type(), *m_migrated_agents[0]);
                    m_migrated_agents[0]->migrate_to_other_world(target_location, true);
                    abm_node_to_ode_node.get_simulation().get_graph_world().add_existing_person(
                        std::move(m_migrated_agents[0]));
                    m_migrated_agents.erase(m_migrated_agents.begin());
                }
            }
        }
    }

private:
    boost::optional<MigrationParameters>
        m_parameters; ///< Parameters to calculated migrating compartments is start node is ode node
    boost::optional<std::vector<std::unique_ptr<Agent>>>
        m_migrated_agents; ///< Agents that have migrated to target node
    TimeSeries<double> m_migrated_compartments; ///< compartments that have migrated to target node
    bool m_start_node_is_abm; ///< whether node from is an abm node
    AbmToOdeConversion
        m_agents_to_compartments_fct; ///< Gets a vector with agents and returns the corresponding compartments
    OdeToAbmConversion
        m_ode_to_abm_fct; ///< Gets an ode node and the corresponding abm node. Updates the agents from the abm node (status and location).
    OdeToAbmMapping
        m_ode_to_abm_mapping; ///< Gets the migrating compartments and the abm node from which they should migrate. Picks agents according to migrating compartments from abm node,
    ///< deletes them from abm node, sets their target location and adds them to m_migrated_agents
};

/**
 * edge functor for hybrid migration simulation.
 * @see MigrationEdgeHybrid::apply_migration
 */
template <class NodeOde, class NodeAbm, class... EdgeParams>
void apply_migration(double t, double dt, MigrationEdgeHybrid<EdgeParams...>& migrationEdge, NodeOde& node_from,
                     NodeAbm& node_to, NodeAbm& abm_to_ode_node)
{
    migrationEdge.apply_migration(t, dt, node_from, node_to, abm_to_ode_node);
}

} // namespace mio

#endif //METAPOPULATION_MOBILITY_HYBRID_H
