#ifndef MIGRATION_H
#define MIGRATION_H

#include "epidemiology/migration/graph_simulation.h"
#include "epidemiology/utils/time_series.h"
#include "epidemiology/utils/eigen.h"
#include <vector>

namespace epi
{

{


/**
 * represents the simulation in one node of the graph.
 */
template <class Model>
class ModelNode
{
public:
    template <class... Args, typename = std::enable_if_t<std::is_constructible<Model, Args...>::value, void>>
    ModelNode(Args&&... args)
        : model(std::forward<Args>(args)...)
    {
    }

    /**
     * get the result of the simulation in this node.
     */
    decltype(auto) get_result() const
    {
        return model.get_result();
    }
    
    /**
     * get the parameters of the simulation in this node.
     */
    decltype(auto) get_params() const
    {
        return model.get_params();
    }

    Model model;
    Eigen::VectorXd last_state;
};

class MigrationEdge
{
public:
    MigrationEdge(const Eigen::VectorXd& coeffs)
        : coefficients(coeffs)
    {
    }
    //per group and compartment
    Eigen::VectorXd coefficients;

    bool operator==(const MigrationEdge& other) const 
    {
        return coefficients == other.coefficients;
    }
};

template <class Model>
void evolve_model(double t, double dt, ModelNode<Model>& model)
{
    model.model.advance(t + dt);
    model.last_state = model.model.get_result().get_last_value();
}

template <class Model>
void apply_migration(double /*t*/, double dt, MigrationEdge& migrationEdge, ModelNode<Model>& model1,
                     ModelNode<Model>& model2)
{
    auto migration = (dt * model1.last_state.array() * migrationEdge.coefficients.array()).matrix();
    model2.model.get_result().get_last_value() += migration;
    model1.model.get_result().get_last_value() -= migration;
}

template <typename Model>
GraphSimulation<Graph<ModelNode<Model>, MigrationEdge>>
make_migration_sim(double t0, double dt, const Graph<ModelNode<Model>, MigrationEdge>& graph)
{
    return make_graph_sim(t0, dt, graph, &evolve_model<Model>, &apply_migration<Model>);
}

} // namespace epi

#endif //MIGRATION_H
