#ifndef MIGRATION_H
#define MIGRATION_H

#include <Eigen/Core>
#include <vector>

#include <epidemiology/graph.h>

namespace epi
{

template <class Model>
class ModelNode
{
public:
    template <class... Args, typename = std::enable_if_t<std::is_constructible<Model, Args...>::value, void>>
    ModelNode(Args&&... args)
        : model(std::forward<Args>(args)...)
    {
    }
    Model model;
    Eigen::VectorXd last_state;
};

class MigrationEdge
{
public:
    MigrationEdge(const Eigen::VectorXd& coefficients)
        : coefficients(coefficients)
    {
    }
    //per group and compartment
    Eigen::VectorXd coefficients;
};

template <class Model>
void evolve_model(double t, double dt, ModelNode<Model>& model)
{
    model.model.advance(t + dt);
    model.last_state = model.model.get_y().back();
}

template <class Model>
void apply_migration(double t, double dt, MigrationEdge& migrationEdge, ModelNode<Model>& model1,
                     ModelNode<Model>& model2)
{
    auto migration = (dt * model1.last_state.array() * migrationEdge.coefficients.array()).matrix();
    model2.model.get_y().back() += migration;
    model1.model.get_y().back() -= migration;
}

template <typename Model>
auto make_migration_sim(double t0, double dt, const Graph<ModelNode<Model>, MigrationEdge>& graph)
{
    return make_graph_sim(t0, dt, graph, &evolve_model<Model>, &apply_migration<Model>);
}

} // namespace epi

#endif //MIGRATION_H
