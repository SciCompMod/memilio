#ifndef MIGRATION_H
#define MIGRATION_H

#include <epidemiology/graph.h>
#include <Eigen/Core>
#include <vector>

namespace epi
{

template <class Params>
class CompartmentModel
{
public:
    CompartmentModel(Params p)
        : params(p)
        , time_series(1) //TODO: Init
    {
    }

    Params params;
    std::vector<Eigen::VectorXd> time_series;
    Eigen::VectorXd last_evolve_state;
};

class Migration
{
public:
    Migration(double coefficient)
        : coefficient(coefficient)
    {
    }
    double coefficient;
};

template <class Params>
void evolve_model(double t, double dt, CompartmentModel<Params>& model)
{
    //TODO: this doesn't quite work yet because simulate can't be continued
    std::vector<Eigen::VectorXd> step(1, model.time_series.back());
    simulate(t, t + dt, 0.1 * dt, model.params, step);
    model.time_series.push_back(step.back());
    model.last_evolve_state = model.time_series.back(); //store for reading during migration
}

template <class Params>
void apply_migration(double t, double dt, Migration& migration, CompartmentModel<Params>& model1,
                     CompartmentModel<Params>& model2)
{
    //TODO: this model is almost certainly too coarse, all compartments don't migrate at the same rate (i.e. same coefficient)
    auto dim           = model1.time_series.back().size();
    auto abs_migration = dt * migration.coefficient * model1.last_evolve_state;
    model2.time_series.back() += abs_migration;
    model1.time_series.back() -= abs_migration;
}

template <typename Params>
void simulate_migration(double t0, double tmax, double dt, Graph<CompartmentModel<Params>, Migration>& graph)
{
    simulate_graph(t0, tmax, dt, graph, &evolve_model<Params>, &apply_migration<Params>);
}

} // namespace epi

#endif //MIGRATION_H