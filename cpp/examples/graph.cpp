#include <epidemiology/migration/migration.h>
#include <epidemiology/secir/seir.h>
#include <epidemiology/model/simulation.h>

#include <iostream>

int main(int argc, char** argv)
{
    const auto t0   = 0.;
    const auto tmax = 10.;
    const auto dt   = 1.; //time step of migration, not integration

    epi::SeirModel model = epi::create_seir_model();
    model.populations.set(10000, epi::InfectionType::S);
    model.parameters.times.set_incubation(1);
    model.parameters.times.set_cont_freq(2.7);
    model.parameters.times.set_infectious(1);

    //two mostly identical groups
    auto model_group1 = model;
    auto model_group2 = model;
    //some contact restrictions in group 1
    model_group1.parameters.dampings.add({5, 0.5});
    //infection starts in group 1
    model_group1.populations.set(9990, epi::InfectionType::S);
    model_group1.populations.set(10, epi::InfectionType::E);

    epi::Graph<epi::ModelNode<epi::Simulation<epi::SeirModel>>, epi::MigrationEdge> g;
    g.add_node(model_group1, t0);
    g.add_node(model_group2, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant((size_t)epi::InfectionType::Count, 0.01));
    g.add_edge(1, 0, Eigen::VectorXd::Constant((size_t)epi::InfectionType::Count, 0.01));

    auto sim = epi::make_migration_sim(t0, dt, g);

    sim.advance(10);

    return 0;
}
