#include <epidemiology/migration/migration.h>
#include <epidemiology/secir/seir.h>

#include <iostream>

int main(int argc, char** argv)
{
    const auto t0   = 0.;
    const auto tmax = 10.;
    const auto dt   = 1.; //time step of migration, not integration

    epi::SeirParams params;
    params.populations.set({epi::SeirCompartments::S}, 10000);
    params.times.set_incubation(1);
    params.times.set_cont_freq(2.7);
    params.times.set_infectious(1);

    //two mostly identical groups
    auto params_group1 = params;
    auto params_group2 = params;
    //some contact restrictions in group 1
    params_group1.dampings.add({5, 0.5});
    //infection starts in group 1
    params_group1.populations.set({epi::SeirCompartments::S}, 9990);
    params_group1.populations.set({epi::SeirCompartments::E}, 10);

    epi::Graph<epi::ModelNode<epi::SeirSimulation>, epi::MigrationEdge> g;
    g.add_node(params_group1, t0);
    g.add_node(params_group2, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant(epi::SeirCompartments::SeirCount, 0.01));
    g.add_edge(1, 0, Eigen::VectorXd::Constant(epi::SeirCompartments::SeirCount, 0.01));

    auto sim = epi::make_migration_sim(t0, dt, g);

    sim.advance(10);

    return 0;
}
