#include <epidemiology/migration/migration.h>
#include <epidemiology/secir/seir.h>
#include <epidemiology/model/simulation.h>

#include <iostream>

int main(int argc, char** argv)
{
    const auto t0   = 0.;
    const auto tmax = 10.;
    const auto dt   = 0.5; //time step of migration, daily migration every second step

    epi::SeirModel model;
    model.populations[{epi::Index<epi::SeirInfType>(epi::SeirInfType::S)}] = 10000;
    model.parameters.set<epi::StageTimeIncubationInv>(1);
    model.parameters.get<epi::ContactFrequency>().get_baseline()(0, 0) = 2.7;
    model.parameters.set<epi::StageTimeInfectiousInv>(1);

    //two mostly identical groups
    auto model_group1 = model;
    auto model_group2 = model;
    //some contact restrictions in group 1
    model_group1.parameters.get<epi::ContactFrequency>().add_damping(0.5, epi::SimulationTime(5));
    //infection starts in group 1
    model_group1.populations[{epi::Index<epi::SeirInfType>(epi::SeirInfType::S)}] = 9990;
    model_group1.populations[{epi::Index<epi::SeirInfType>(epi::SeirInfType::E)}] = 10;

    epi::Graph<epi::SimulationNode<epi::Simulation<epi::SeirModel>>, epi::MigrationEdge> g;
    g.add_node(1001, model_group1, t0);
    g.add_node(1002, model_group2, t0);
    g.add_edge(0, 1, Eigen::VectorXd::Constant((size_t)epi::SeirInfType::Count, 0.01));
    g.add_edge(1, 0, Eigen::VectorXd::Constant((size_t)epi::SeirInfType::Count, 0.01));

    auto sim = epi::make_migration_sim(t0, dt, std::move(g));

    sim.advance(tmax);

    return 0;
}
