#include <epidemiology/migration.h>
#include <epidemiology/seir.h>
#include <iostream>

int main(int argc, char** argv)
{
    const auto t0   = 0.;
    const auto tmax = 10.;
    const auto dt   = 1.; //time step of migration, not integration

    epi::SeirParams params;
    params.populations.set_exposed_t0(0);
    params.populations.set_infectious_t0(0);
    params.populations.set_total_t0(10000);
    params.populations.set_recovered_t0(0);
    params.times.set_incubation(1);
    params.times.set_cont_freq(2.7);
    params.times.set_infectious(1);

    //two mostly identical groups
    auto params_group1 = params;
    auto params_group2 = params;
    //some contact restrictions in group 1
    params_group1.dampings.add({5, 0.5});
    //infection starts in group 1
    params_group1.populations.set_exposed_t0(10);
    
    epi::Graph<epi::CompartmentModel<epi::SeirParams>, epi::Migration> g;
    g.add_node(params_group1);
    g.add_node(params_group2);
    g.add_edge(0, 1, 0.01);
    g.add_edge(1, 0, 0.01);

    simulate_migration(t0, tmax, dt, g);

    return 0;
}