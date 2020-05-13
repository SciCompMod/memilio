#include "epidemiology/seir.h"

int main()
{
    auto t0   = 0.;
    auto tmax = 10.;
    auto dt   = 0.1002004008016032;

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
    //some contact restrictions in one group
    params_group1.dampings.add({5, 0.5});
    //infection starts in group 1
    params_group1.populations.set_exposed_t0(10);

    std::vector<Eigen::VectorXd> result_groups;
    auto migration_function = [](epi::SeirCompartment c, double t, Eigen::MatrixXd& migration) {
        if (c == epi::SeirCompartment::Infectious) {
            migration = Eigen::MatrixXd::Identity(2, 2); //no migration of sick people
        }
        else {
            migration = Eigen::MatrixXd::Constant(2, 2, 0.5); //complete mixing of everyone else
        }
    };
    auto t_groups = simulate(0, tmax, dt, std::vector<epi::SeirParams>{params_group1, params_group2},
                             migration_function, result_groups);

    return 0;
}