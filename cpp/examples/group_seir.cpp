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
    auto migration_function = [](size_t var_idx, double t, Eigen::MatrixXd& migration) {
        if (var_idx == 2) {
            migration = Eigen::MatrixXd::Identity(2, 2); //no migration of infectious people
        }
        else {
            migration = Eigen::MatrixXd::Constant(2, 2, 0.5); //complete mixing of everyone else
        }
    };
    auto t_groups = epi::simulate_groups(0, tmax, dt, {params_group1, params_group2}, migration_function, result_groups);

    return 0;
}