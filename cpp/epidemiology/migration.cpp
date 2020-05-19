#include "epidemiology/migration.h"
#include "epidemiology/adapt_rk.h"
#include "epidemiology/eigen_util.h"

namespace epi
{

std::vector<double> ode_integrate_with_migration(double t0, double tmax, double dt,
                                                 const std::vector<std::unique_ptr<IntegratorBase>>& integrators,
                                                 MigrationFunction migration_function,
                                                 std::vector<Eigen::VectorXd>& result)
{
    assert(!integrators.empty());
    assert(result.size() == 1);

    auto num_groups     = integrators.size();
    auto num_vars_total = result[0].size();
    auto num_vars       = num_vars_total / num_groups;
    auto num_steps      = static_cast<size_t>(std::ceil(tmax - t0));
    auto t              = std::vector<double>(1, t0);
    t.reserve(num_steps);
    result.reserve(num_steps);

    auto migration     = Eigen::MatrixXd(num_groups, num_groups); //reusable, matrix of migration between groups
    auto result_single_group = std::vector<Eigen::VectorXd>(1, Eigen::VectorXd(num_vars)); //reusable, result of a single group and migration step

    while (t.back() < tmax) {
        //TODO: avoid copying the results of a single step/group
        //TODO: variable migration step size
        //TODO: report all steps, not just once per migration interval

        //simulate the same timestep for every group
        auto t0_step   = t.back();
        auto tmax_step = std::min(t0_step + 1, tmax);
        result.emplace_back(num_vars_total);
        const auto& init_step = result[result.size() - 2];
        auto& result_step     = result[result.size() - 1];
        for (size_t group_idx = 0; group_idx < num_groups; group_idx++) {
            result_single_group.resize(
                1, slice(init_step, {(Eigen::Index)group_idx, (Eigen::Index)num_vars, (Eigen::Index)num_groups}));
            ode_integrate(t0_step, tmax_step, dt, *integrators[group_idx], result_single_group);
            slice(result_step, {(Eigen::Index)group_idx, (Eigen::Index)num_vars, (Eigen::Index)num_groups}) =
                result_single_group.back();
        }

        //migration for each variable
        for (size_t var_idx = 0; var_idx < num_vars; var_idx++) {
            migration_function(var_idx, t.back(), migration);
            auto one_var_all_groups = result_step.segment(num_groups * var_idx, num_groups);
            one_var_all_groups      = migration * one_var_all_groups;
        }

        t.push_back(tmax_step);
    }

    return t;
}

} // namespace epi