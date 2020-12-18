#include "epidemiology/migration/migration.h"
#include "epidemiology/secir/secir.h"
#include "epidemiology/secir/seir.h"

namespace epi
{

void MigrationEdge::calculate_returns_ode(Eigen::Ref<TimeSeries<double>::Vector> migrated, const SecirParams& params,
                                          Eigen::Ref<const TimeSeries<double>::Vector> total, double t, double dt)
{
    auto y0 = migrated.eval();
    auto y1 = migrated;
    EulerIntegratorCore().step(
        [&](auto&& y, auto&& t_, auto&& dydt) {
            secir_get_derivatives(params, total, y, t_, dydt);
        },
        y0, t, dt, y1);
}

void MigrationEdge::calculate_returns_ode(Eigen::Ref<TimeSeries<double>::Vector> migrated, const SeirParams& params,
                                          Eigen::Ref<const TimeSeries<double>::Vector> total, double t, double dt)
{
    auto y0 = migrated.eval();
    auto y1 = migrated;
    EulerIntegratorCore().step(
        [&](auto&& y, auto&& t_, auto&& dydt) {
            seir_get_derivatives(params, total, y, t_, dydt);
        },
        y0, t, dt, y1);
}

} // namespace epi