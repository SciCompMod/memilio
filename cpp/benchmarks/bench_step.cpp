#include "bench_core.h"
#include "bench_model.h"

// google's benchmark makros require a function in the global scope, hence the following "using" statement
using mio::benchmark::integrator_step;

// dummy runs to avoid large effects of cpu scaling on times of actual benchmarks
BENCHMARK_TEMPLATE(integrator_step, mio::RKIntegratorCore)
    ->Name("Dummy 1/3");
BENCHMARK_TEMPLATE(integrator_step, mio::RKIntegratorCore)
    ->Name("Dummy 2/3");
BENCHMARK_TEMPLATE(integrator_step, mio::RKIntegratorCore)
    ->Name("Dummy 3/3");
// register functions as a benchmarks and set a name
BENCHMARK_TEMPLATE(integrator_step, mio::RKIntegratorCore)
    ->Name("simulate SecirModel adapt_rk");
BENCHMARK_TEMPLATE(integrator_step, mio::VRKIntegratorCore)
    ->Name("simulate SecirModel vadapt_rk");
BENCHMARK_TEMPLATE(integrator_step, mio::RKAllIntegratorCore)
    ->Name("simulate SecirModel adapt_rk_all");
BENCHMARK_TEMPLATE(integrator_step, mio::VRKAllIntegratorCore)
    ->Name("simulate SecirModel vadapt_rk_all");
BENCHMARK_TEMPLATE(integrator_step, mio::VRKOptIntegratorCore)
    ->Name("simulate SecirModel vadapt_rk_opt");
BENCHMARK_TEMPLATE(integrator_step, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>)
    ->Name("simulate SecirModel boost rk_ck54");
BENCHMARK_TEMPLATE(integrator_step, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_dopri5>)
    ->Name("simulate SecirModel boost rk_dopri5");
BENCHMARK_TEMPLATE(integrator_step, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_fehlberg78>)
    ->Name("simulate SecirModel boost rkf78");
// run all benchmarks
BENCHMARK_MAIN();