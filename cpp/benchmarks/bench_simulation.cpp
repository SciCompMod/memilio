#include "bench_core.h"
#include "bench_model.h"

// google's benchmark templates require a function in the global scope, hence the following "using" statement
using mio::benchmark::simulation;
// use a global string since string literals are not allowed as template arguments
const char configpath[] = "benchmarks/simulation_benchmark.config";
using init = mio::benchmark::simulation_file_init<configpath, mio::benchmark::model::SecirAgeresAbsurdDampings>;

// dummy runs to avoid large effects of cpu scaling on times of actual benchmarks
BENCHMARK_TEMPLATE(simulation, mio::RKIntegratorCore, init)
    ->Name("Dummy 1/3");
BENCHMARK_TEMPLATE(simulation, mio::RKIntegratorCore, init)
    ->Name("Dummy 2/3");
BENCHMARK_TEMPLATE(simulation, mio::RKIntegratorCore, init)
    ->Name("Dummy 3/3");
// register functions as a benchmarks and set a name
BENCHMARK_TEMPLATE(simulation, mio::RKIntegratorCore, init)
    ->Name("simulate SecirModel adapt_rk");
BENCHMARK_TEMPLATE(simulation, mio::VRKIntegratorCore, init)
    ->Name("simulate SecirModel vadapt_rk");
BENCHMARK_TEMPLATE(simulation, mio::RKAllIntegratorCore, init)
    ->Name("simulate SecirModel adapt_rk_all");
BENCHMARK_TEMPLATE(simulation, mio::VRKAllIntegratorCore, init)
    ->Name("simulate SecirModel vadapt_rk_all");
BENCHMARK_TEMPLATE(simulation, mio::VRKOptIntegratorCore, init)
    ->Name("simulate SecirModel vadapt_rk_opt");
BENCHMARK_TEMPLATE(simulation, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>, init)
    ->Name("simulate SecirModel boost rk_ck54");
BENCHMARK_TEMPLATE(simulation, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_dopri5>, init)
    ->Name("simulate SecirModel boost rk_dopri5");
BENCHMARK_TEMPLATE(simulation, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_fehlberg78>, init)
    ->Name("simulate SecirModel boost rkf78");
// run all benchmarks
BENCHMARK_MAIN();