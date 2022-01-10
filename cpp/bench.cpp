#include "bench.h"

// eps_abs, eps_rel, dt_min, dt_max
#define BENCH_INTEGRATOR_ARGS 1e-10, 1e-5, std::numeric_limits<double>::min(), std::numeric_limits<double>::max()
// t0, tmax, dt
#define BENCH_SIMULATION_TIME 0, 100, 0.1

#define BENCH_MODEL make_model(10)


template <class Integrator>
void BM_Simulate(benchmark::State& state) {
   mio::set_log_level(mio::LogLevel::off);
    auto model = BENCH_MODEL;

    for (auto _ : state) {
        // This code gets timed
        std::shared_ptr<mio::IntegratorCore> I = std::make_shared<Integrator>(BENCH_INTEGRATOR_ARGS);
        simulate(BENCH_SIMULATION_TIME, model, I);
    }
}

// Register the function as a benchmark
BENCHMARK_TEMPLATE(BM_Simulate, mio::RKIntegratorCore)
    ->Name("Dummy 1/3");
BENCHMARK_TEMPLATE(BM_Simulate, mio::RKIntegratorCore)
    ->Name("Dummy 2/3");
BENCHMARK_TEMPLATE(BM_Simulate, mio::RKIntegratorCore)
    ->Name("Dummy 3/3");
BENCHMARK_TEMPLATE(BM_Simulate, mio::RKIntegratorCore)
    ->Name("simulate SecirModel adapt_rk");
BENCHMARK_TEMPLATE(BM_Simulate, mio::RKIntegratorCore2)
    ->Name("simulate SecirModel vadapt_rk");
BENCHMARK_TEMPLATE(BM_Simulate, mio::RKIntegratorCoreAll)
    ->Name("simulate SecirModel adapt_rk_all");
BENCHMARK_TEMPLATE(BM_Simulate, mio::RKIntegratorCoreAll2)
    ->Name("simulate SecirModel vadapt_rk_all");
BENCHMARK_TEMPLATE(BM_Simulate, mio::RKIntegratorCore3)
    ->Name("simulate SecirModel vadapt_rk_opt");
BENCHMARK_TEMPLATE(BM_Simulate, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>)
    ->Name("simulate SecirModel boost rk_ck54");
BENCHMARK_TEMPLATE(BM_Simulate, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_dopri5>)
    ->Name("simulate SecirModel boost rk_dopri5");
BENCHMARK_TEMPLATE(BM_Simulate, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_fehlberg78>)
    ->Name("simulate SecirModel boost rkf78");
// Run the benchmark
BENCHMARK_MAIN();