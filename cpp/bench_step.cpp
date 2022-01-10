#include "bench.h"

#define BENCH_STEP_MODEL make_model(1)
#define BENCH_STEP_MODEL_TYPE mio::SecirModel

void BM_STEP_reset_problem(double& t, double& dt) {
    t       = 50;
    dt      = 1; // adapt_rk says 3.2, vadapt_rk_opt says ~2.12
}

void BM_STEP_set_problem(
    mio::SecirModel& model,
    double& abs_tol, double& rel_tol, double& dt_min, double& dt_max,
    double& t, double& dt,
    mio::DerivFunction& f, Eigen::VectorXd& yt, Eigen::VectorXd& ytp1
) {
    BM_STEP_reset_problem(t, dt);
    abs_tol = 1e-10;
    rel_tol = 1e-5;
    dt_min  = std::numeric_limits<double>::min();
    dt_max  = std::numeric_limits<double>::max();
    f       = [model](Eigen::Ref<const Eigen::VectorXd> y_f, double t_f, Eigen::Ref<Eigen::VectorXd> dydt_f) { model.eval_right_hand_side(y_f, y_f, t_f, dydt_f); };
    // [](Eigen::Ref<const Eigen::VectorXd> y_f, double t_f, Eigen::Ref<Eigen::VectorXd> dydt_f)){model.get_derivatives();}
    // get_derivatives(Eigen::Ref<const Eigen::VectorXd> pop,
    //                  Eigen::Ref<const Eigen::VectorXd> y, double t,
    //                  Eigen::Ref<Eigen::VectorXd> dydt)
    //std::function<void(Eigen::Ref<const Eigen::VectorXd> y, double t, Eigen::Ref<Eigen::VectorXd> dydt)>;
    yt      = Eigen::Matrix<double, 8, 1>::Zero();
    double yt_vals[8] = {6377.873644, 35.249156, 30.029611, 182.145865, 66.153059, 79.530621, 3069.383604, 159.634440};
    for (int i = 0; i < 8; i++) {
        yt[i] = yt_vals[i];
    }
    ytp1    = Eigen::Matrix<double, 8, 1>::Zero();
}

// const DerivFunction& f, Eigen::Ref<Eigen::VectorXd const> yt, double& t, double& dt,
//              Eigen::Ref<Eigen::VectorXd> ytp1
    
template <class Integrator>
static void BM_STEP_adapt_rk(benchmark::State& state) {
    mio::set_log_level(mio::LogLevel::off);
    auto model = BENCH_STEP_MODEL;
    double abs_tol, rel_tol, dt_min, dt_max, t, dt;
    mio::DerivFunction f;
    Eigen::VectorXd yt, ytp1;


    BM_STEP_set_problem(model, abs_tol, rel_tol, dt_min, dt_max, t, dt, f, yt, ytp1);
    
    auto I = Integrator(abs_tol, rel_tol, dt_min, dt_max);

    for (auto _ : state) {
        // This code gets timed
        I.step(f, yt, t, dt, ytp1);
        //BM_STEP_reset_problem(t, dt);
    }
}

BENCHMARK_TEMPLATE(BM_STEP_adapt_rk, mio::RKIntegratorCore)
    ->Name("Dummy 1/3");
BENCHMARK_TEMPLATE(BM_STEP_adapt_rk, mio::RKIntegratorCore)
    ->Name("Dummy 2/3");
BENCHMARK_TEMPLATE(BM_STEP_adapt_rk, mio::RKIntegratorCore)
    ->Name("Dummy 3/3");
BENCHMARK_TEMPLATE(BM_STEP_adapt_rk, mio::RKIntegratorCore)
    ->Name("simulate SecirModel adapt_rk");
BENCHMARK_TEMPLATE(BM_STEP_adapt_rk, mio::RKIntegratorCore2)
    ->Name("simulate SecirModel vadapt_rk");
BENCHMARK_TEMPLATE(BM_STEP_adapt_rk, mio::RKIntegratorCoreAll)
    ->Name("simulate SecirModel adapt_rk_all");
BENCHMARK_TEMPLATE(BM_STEP_adapt_rk, mio::RKIntegratorCoreAll2)
    ->Name("simulate SecirModel vadapt_rk_all");
BENCHMARK_TEMPLATE(BM_STEP_adapt_rk, mio::RKIntegratorCore3)
    ->Name("simulate SecirModel vadapt_rk_opt");
BENCHMARK_TEMPLATE(BM_STEP_adapt_rk, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_cash_karp54>)
    ->Name("simulate SecirModel boost rk_ck54");
BENCHMARK_TEMPLATE(BM_STEP_adapt_rk, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_dopri5>)
    ->Name("simulate SecirModel boost rk_dopri5");
BENCHMARK_TEMPLATE(BM_STEP_adapt_rk, mio::ControlledStepperWrapper<boost::numeric::odeint::runge_kutta_fehlberg78>)
    ->Name("simulate SecirModel boost rkf78");
BENCHMARK_MAIN();