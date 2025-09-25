#pragma once

#include "memilio/config.h"
#include "memilio/math/integrator.h"
#include "memilio/math/euler.h"
#include "memilio/math/adapt_rk.h"
#include "boost/numeric/odeint.hpp"


enum class IntegratorType
{
    ExplicitEuler,
    ExplicitCashKarp54,
    ExplicitFehlberg78,
    ControlledCashKarp54,
    ControlledFehlberg78
};

// template <typename FP>
// std::unique_ptr<mio::IntegratorCore<FP>> make_integrator(IntegratorType type, FP dt)
// {
//     FP dt_min  = std::numeric_limits<FP>::min();
//     FP dt_max  = dt;
//     FP abs_tol = 1e-12;
//     FP rel_tol = 1e-8;

//     using namespace boost::numeric::odeint;

//     switch (type) {
//     case IntegratorType::ExplicitEuler: {
//         auto eulerIntegrator = std::make_unique<mio::EulerIntegratorCore<FP>>();
//         return eulerIntegrator;
//     }
//     case IntegratorType::ExplicitCashKarp54: {
//         auto cashKarpStepper          = std::make_unique<mio::ExplicitStepperWrapper<FP, runge_kutta_cash_karp54>>();
//         cashKarpStepper->get_dt_min() = dt_min;
//         cashKarpStepper->get_dt_max() = dt_max;
//         return cashKarpStepper;
//     }
//     case IntegratorType::ExplicitFehlberg78: {
//         auto fehlbergStepper          = std::make_unique<mio::ExplicitStepperWrapper<FP, runge_kutta_fehlberg78>>();
//         fehlbergStepper->get_dt_min() = dt_min;
//         fehlbergStepper->get_dt_max() = dt_max;
//         return fehlbergStepper;
//     }
//     case IntegratorType::ControlledCashKarp54: {
//         auto cashKarpControlled = std::make_unique<mio::ControlledStepperWrapper<FP, runge_kutta_cash_karp54>>();
//         cashKarpControlled->set_dt_min(dt_min);
//         cashKarpControlled->set_dt_max(dt_max);
//         cashKarpControlled->set_abs_tolerance(abs_tol);
//         cashKarpControlled->set_rel_tolerance(rel_tol);
//         return cashKarpControlled;
//     }
//     case IntegratorType::ControlledFehlberg78: {
//         auto fehlbergControlled = std::make_unique<mio::ControlledStepperWrapper<FP, runge_kutta_fehlberg78>>();
//         fehlbergControlled->set_dt_min(dt_min);
//         fehlbergControlled->set_dt_max(dt_max);
//         fehlbergControlled->set_abs_tolerance(abs_tol);
//         fehlbergControlled->set_rel_tolerance(rel_tol);
//         return fehlbergControlled;
//     }
//     default: {
//         throw std::invalid_argument("Unknown IntegratorType.");
    
//     }
//     }
// }
