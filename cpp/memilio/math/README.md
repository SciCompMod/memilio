# Math classes

This module contains different mathematical methods for integration of sets of ordinary differential equations. For the moment, it contains explicit methods for generic right hand sides and a work-in-progress (semi-)implicit method for the particular SECIR model.

## Structure

The model consists of the following classes:
1. Integrator: The integrator module contains an IntegratorCore and an OdeIntegrator class, currently designed for explicit methods only. The IntegratorCore contains a step function that takes as input the right hand side function of the ODE, the current time, the step size etc. It represents a generic integration method from which the other integrators are derived from. The OdeIntegrator stores an IntegratorCore as well as the right hand side function and a time series of results.
2. The following integrators implement IntegratorCore and can be used with OdeIntegratorCore or [mio::Simulation](../compartments/README.md).
    - ControlledStepperWrapper: The ControlledStepperWrapper class allows using integrators from boosts::numeric::odeint in memilio. The integrator is passed as template argument. It must be of the concept "Error Stepper" (for more details, see [boost's own documentation](https://www.boost.org/doc/libs/1_75_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/odeint_in_detail/steppers.html#boost_numeric_odeint.odeint_in_detail.steppers.stepper_overview) ). Currently, the default for mio::Simulation is runge_kutta_cash_karp54. Alternatively, runge_kutta_dopri5 or runge_kutta_fehlberg78 can be used.
    - Euler: The Euler class contains and explicit Euler method, adapted to the Integrator function. It also contains a (semi-)implicit Euler method which is WIP and taylored for the particular SECIR model from ../ode_secir/model.h
    - Adapt_RK: The Adapt_RK module contains the RKIntegratorCore class, which implements adaptive Runge-Kutta integrators, where different pairs of methods (in form of combined Butcher tableaus) can be added. Absolute and relative tolerances can be set and the Tableau in use is that of an adaptive Runge-Kutta-Fehlberg (45) method; see, e.g., https://www.johndcook.com/blog/2020/02/19/fehlberg/. Steps where the mixed criterion on absolute and relative values (m_abs_tol + max_val * m_rel_tol) are not satisfied are directly discarded and never used. If the minimal step size (set) is reached and the criterion cannot be satisfied, it is returned that the adaptive step sizing failed.
3. Smoother: The smoother classes smoothes discrete jumps of function values y0 and y1 on the interval [x0,x1] by a continuously differentiable function

## Example

Example can be found at 
- examples/euler_test.cpp
- examples/adapt_rk_test.cpp
