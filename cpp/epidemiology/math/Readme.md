# Math classes

This module contains different mathematical methods for integration of sets of ordinary differential equations. For the moment, it contains explicit methods for generic right hand sides and a work-in-progress (semi-)implicit method for the particular SECIR model.

## Structure

The model consists of the following classes:
1. Integrator: The integrator module contains an IntegratorCore and an OdeIntegrator class, currently designed for explicit methods only. The integrator core contains a step function that takes as input the right hand side function of the ODE, the current time, the step size etc. It represents a generic integration method from which EulerIntegratorCore oder RKIntegratorCore are derived from. The OdeIntegrator stores an IntegratorCore as well as the right hand side function and a time series of results.
2. Euler: The Euler class contains and explicit Euler method, adapted to the Integrator function. It also contains a (semi-)implicit Euler method which is WIP and taylored for the particular SECIR model from ../secir/secir.cpp
3. Adapt_RK: The Adapt_RK class implements adaptive Runge-Kutta integrators where different pairs of methods (in form of combined Butcher tableaus) can be added. Absolute and relative tolerances can be set and the Tableau in use is that of an adaptive Runge-Kutta-Fehlberg (45) method; see, e.g., https://www.johndcook.com/blog/2020/02/19/fehlberg/. Steps where the mixed criterion on absolute and relative values (m_abs_tol + max_val * m_rel_tol) are not satisfied are directly discarded and never used. If the minimal step size (set) is reached and the criterion cannot be satisfied, it is returned that the adaptive step sizing failed.

## Example

Example can be found at 
- examples/euler_test.cpp
- examples/adapt_rk_test.cpp
