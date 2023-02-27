#include "ad/ad.hpp"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "boost/numeric/odeint/stepper/controlled_runge_kutta.hpp"
//#include "memilio/math/eigen.h"
#include <iostream>

using ad_type = typename ad::gt1s<double>::type;

/* The type of container used to hold the state vector */
typedef std::vector< ad_type > state_type;

double gam = 0.15;

/* The rhs of x' = f(x) */
void harmonic_oscillator( const state_type &x , state_type &dxdt , const double /* t */ )
{
    dxdt[0] = x[1];
    dxdt[1] = -x[0] - gam*x[1];
}

using error_stepper_type = boost::numeric::odeint::runge_kutta_cash_karp54< state_type >;
using controlled_stepper_type = boost::numeric::odeint::controlled_runge_kutta<error_stepper_type>;

int main() {
  state_type x(2);
  x[0] = 1.0; // start at x=1.0, p=0.0
  x[1] = 0.0;
  controlled_stepper_type controlled_stepper;
  //boost::numeric::odeint::integrate_adaptive( controlled_stepper , harmonic_oscillator , x , 0.0 , 10.0 , 0.01 );
  return 0;
}
