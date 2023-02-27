#include "ad/ad.hpp"
#include "boost/numeric/odeint/stepper/runge_kutta_cash_karp54.hpp"
#include "memilio/math/eigen.h"
#include <iostream>

using ad_type = typename ad::gt1s<double>::type;

/* The type of container used to hold the state vector */
typedef std::vector< ad_type > state_type;

ad_type gam;

/* The rhs of x' = f(x) */
void harmonic_oscillator( const state_type &x , state_type &dxdt , const double /* t */ )
{
    dxdt[0] = x[1];
    dxdt[1] = -x[0] - gam*x[1];
}

int main() {

  return 0;
}
