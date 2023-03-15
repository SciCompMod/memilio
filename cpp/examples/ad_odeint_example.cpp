#include "ad/ad.hpp"

#include "memilio/utils/compiler_diagnostics.h"

/* disable warnings of odint.hpp */

GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wdeprecated-declarations")

#include "boost/numeric/odeint.hpp" // IWYU pragma: keep

GCC_CLANG_DIAGNOSTIC(pop)


using ad_forward_t = typename ad::gt1s<double>::type;
using ad_adjoint_t = typename ad::ga1s<double>::type;

/* The type of container used to hold the state vector */
using value_type = ad_forward_t;
using time_type = value_type;
typedef std::vector< value_type > state_type;



double gam = 0.15;

/* The rhs of x' = f(x) */
void harmonic_oscillator( const state_type &x , state_type &dxdt , const time_type /* t */ )
{
    dxdt[0] = x[1];
    dxdt[1] = -x[0] - gam*x[1];
}

using error_stepper_type = boost::numeric::odeint::runge_kutta_cash_karp54< state_type, value_type,
    state_type, value_type >;
using controlled_stepper_type = boost::numeric::odeint::controlled_runge_kutta<error_stepper_type>;

using boost::numeric::odeint::make_controlled;

int main() {
  state_type x(2);
  x[0] = 1.0; // start at x=1.0, p=0.0
  x[1] = 0.0;
  ad::derivative(x[0]) = 1.0;
  controlled_stepper_type controlled_stepper;
  boost::numeric::odeint::integrate_adaptive(make_controlled<error_stepper_type>(1e-6,1e-6),
                                             harmonic_oscillator , x , time_type( 0.0) , time_type(10.0) , time_type(0.01) );
  std::cout << "ad derivative of x is (" << ad::derivative(x[0]) << ", " << ad::derivative(x[1]) << ")\n";
  const double h = 1e-3;
  std::vector<double> y = {ad::value(x[0]), ad::value(x[1])};
  x[0] = 1.0 + h;
  x[1] = 0.0;
  boost::numeric::odeint::integrate_adaptive(make_controlled<error_stepper_type>(1e-6,1e-6),
                                             harmonic_oscillator , x , time_type( 0.0) , time_type(10.0) , time_type(0.01) );


  std::cout << "finite differences derivative of x is (" << (ad::value(x[0])-y[0])/h << ", "
            << (ad::value(x[1])-y[1])/h << ")\n";

  return 0;
}
