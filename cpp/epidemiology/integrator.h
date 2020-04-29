#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <functional>

/**
 * Function template to be integrated
 */
using DerivFunction = std::function<void(std::vector<double> const& y, const double t, std::vector<double>& dydt)>;

#endif // INTEGRATOR_H
