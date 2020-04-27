#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include <vector>
#include <functional>

/**
 * Function template to be integrated
 */
template <typename T>
using DerivFunction = std::function<void(std::vector<T> const &y, const T t, std::vector<T> &dydt)>;


#endif // INTEGRATOR_H
