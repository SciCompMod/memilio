#ifndef SECIR_H
#define SECIR_H

#include <vector>

#include <epidemiology/seir_param.h>

/**
 * Returns the damping factor
 *
 * @param[in] damping_array Array of dampings
 * @param[in] day Current day
 */
double getDampingFactor(std::vector<Damping> const& damping_array, double day);

/**
 * Computes the current time-derivative of S, E, I, and R in the SEIR model
 * @tparam T the datatype of the cases
 * @param[in] params SEIR Model parameters, created by seir_param
 * @param[in] y current  S, E, I, and R values at t; y: [0:S, 1:E, 2:I, 3:R]
 * @param[in] t time / current day
 * @param[out] dydt the values of the time derivatices of S, E, I, and R
 */
void seir_getDerivatives(SeirParams const& params, std::vector<double> const& y, double t, std::vector<double>& dydt);

/**
 * Computes the current time-derivative of S, E, I, and R in the SEIR model
 * @param[in] params SEIR Model parameters, created by seir_param
 * @tparam T the datatype of the cases
 * @param[in] y current  S, E, I, and R values at t; y: [0:S, 1:E, 2:I, 3:R]
 * @param[in] t time / current day
 * @param[out] dydt the values of the time derivatices of S, E, I, and R
 */
void secir_getDerivatives(SeirParams const& params, std::vector<double> const& y, double t, std::vector<double>& dydt);

/**
 * Computes the seir curve by integration
 * @param[in] seir_0 Initial S, E, I, and R values at t0
 * @param[in] t0 start time of simulation
 * @param[in] tmax end time of simulation
 * @param[in] dt initial time step
 * @param[in] params SEIR model parameters
 *
 * @returns Vector of times t
 */
std::vector<double> simulate(double t0, double tmax, double dt, SeirParams const& params,
                             std::vector<std::vector<double>>& seir);
#endif // SECIR_H
