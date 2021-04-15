#ifndef PARAMETER_SPACE_H
#define PARAMETER_SPACE_H

#include "epidemiology/utils/memory.h"
#include "epidemiology/utils/logging.h"
#include "epidemiology/utils/parameter_distributions.h"
#include "epidemiology/secir/secir.h"

#include <assert.h>
#include <string>
#include <vector>
#include <random>
#include <memory>

namespace epi
{
/* Sets alls SecirParams parameters normally distributed, 
*  using the current value and a given standard deviation
* @param[inout] params SecirParams including contact patterns for alle age groups
* @param[in] t0 start time
* @param[in] tmax end time
* @param[in] dev_rel maximum relative deviation from particular value(s) given in params
*/
void set_params_distributions_normal(
    SecirModel& model, double t0,
    double tmax, double dev_rel);

/**
 * draws a sample from the specified distributions for all parameters related to the demographics, e.g. population.
 * @param[inout] model SecirModel including contact patterns for alle age groups
 */
void draw_sample_demographics(SecirModel& model);

/**
 * draws a sample from the specified distributions for all parameters related to the infection.
 * @param[inout] model SecirModel including contact patterns for alle age groups
 */
void draw_sample_infection(SecirModel& model);

/** Draws a sample from SecirModel parameter distributions and stores sample values
* as SecirParams parameter values (cf. UncertainValue and SecirParams classes)
* @param[inout] model SecirModel including contact patterns for alle age groups
*/
void draw_sample(SecirModel& model);

} // namespace epi

#endif // PARAMETER_SPACE_H
