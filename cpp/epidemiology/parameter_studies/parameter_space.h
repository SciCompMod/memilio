#ifndef PARAMETER_SPACE_H
#define PARAMETER_SPACE_H

#include <assert.h>
#include <string>
#include <vector>
#include <random>
#include <memory>
#include <epidemiology/memory.h>
#include <epidemiology/secir.h>
#include <epidemiology/logging.h>
#include <epidemiology/parameter_studies/parameter_distributions.h>

namespace epi
{

/* Creates a parameter space from SecirParams values; using normal distributions
* @param[inout] params SecirParams including contact patterns for alle age groups
* @param[in] t0 start time
* @param[in] tmax end time
* @param[in] dev_rel maximum relative deviation from particular value(s) given in params
*/
void create_param_space_normal(SecirParams& params, double t0, double tmax, double dev_rel);

/* Draws a sample from SecirParams parameter distributions and stores sample values
* as SecirParams parameter values (cf. UncertainValue and SecirParams classes)
* @param[inout] params SecirParams including contact patterns for alle age groups
*/
void draw_sample(SecirParams& params);

} // namespace epi

#endif // PARAMETER_SPACE_H
