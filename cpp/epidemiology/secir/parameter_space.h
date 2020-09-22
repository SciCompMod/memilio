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
void set_params_distributions_normal(SecirParams& params, double t0, double tmax, double dev_rel);

/* Draws a sample from SecirParams parameter distributions and stores sample values
* as SecirParams parameter values (cf. UncertainValue and SecirParams classes)
* @param[inout] params SecirParams including contact patterns for alle age groups
*/
void draw_sample(SecirParams& params);

} // namespace epi

#endif // PARAMETER_SPACE_H
