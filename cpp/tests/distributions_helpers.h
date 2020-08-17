#ifndef DISTRIBUTIONS_HELPERS_H
#define DISTRIBUTIONS_HELPERS_H

#include <epidemiology/parameter_studies/parameter_studies.h>
#include <gtest/gtest.h>

void check_dist(const epi::ParameterDistribution& dist, const epi::ParameterDistribution& dist_read);

#endif