#include <epidemiology/parameter_studies/parameter_space.h>
#include <epidemiology/secir.h>
#include <gtest/gtest.h>

TEST(parameterspace, init_from_seir_params) {
  epi::SecirParams params;
  epi::parameter_space_t parameter_space(params, 0.1);
}