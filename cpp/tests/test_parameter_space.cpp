#include <gtest/gtest.h>
#include <epidemiology/seirParam.h>
#include <epidemiology/parameter_studies/parameter_space.h>

TEST (parameterspace, init_from_seir_params)
{
    struct seirParam<double> params;
    epi::parameter_space_t parameter_space(params, 0.1);
}