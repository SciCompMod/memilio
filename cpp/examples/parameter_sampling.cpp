#include <epidemiology/parameter_studies/parameter_space.h>
#include <epidemiology/secir.h>
#include <stdio.h>

// verify manually that the correct distribution is chosen
int main()
{
    /*
     * Real valued variable element
     */
    double mean   = 5;
    double stddev = 1.5;
    double min    = 1;
    double max    = 10;
    // check if constructor works correctly
    epi::RealVariableElement some_parameter{"some parameter",
                                            new epi::ParameterDistributionNormal(min, max, mean, stddev)};

    printf("\n N(%.0f,%.0f)-distribution with sampling only in [%.0f,%.0f]", mean, stddev, min, max);
    int counter[10] = {0};
    for (int i = 0; i < 1000; i++) {
        int rounded = (int)(some_parameter.get_sample() - 1);
        if (rounded >= 0 && rounded < 10) {
            counter[rounded]++;
        }
    }
    double acc = 0;
    for (int i = 0; i < 9; i++) {
        acc += (double)counter[i] / 1000.0;
        printf("\n [%d-%d): %.2f %.2f ", i + 1, i + 2, (double)counter[i] / 1000.0, acc);
    }
    printf("\n");

    // check if constructor works correctly
    printf("\n U(%.0f,%.0f)-distribution", min, max);
    epi::RealVariableElement some_other_parameter{"some parameter", new epi::ParameterDistributionUniform(1.0, 10.0)};

    double counter_unif[10] = {0};
    for (int i = 0; i < 1000; i++) {
        int rounded = (int)(some_other_parameter.get_sample() - 1);
        if (rounded >= 0 && rounded < 10) {
            counter_unif[rounded]++;
        }
    }
    acc = 0;
    for (int i = 0; i < 9; i++) {
        acc += (double)counter_unif[i] / 1000.0;
        printf("\n [%d-%d): %.2f %.2f ", i + 1, i + 2, (double)counter_unif[i] / 1000.0, acc);
    }

    /*
     * Contact frequency and dampings variable element
     */
    size_t nb_groups = 3;
    epi::ContactFrequencyMatrix contact_freq_matrix{nb_groups};
    for (size_t i = 0; i < nb_groups; i++) {
        for (size_t j = i; i < nb_groups; i++) {
            contact_freq_matrix.set_cont_freq(0.5, i, j);
        }
    }

    double t0   = 0;
    double tmax = 10;
    epi::ContactFrequencyVariableElement contact_varel{contact_freq_matrix,
                                                       epi::ParameterDistributionUniform(1, (tmax - t0) / 10),
                                                       epi::ParameterDistributionUniform(t0, tmax),
                                                       epi::ParameterDistributionUniform(0.1, 1),
                                                       epi::ParameterDistributionUniform(0.6, 1.4),
                                                       epi::ParameterDistributionUniform(0.7, 1.1)};

    epi::ContactFrequencyMatrix cfmat_sample = contact_varel.get_sample();

    for (size_t i = 0; i < nb_groups; i++) {
        for (size_t j = i; i < nb_groups; i++) {
            cfmat_sample.get_dampings(i, j);
        }
    }
}