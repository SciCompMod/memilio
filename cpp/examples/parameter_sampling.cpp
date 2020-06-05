#include <epidemiology/parameter_studies/parameter_space.h>
#include <epidemiology/secir.h>
#include <stdio.h>

// verify manually that the correct distribution is chosen
int main()
{
    double mean   = 5;
    double stddev = 1.5;
    double min    = 1;
    double max    = 10;
    // check if full argument constructor works correctly
    epi::ParameterDistributionNormal parameter_dist_normal("dummy", min, max, mean, stddev);

    printf("\n N(%.0f,%.0f)-distribution with sampling only in [%.0f,%.0f]", mean, stddev, min, max);
    int counter[10] = {0};
    for (int i = 0; i < 1000; i++) {
        int rounded = (int)(parameter_dist_normal.get_sample_point() - 1);
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

    // check if full argument constructor works correctly
    printf("\n U(%.0f,%.0f)-distribution", min, max);
    epi::ParameterDistributionUniform parameter_dist_unif("dummy", 1.0, 10.0);

    double counter_unif[10] = {0};
    for (int i = 0; i < 1000; i++) {
        int rounded = (int)(parameter_dist_unif.get_sample_point() - 1);
        if (rounded >= 0 && rounded < 10) {
            counter_unif[rounded]++;
        }
    }
    acc = 0;
    for (int i = 0; i < 9; i++) {
        acc += (double)counter_unif[i] / 1000.0;
        printf("\n [%d-%d): %.2f %.2f ", i + 1, i + 2, (double)counter_unif[i] / 1000.0, acc);
    }
}