#include <epidemiology/utils/parameter_distributions.h>
#include <epidemiology/secir/parameter_space.h>
#include <epidemiology/secir/secir.h>

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
    epi::ParameterDistributionNormal some_parameter(min, max, mean, stddev);
    // "some parameter",
    // std::make_unique<epi::ParameterDistributionNormal>(epi::ParameterDistributionNormal(min, max, mean, stddev))};

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
    epi::ParameterDistributionUniform some_other_parameter(1.0, 10.0);

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
    epi::SecirModel<epi::AgeGroup3> model;
    auto& params = model.parameters;

    size_t nb_groups = params.get_num_groups();
    epi::ContactMatrixGroup cm_group{epi::ContactMatrix(Eigen::MatrixXd::Constant(nb_groups, nb_groups, 0.5))};
    params.get_contact_patterns() = cm_group;

    double t0   = 0;
    double tmax = 100;

    params.get_contact_patterns().set_distribution_damp_nb(epi::ParameterDistributionUniform(2, (tmax - t0) / 10));
    params.get_contact_patterns().set_distribution_damp_days(epi::ParameterDistributionUniform(t0, tmax));
    params.get_contact_patterns().set_distribution_damp_diag_base(epi::ParameterDistributionUniform(0.1, 1));
    params.get_contact_patterns().set_distribution_damp_diag_rel(epi::ParameterDistributionUniform(0.6, 1.4));
    params.get_contact_patterns().set_distribution_damp_offdiag_rel(epi::ParameterDistributionUniform(0.7, 1.1));

    draw_sample(model);
    auto& cfmat_sample = params.get_contact_patterns().get_cont_freq_mat();

    printf("\n\n Number of dampings: %zu\n", cfmat_sample[0].get_dampings().size());

    printf("\n First damping at %.2f with factor %.2f\n",
           double(cfmat_sample[0].get_dampings()[0].get_time()),
           cfmat_sample[0].get_dampings()[0].get_coeffs()(0, 0));

    // printout the second damping
    printf("\n Damping at day %.2f\n\t", double(cfmat_sample[0].get_dampings()[1].get_time()));
    std::cout << cfmat_sample[0].get_dampings()[1].get_coeffs() << std::endl;
}
