#include "matchers.h"
#include "test_data_dir.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/mobility/mobility.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/parameter_distributions.h"
#include "memilio/utils/parameter_set.h"
#include "memilio/utils/random_number_generator.h"
#include "memilio/utils/uncertain_value.h"
#include "ode_secirvvs/infection_state.h"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/parameter_space.h"
#include "ode_secirvvs/parameters.h"
#include "ode_secirvvs/parameters_io.h"

#include "gtest/gtest.h"
#include "gmock/gmock-matchers.h"
#include <gmock/gmock-generated-matchers.h>
#include <limits>

void assign_uniform_distribution(mio::UncertainValue& p, double min, double max, bool set_invalid_initial_value)
{
    auto invalid_initial = max == 0 ? 1.0 : max * 1.001;
    auto valid_initial   = (max + min) * 0.5;
    auto initial         = set_invalid_initial_value ? invalid_initial : valid_initial;
    p                    = mio::UncertainValue(initial);
    p.set_distribution(mio::ParameterDistributionUniform(min, max));
}

template <size_t N>
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>& array,
                                       const double (&min)[N], const double (&max)[N], bool set_invalid_initial_value)
{
    assert(N == array.numel());
    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(N); ++i) {
        assign_uniform_distribution(array[i], min[size_t(i)], max[size_t(i)], set_invalid_initial_value);
    }
}

void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue, mio::AgeGroup>& array, double min,
                                       double max, bool set_invalid_initial_value)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max, set_invalid_initial_value);
    }
}

void set_synthetic_population_data(mio::osecirvvs::Model::Populations& populations, bool set_invalid_initial_value)
{
    for (mio::AgeGroup i = 0; i < mio::get<mio::AgeGroup>(populations.size()); i++) {
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::ExposedNaive}], 10, 20,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::ExposedImprovedImmunity}], 10, 20,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::ExposedPartialImmunity}], 10, 20,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::CarrierNaive}], 5, 10,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::CarrierPartialImmunity}], 5, 10,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::CarrierImprovedImmunity}], 5, 10,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedNaive}], 5, 10,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedPartialImmunity}], 5, 10,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::InfectedImprovedImmunity}], 5, 10,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::HospitalizedNaive}], 1, 2,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::HospitalizedImprovedImmunity}], 1,
                                    2, set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::HospitalizedPartialImmunity}], 1, 2,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::ICUNaive}], 1, 2,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::ICUPartialImmunity}], 1, 2,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::ICUImprovedImmunity}], 1, 2,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::Recovered}], 200, 300,
                                    set_invalid_initial_value);
        assign_uniform_distribution(populations[{i, mio::osecirvvs::InfectionState::SusceptiblePartialImmunity}], 200,
                                    300, set_invalid_initial_value);
        populations.set_difference_from_group_total<mio::AgeGroup>(
            {i, mio::osecirvvs::InfectionState::SusceptibleNaive}, 1000);
    }
}

void set_demographic_parameters(mio::osecirvvs::Model::ParameterSet& parameters, bool set_invalid_initial_value)
{
    assign_uniform_distribution(parameters.get<mio::osecirvvs::ICUCapacity>(), 20, 50, set_invalid_initial_value);
    assign_uniform_distribution(parameters.get<mio::osecirvvs::TestAndTraceCapacity>(), 100, 200,
                                set_invalid_initial_value);
}

void set_covid_parameters(mio::osecirvvs::Model::ParameterSet& params, bool set_invalid_initial_value)
{
    //times
    const double tinc             = 5.2; // R_2^(-1)+R_3^(-1)
    const double tserint_min      = 0.5 * 2.67 + 0.5 * 5.2; // R_2^(-1)+0.5*R_3^(-1)
    const double tserint_max      = 0.5 * 4.00 + 0.5 * 5.2;
    const double t_inf_rec_min    = 5.6; // R4^(-1) = T_I^R
    const double t_inf_rec_max    = 8.4;
    const double t_inf_hosp_min[] = {9, 9, 9, 5, 5, 5}; // R6^(-1) = T_I^H
    const double t_inf_hosp_max[] = {12, 12, 12, 7, 7, 7};
    const double t_hosp_rec_min[] = {4, 4, 5, 7, 9, 13}; // R5^(-1) = T_H^R
    const double t_hosp_rec_max[] = {6, 6, 7, 9, 11, 17};
    const double t_hosp_icu_min   = 3; // R7^(-1) = T_H^U
    const double t_hosp_icu_max   = 7;
    const double t_icu_rec_min[]  = {5, 5, 5, 14, 14, 10}; // R8^(-1) = T_U^R
    const double t_icu_rec_max[]  = {9, 9, 9, 21, 21, 15};
    const double t_icu_dead_min[] = {4, 4, 4, 15, 15, 10}; // 5-16 (=R8^(-1) = T_U^R)
    const double t_icu_dead_max[] = {8, 8, 8, 18, 18, 12};

    array_assign_uniform_distribution(params.get<mio::osecirvvs::IncubationTime>(), tinc, tinc,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::SerialInterval>(), tserint_min, tserint_max,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::InfectiousTimeMild>(), t_inf_rec_min, t_inf_rec_max,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::HomeToHospitalizedTime>(), t_inf_hosp_min,
                                      t_inf_hosp_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::HospitalizedToHomeTime>(), t_hosp_rec_min,
                                      t_hosp_rec_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::HospitalizedToICUTime>(), t_hosp_icu_min,
                                      t_hosp_icu_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ICUToHomeTime>(), t_icu_rec_min, t_icu_rec_max,
                                      set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ICUToDeathTime>(), t_icu_dead_min, t_icu_dead_max,
                                      set_invalid_initial_value);

    //probabilities
    double fac_variant                   = 1.4;
    const double transmission_risk_min[] = {0.02 * fac_variant, 0.05 * fac_variant, 0.05 * fac_variant,
                                            0.05 * fac_variant, 0.08 * fac_variant, 0.1 * fac_variant};

    const double transmission_risk_max[] = {0.04 * fac_variant, 0.07 * fac_variant, 0.07 * fac_variant,
                                            0.07 * fac_variant, 0.10 * fac_variant, 0.15 * fac_variant};
    const double carr_infec_min          = 0.5;
    const double carr_infec_max          = 0.5;
    const double beta_low_incidenc_min   = 0.0; // beta (depends on incidence and test and trace capacity)
    const double beta_low_incidenc_max   = 0.2;
    const double beta_high_incidence_min = 0.4;
    const double beta_high_incidence_max = 0.5;
    const double prob_car_rec_min[]      = {0.2, 0.2, 0.15, 0.15, 0.15, 0.15}; // alpha
    const double prob_car_rec_max[]      = {0.3, 0.3, 0.25, 0.25, 0.25, 0.25};
    const double prob_inf_hosp_min[]     = {0.006, 0.006, 0.015, 0.049, 0.15, 0.20}; // rho
    const double prob_inf_hosp_max[]     = {0.009, 0.009, 0.023, 0.074, 0.18, 0.25};
    const double prob_hosp_icu_min[]     = {0.05, 0.05, 0.05, 0.10, 0.25, 0.35}; // theta
    const double prob_hosp_icu_max[]     = {0.10, 0.10, 0.10, 0.20, 0.35, 0.45};
    const double prob_icu_dead_min[]     = {0.00, 0.00, 0.10, 0.10, 0.30, 0.5}; // delta
    const double prob_icu_dead_max[]     = {0.10, 0.10, 0.18, 0.18, 0.50, 0.7};

    const double reduc_vacc_exp_min    = 0.75;
    const double reduc_vacc_exp_max    = 0.85;
    const double reduc_immune_exp_min  = 0.281;
    const double reduc_immune_exp_max  = 0.381;
    const double reduc_vacc_inf_min    = 0.6;
    const double reduc_vacc_inf_max    = 0.7;
    const double reduc_immune_inf_min  = 0.193;
    const double reduc_immune_inf_max  = 0.293;
    const double reduc_vacc_hosp_min   = 0.05;
    const double reduc_vacc_hosp_max   = 0.15;
    const double reduc_immune_hosp_min = 0.041;
    const double reduc_immune_hosp_max = 0.141;
    const double reduc_mild_rec_time   = 1.0;

    array_assign_uniform_distribution(params.get<mio::osecirvvs::InfectionProbabilityFromContact>(),
                                      transmission_risk_min, transmission_risk_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RelativeCarrierInfectability>(), carr_infec_min,
                                      carr_infec_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RiskOfInfectionFromSympomatic>(),
                                      beta_low_incidenc_min, beta_low_incidenc_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::MaxRiskOfInfectionFromSympomatic>(),
                                      beta_high_incidence_min, beta_high_incidence_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::AsymptoticCasesPerInfectious>(), prob_car_rec_min,
                                      prob_car_rec_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::HospitalizedCasesPerInfectious>(), prob_inf_hosp_min,
                                      prob_inf_hosp_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ICUCasesPerHospitalized>(), prob_hosp_icu_min,
                                      prob_hosp_icu_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::DeathsPerICU>(), prob_icu_dead_min, prob_icu_dead_max,
                                      set_invalid_initial_value);

    array_assign_uniform_distribution(params.get<mio::osecirvvs::ExposedFactorPartialImmunity>(), reduc_vacc_exp_min,
                                      reduc_vacc_exp_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ExposedFactorImprovedImmunity>(), reduc_immune_exp_min,
                                      reduc_immune_exp_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::InfectedFactorPartialImmunity>(), reduc_vacc_inf_min,
                                      reduc_vacc_inf_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::InfectedFactorImprovedImmunity>(),
                                      reduc_immune_inf_min, reduc_immune_inf_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::HospitalizedFactorPartialImmunity>(),
                                      reduc_vacc_hosp_min, reduc_vacc_hosp_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::HospitalizedFactorImprovedImmunity>(),
                                      reduc_immune_hosp_min, reduc_immune_hosp_max, set_invalid_initial_value);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::InfectiousTimeFactorImmune>(), reduc_mild_rec_time,
                                      reduc_mild_rec_time, set_invalid_initial_value);

    //sasonality
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.3;

    assign_uniform_distribution(params.get<mio::osecirvvs::Seasonality>(), seasonality_min, seasonality_max,
                                set_invalid_initial_value);
}

mio::osecirvvs::Model make_model(int num_age_groups, bool set_invalid_initial_value = false)
{
    mio::osecirvvs::Model model(num_age_groups);
    set_covid_parameters(model.parameters, set_invalid_initial_value);
    set_synthetic_population_data(model.populations, set_invalid_initial_value);
    set_demographic_parameters(model.parameters, set_invalid_initial_value);
    return model;
}

TEST(TestOdeSECIRVVS, draw_sample)
{
    mio::log_thread_local_rng_seeds(mio::LogLevel::warn);

    mio::Graph<mio::osecirvvs::Model, mio::MigrationParameters> graph;

    auto num_age_groups = 6;
    //create model with invalid initials so the test fails if no sampling is done
    graph.add_node(0, make_model(num_age_groups, /*set_invalid_initial_value*/true));
    graph.add_node(0, make_model(num_age_groups, /*set_invalid_initial_value*/true));
    graph.add_edge(0, 1, Eigen::VectorXd::Constant(num_age_groups, num_age_groups));

    auto sampled_graph = mio::osecirvvs::draw_sample(graph, true);

    ASSERT_EQ(sampled_graph.nodes().size(), graph.nodes().size());
    ASSERT_EQ(sampled_graph.edges().size(), graph.edges().size());

    // spot check for sampling
    auto& parameters0 = sampled_graph.nodes()[0].property.parameters;
    auto& populations0 = sampled_graph.nodes()[0].property.populations;
    ASSERT_GE(double(parameters0.get<mio::osecirvvs::ICUToHomeTime>()[mio::AgeGroup(1)]), 5.0);
    ASSERT_LE(double(parameters0.get<mio::osecirvvs::ICUToHomeTime>()[mio::AgeGroup(1)]), 9.0);
    ASSERT_GE(double(parameters0.get<mio::osecirvvs::ExposedFactorPartialImmunity>()[mio::AgeGroup(0)]), 0.75);
    ASSERT_LE(double(parameters0.get<mio::osecirvvs::ExposedFactorPartialImmunity>()[mio::AgeGroup(0)]), 0.85);
    ASSERT_GE(double(populations0[{mio::AgeGroup(2), mio::osecirvvs::InfectionState::InfectedPartialImmunity}]), 5.0);
    ASSERT_LE(double(populations0[{mio::AgeGroup(2), mio::osecirvvs::InfectionState::InfectedPartialImmunity}]), 10.0);

    // special cases
    ASSERT_DOUBLE_EQ(populations0.get_total(), 1000 * num_age_groups);
    ASSERT_TRUE((parameters0.get<mio::osecirvvs::BaseInfectiousnessB161>().array(),
                 parameters0.get<mio::osecirvvs::InfectionProbabilityFromContact>().array() * 1.6) //using high variant
                    .all());

    // spot check for parameters that should be equal or different between nodes
    auto& parameters1  = sampled_graph.nodes()[1].property.parameters;
    auto& populations1 = sampled_graph.nodes()[1].property.populations;
    ASSERT_TRUE((parameters1.get<mio::osecirvvs::InfectiousTimeMild>().array() ==
                 parameters0.get<mio::osecirvvs::InfectiousTimeMild>().array())
                    .all());
    ASSERT_NE(parameters1.get<mio::osecirvvs::ICUCapacity>(), parameters0.get<mio::osecirvvs::ICUCapacity>());
    ASSERT_FALSE((populations1.array() == populations0.array()).all());
}

TEST(TestOdeSECIRVVS, read_data)
{
    auto num_age_groups = 6;
    auto model = std::vector<mio::osecirvvs::Model>({make_model(num_age_groups)});

    auto read_result = mio::osecirvvs::read_input_data_county(
        model, {2020, 12, 01}, {1002}, std::vector<double>(size_t(num_age_groups), 1.0), 1.0, TEST_DATA_DIR, 10);
    ASSERT_THAT(read_result, IsSuccess());

    auto expected_values =
        (Eigen::ArrayXd(num_age_groups * Eigen::Index(mio::osecirvvs::InfectionState::Count)) << 8796.3, 499.898,
         2.98826, 0.167206, 0.155258, 1.13759, 0.105142, 0.259947, 0, 0, 0, 6.86457, 0.252844, 0.191817, 0, 0, 0,
         0.0158157, 1.38274e-05, 6.21083e-05, 0.0292738, 0.000163712, 0.000302092, 1012.63, 0, 0, 17075.8, 482.905,
         7.22816, 0.201271, 0.309686, 4.32397, 0.185362, 0.743502, 0, 0, 0, 11.684, 0.214395, 0.269443, 0, 0, 0,
         0.0497858, 2.1661e-05, 0.000161222, 0.0294097, 8.22282e-05, 0.000251431, 1621.04, 0, 0, 72508.6, 314.79,
         43.4754, 0.185779, 0.411666, 18.0428, 0.122609, 0.715364, 0, 0, 0, 65.2366, 0.182585, 0.330993, 0, 0, 0,
         1.26041, 8.41558e-05, 0.000902068, 0.0751831, 3.25529e-05, 0.000143349, 1533.58, 0, 0, 80801.2, 375.849,
         41.674, 0.190853, 0.427392, 17.9371, 0.131189, 0.774526, 0, 0, 0, 58.4369, 0.17792, 0.324688, 0, 0, 0, 3.66709,
         0.000262407, 0.00284255, 0.48662, 0.000225766, 0.00100471, 1829.35, 0, 0, 41585.3, 508.615, 9.2403, 0.111291,
         0.181928, 4.1603, 0.0790557, 0.339463, 0, 0, 0, 12.977, 0.105515, 0.140007, 0, 0, 0, 2.98513, 0.000561766,
         0.00444226, 2.60421, 0.00316907, 0.0102951, 1803.31, 0, 0, 13235.8, 413.078, 9.7863, 0.300906, 0.687777,
         1.0229, 0.0471198, 0.278325, 0, 0, 0, 3.41896, 0.0709641, 0.131662, 0, 0, 0, 0.669835, 0.000321811, 0.00355815,
         4.67857, 0.0143617, 0.0652351, 2049.13, 0, 0)
            .finished();

    ASSERT_THAT(print_wrap(model[0].populations.array().cast<double>()), MatrixNear(print_wrap(expected_values), 1e-5, 1e-5));
}
