/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Henrik Zunker
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include "ode_secirvvs/analyze_result.h"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/parameters.h"
#include "ode_secirvvs/parameters_io.h"
#include "ode_secirvvs/parameter_space.h"
#include "memilio/io/result_io.h"
#include "memilio/io/mobility_io.h"
#include "memilio/compartments/parameter_studies.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/compartments/simulation.h"
#include "memilio/utils/logging.h"

#include <iostream>
#include <random>

/**
 * Set a value and distribution of an UncertainValue.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution.
 * @param p uncertain value to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void assign_uniform_distribution(mio::UncertainValue<double>& p, double min, double max)
{
    p = mio::UncertainValue<double>(0.5 * (max + min));
    p.set_distribution(mio::ParameterDistributionUniform(min, max));
}

/**
 * Set a value and distribution of an array of UncertainValues.
 * Assigns average of min[i] and max[i] as a value and UNIFORM(min[i], max[i]) as a distribution for
 * each element i of the array.
 * @param array array of UncertainValues to set.
 * @param min minimum of distribution for each element of array.
 * @param max minimum of distribution for each element of array.
 */
template <size_t N>
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>& array,
                                       const double (&min)[N], const double (&max)[N])
{
    assert(N == array.numel());
    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(N); ++i) {
        assign_uniform_distribution(array[i], min[size_t(i)], max[size_t(i)]);
    }
}

/**
 * Set a value and distribution of an array of UncertainValues.
 * Assigns average of min and max as a value and UNIFORM(min, max) as a distribution to every element of the array.
 * @param array array of UncertainValues to set.
 * @param min minimum of distribution.
 * @param max minimum of distribution.
 */
void array_assign_uniform_distribution(mio::CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>& array,
                                       double min, double max)
{
    for (auto i = mio::AgeGroup(0); i < array.size<mio::AgeGroup>(); ++i) {
        assign_uniform_distribution(array[i], min, max);
    }
}

/**
 * Set epidemiological parameters of Covid19.
 * @param params Object that the parameters will be added to.
 * @returns Currently generates no errors.
 */
mio::IOResult<void> set_covid_parameters(mio::osecirvvs::Parameters<double>& params, ScalarType fac_variant,
                                         ScalarType red_fact_exposed)
{
    //times
    // doi.org/10.1016/j.lanepe.2022.100446 , doi.org/10.3201/eid2806.220158
    const double timeExposedMin            = 1.66;
    const double timeExposedMax            = 1.66;
    const double timeInfectedNoSymptomsMin = 1.44;
    const double timeInfectedNoSymptomsMax = 1.44;

    const double timeInfectedSymptomsMin = 6.58; //https://doi.org/10.1016/S0140-6736(22)00327-0
    const double timeInfectedSymptomsMax = 7.16; //https://doi.org/10.1016/S0140-6736(22)00327-0
    const double timeInfectedSevereMin[] = {1.8, 1.8, 1.8, 2.5, 3.5, 4.91}; // doi.org/10.1186/s12879-022-07971-6
    const double timeInfectedSevereMax[] = {2.3, 2.3, 2.3, 3.67, 5, 7.01}; // doi.org/10.1186/s12879-022-07971-6

    const double timeInfectedCriticalMin[] = {9.29,   9.29,  9.29,
                                              10.842, 11.15, 11.07}; // https://doi.org/10.1186/s12879-022-07971-6
    const double timeInfectedCriticalMax[] = {10.57, 10.57, 10.57,
                                              12.86, 13.23, 13.25}; // https://doi.org/10.1186/s12879-022-07971-6

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeExposed<double>>(), timeExposedMin,
                                      timeExposedMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedNoSymptoms<double>>(),
                                      timeInfectedNoSymptomsMin, timeInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSymptoms<double>>(),
                                      timeInfectedSymptomsMin, timeInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedSevere<double>>(), timeInfectedSevereMin,
                                      timeInfectedSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::TimeInfectedCritical<double>>(),
                                      timeInfectedCriticalMin, timeInfectedCriticalMax);

    //probabilities
    const double transmissionProbabilityOnContactMin[] = {0.02 * fac_variant, 0.05 * fac_variant, 0.05 * fac_variant,
                                                          0.05 * fac_variant, 0.08 * fac_variant, 0.1 * fac_variant};

    const double transmissionProbabilityOnContactMax[] = {0.04 * fac_variant, 0.07 * fac_variant, 0.07 * fac_variant,
                                                          0.07 * fac_variant, 0.10 * fac_variant, 0.15 * fac_variant};

    const double relativeTransmissionNoSymptomsMin = 0.5;

    //{0.6, 0.55, 0.65,0.7, 0.75, 0.85}; // DOI: 10.1097/INF.0000000000003791
    const double relativeTransmissionNoSymptomsMax = 0.5;
    // {0.8, 0.75,  0.8,0.8, 0.825, 0.9}; // DOI: 10.1097/INF.0000000000003791
    const double riskOfInfectionFromSymptomaticMin    = 0.0; // beta (depends on incidence and test and trace capacity)
    const double riskOfInfectionFromSymptomaticMax    = 0.2;
    const double maxRiskOfInfectionFromSymptomaticMin = 0.4;
    const double maxRiskOfInfectionFromSymptomaticMax = 0.5;

    // DOI: 10.1097/INF.0000000000003791 geht hier auch. aber ähnliche werte
    const double recoveredPerInfectedNoSymptomsMin[] = {0.2, 0.25,  0.2,
                                                        0.2, 0.175, 0.1}; // doi.org/10.1101/2022.05.05.22274697
    const double recoveredPerInfectedNoSymptomsMax[] = {0.4, 0.45, 0.35,
                                                        0.3, 0.25, 0.15}; // doi.org/10.1101/2022.05.05.22274697

    // 56% weniger risiko ins krankenhaus doi:10.1136/bmjgh-2023-0123
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10347449/pdf/bmjgh-2023-012328.pdf
    // alternativ: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9321237/pdf/vaccines-10-01049.pdf

    // Faktoren aus https://www.thelancet.com/journals/lancet/article/PIIS0140-6736(22)00462-7/fulltext
    const double severePerInfectedSymptomsMin[] = {1 * 0.006,   0.8 * 0.006, 0.4 * 0.015,
                                                   0.3 * 0.049, 0.25 * 0.15, 0.35 * 0.2}; // 2021 paper
    const double severePerInfectedSymptomsMax[] = {1 * 0.009,   0.8 * 0.009, 0.4 * 0.023, 0.3 * 0.074,
                                                   0.25 * 0.18, 0.35 * 0.25}; // quelle 2021 paper + factors

    // const double criticalPerSevereMin[] = {
    //     0.0511, 0.0686, 0.0491, 0.114,
    //     0.1495, 0.0674}; // www.sozialministerium.at/dam/jcr:f472e977-e1bf-415f-95e1-35a1b53e608d/Factsheet_Coronavirus_Hospitalisierungen.pdf
    // const double criticalPerSevereMax[] = {
    //     0.0511, 0.0686, 0.0491, 0.114,
    //     0.1495, 0.0674}; // www.sozialministerium.at/dam/jcr:f472e977-e1bf-415f-95e1-35a1b53e608d/Factsheet_Coronavirus_Hospitalisierungen.pdf

    // delta paper
    // risk of icu admission https://doi.org/10.1177/14034948221108548
    const double fac_icu                = 0.05;
    const double criticalPerSevereMin[] = {0.05 * fac_icu, 0.05 * fac_icu, 0.05 * fac_icu,
                                           0.10 * fac_icu, 0.25 * fac_icu, 0.35 * fac_icu};
    const double criticalPerSevereMax[] = {0.10 * fac_icu, 0.10 * fac_icu, 0.10 * fac_icu,
                                           0.20 * fac_icu, 0.35 * fac_icu, 0.45 * fac_icu};

    // 61% weniger risiko zu sterben doi:10.1136/bmjgh-2023-0123
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10347449/pdf/bmjgh-2023-012328.pdf
    const double fac_dead               = 0.39;
    const double deathsPerCriticalMin[] = {fac_dead * 0.00, fac_dead * 0.00, fac_dead * 0.10,
                                           fac_dead * 0.10, fac_dead * 0.30, fac_dead * 0.5}; // 2021 paper
    const double deathsPerCriticalMax[] = {fac_dead * 0.10, fac_dead * 0.10, fac_dead * 0.18,
                                           fac_dead * 0.18, fac_dead * 0.50, fac_dead * 0.7};

    // alternative https://doi.org/10.1136/bmj-2022-070695
    // const double fac_dead_u59           = 0.14;
    // const double fac_dead_p59           = 0.44;
    // const double deathsPerCriticalMin[] = {fac_dead_u59 * 0.00, fac_dead_u59 * 0.00, fac_dead_u59 * 0.10,
    //                                        fac_dead_u59 * 0.10, fac_dead_p59 * 0.30, fac_dead_p59 * 0.5}; // 2021 paper
    // const double deathsPerCriticalMax[] = {fac_dead_u59 * 0.10, fac_dead_u59 * 0.10, fac_dead_u59 * 0.18,
    //                                        fac_dead_u59 * 0.18, fac_dead_p59 * 0.50, fac_dead_p59 * 0.7};

    const double reducExposedPartialImmunityMin = red_fact_exposed; //0.5; //0.569; // doi.org/10.1136/bmj-2022-071502
    const double reducExposedPartialImmunityMax = red_fact_exposed; //0.5; // 0.637; // doi.org/10.1136/bmj-2022-071502
    // const double reducExposedImprovedImmunityMin = 0.36; // https://doi.org/10.1038/s41591-023-02219-5
    // const double reducExposedImprovedImmunityMax = 0.66; // https://doi.org/10.1038/s41591-023-02219-5

    const double reducExposedImprovedImmunityMin = red_fact_exposed; //0.5;
    //0.34 * reducExposedPartialImmunityMin; // https://jamanetwork.com/journals/jama/fullarticle/2788487 0.19346
    const double reducExposedImprovedImmunityMax = red_fact_exposed; //0.5;
    //0.34 * reducExposedPartialImmunityMax; // https://jamanetwork.com/journals/jama/fullarticle/2788487 0.21658

    const double reducInfectedSymptomsPartialImmunityMin  = 0.746; // doi.org/10.1056/NEJMoa2119451
    const double reducInfectedSymptomsPartialImmunityMax  = 0.961; // doi.org/10.1056/NEJMoa2119451
    const double reducInfectedSymptomsImprovedImmunityMin = 0.295; // doi.org/10.1056/NEJMoa2119451
    const double reducInfectedSymptomsImprovedImmunityMax = 0.344; // doi.org/10.1056/NEJMoa2119451

    const double reducInfectedSevereCriticalDeadPartialImmunityMin =
        0.52; // www.assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1050721/Vaccine-surveillance-report-week-4.pdf
    const double reducInfectedSevereCriticalDeadPartialImmunityMax =
        0.82; // www.assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/1050721/Vaccine-surveillance-report-week-4.pdf
    const double reducInfectedSevereCriticalDeadImprovedImmunityMin = 0.1; // doi.org/10.1136/bmj-2022-071502
    const double reducInfectedSevereCriticalDeadImprovedImmunityMax = 0.19; // doi.org/10.1136/bmj-2022-071502
    const double reducTimeInfectedMild                              = 0.5; // doi.org/10.1101/2021.09.24.21263978

    array_assign_uniform_distribution(params.get<mio::osecirvvs::TransmissionProbabilityOnContact<double>>(),
                                      transmissionProbabilityOnContactMin, transmissionProbabilityOnContactMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RelativeTransmissionNoSymptoms<double>>(),
                                      relativeTransmissionNoSymptomsMin, relativeTransmissionNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RiskOfInfectionFromSymptomatic<double>>(),
                                      riskOfInfectionFromSymptomaticMin, riskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<double>>(),
                                      maxRiskOfInfectionFromSymptomaticMin, maxRiskOfInfectionFromSymptomaticMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::RecoveredPerInfectedNoSymptoms<double>>(),
                                      recoveredPerInfectedNoSymptomsMin, recoveredPerInfectedNoSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::SeverePerInfectedSymptoms<double>>(),
                                      severePerInfectedSymptomsMin, severePerInfectedSymptomsMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::CriticalPerSevere<double>>(), criticalPerSevereMin,
                                      criticalPerSevereMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::DeathsPerCritical<double>>(), deathsPerCriticalMin,
                                      deathsPerCriticalMax);

    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedPartialImmunity<double>>(),
                                      reducExposedPartialImmunityMin, reducExposedPartialImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducExposedImprovedImmunity<double>>(),
                                      reducExposedImprovedImmunityMin, reducExposedImprovedImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<double>>(),
                                      reducInfectedSymptomsPartialImmunityMin, reducInfectedSymptomsPartialImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<double>>(),
                                      reducInfectedSymptomsImprovedImmunityMin,
                                      reducInfectedSymptomsImprovedImmunityMax);
    array_assign_uniform_distribution(
        params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<double>>(),
        reducInfectedSevereCriticalDeadPartialImmunityMin, reducInfectedSevereCriticalDeadPartialImmunityMax);
    array_assign_uniform_distribution(
        params.get<mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<double>>(),
        reducInfectedSevereCriticalDeadImprovedImmunityMin, reducInfectedSevereCriticalDeadImprovedImmunityMax);
    array_assign_uniform_distribution(params.get<mio::osecirvvs::ReducTimeInfectedMild<double>>(),
                                      reducTimeInfectedMild, reducTimeInfectedMild);

    //sasonality
    const double seasonality_min = 0.1;
    const double seasonality_max = 0.3;

    assign_uniform_distribution(params.get<mio::osecirvvs::Seasonality<double>>(), seasonality_min, seasonality_max);

    // Delta specific parameter
    params.get<mio::osecirvvs::StartDayNewVariant>() = mio::get_day_in_year(mio::Date(2021, 6, 6));

    return mio::success();
}

/**
 * indices of contact matrix corresponding to locations where contacts occur.
 */
enum class ContactLocation
{
    Home = 0,
    School,
    Work,
    Other,
    Count,
};

static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}};

/**
 * Set contact matrices.
 * Reads contact matrices from files in the data directory.
 * @param data_dir data directory.
 * @param params Object that the contact matrices will be added to.
 * @returns any io errors that happen during reading of the files.
 */
mio::IOResult<void> set_contact_matrices(const fs::path& data_dir, mio::osecirvvs::Parameters<double>& params)
{
    //TODO: io error handling
    auto contact_matrices = mio::ContactMatrixGroup(contact_locations.size(), size_t(params.get_num_groups()));
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(
                              (data_dir / "contacts" / ("baseline_" + contact_location.second + ".txt")).string()));

        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline;
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::MatrixXd::Zero(6, 6);
    }
    params.get<mio::osecirvvs::ContactPatterns<double>>() = mio::UncertainContactMatrix<double>(contact_matrices);

    return mio::success();
}

/**
 * Create the input graph for the parameter study.
 * Reads files from the data directory.
 * @param start_date start date of the simulation.
 * @param end_date end date of the simulation.
 * @param data_dir data directory.
 * @returns created graph or any io errors that happen during reading of the files.
 */
mio::IOResult<Eigen::VectorXd> init_condition_germany(mio::Date start_date, mio::Date end_date,
                                                      const fs::path& data_dir)
{
    // global parameters
    const int num_age_groups = 6;
    mio::osecirvvs::Parameters<double> params(num_age_groups);
    params.get<mio::osecirvvs::StartDay>() = mio::get_day_in_year(start_date);
    BOOST_OUTCOME_TRY(set_covid_parameters(params, 1.0, 1.0));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, params));

    auto scaling_factor_infected = std::vector<double>(size_t(params.get_num_groups()), 1.0);
    auto scaling_factor_icu      = 1.0;
    auto tnt_capacity_factor     = 1.43 / 100000.;

    // graph of counties with populations and local parameters
    // and mobility between counties
    mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>> params_graph;
    const auto& read_function_nodes = mio::osecirvvs::read_input_data_county<mio::osecirvvs::Model<double>>;
    const auto& node_id_function    = mio::get_node_ids;

    const auto& set_node_function =
        mio::set_nodes<mio::osecirvvs::TestAndTraceCapacity<double>, mio::osecirvvs::ContactPatterns<double>,
                       mio::osecirvvs::Model<double>, mio::MobilityParameters<double>,
                       mio::osecirvvs::Parameters<double>, decltype(read_function_nodes), decltype(node_id_function)>;
    BOOST_OUTCOME_TRY(set_node_function(
        params, start_date, end_date, data_dir,
        mio::path_join((data_dir / "pydata" / "Germany").string(), "county_current_population.json"), true,
        params_graph, read_function_nodes, node_id_function, scaling_factor_infected, scaling_factor_icu,
        tnt_capacity_factor, mio::get_offset_in_days(end_date, start_date), false, true));

    // iterate over all Edges and add up the initial state which is gained by node.property.get_initial_values() and returns a Vector<double>
    Eigen::VectorXd result = params_graph.nodes()[0].property.get_initial_values();
    // set to zero
    result.setZero();
    for (auto& node : params_graph.nodes()) {
        result += node.property.get_initial_values();
    }

    // print the initial state where each entry is seperated by a comma
    std::cout << "Initial state: ";
    for (int i = 0; i < result.size(); i++) {
        std::cout << result[i];
        if (i < result.size() - 1) {
            std::cout << ",";
        }
    }
    std::cout << std::endl;

    return mio::success(result);
}

mio::IOResult<void> run(const fs::path& data_dir, const fs::path& result_dir)
{
    mio::set_log_level(mio::LogLevel::warn);
    double t0           = 0.0;
    auto t0_day         = mio::Date(2022, 2, 1);
    double tmax         = 1000;
    double dt           = 0.5;
    const auto num_runs = 500;

    mio::Graph<mio::osecirvvs::Model<double>, mio::MobilityParameters<double>> params_graph;
    mio::osecirvvs::Model<double> model(6);

    BOOST_OUTCOME_TRY(set_covid_parameters(model.parameters, 1.0, 1.0));
    BOOST_OUTCOME_TRY(set_contact_matrices(data_dir, model.parameters));

    const auto end_date = mio::offset_date_by_days(t0_day, int(std::ceil(tmax)));
    BOOST_OUTCOME_TRY(auto&& initial_conditions, init_condition_germany(t0_day, end_date, data_dir));

    // std::exit(1);

    // values obtained from function init_condition_germany:
    // auto initial_conditions = std::vector<double>{
    //     3.12341e+06, 7.70638,  1074.96,   0.00310779,  950.896, 589.014,     0.00173645, 521.748, 0,
    //     0,           0,        3075.29,   0.00303951,  129.243, 0,           0,          0,       6.86252,
    //     8.66913e-06, 0.122063, 7.32792,   1.01578e-05, 0.30537, 869520,      0,          0,       0,
    //     1.50414e+06, 21082.5,  1179.33,   21.4345,     16331.7, 658.904,     11.8292,    9091.16, 0,
    //     0,           0,        6340.84,   37.531,      4081.3,  0,           0,          0,       16.9719,
    //     0.127313,    4.68545,  3.67226,   0.035589,    2.3988,  6.36667e+06, 0,          0,       0,
    //     7386.97,     18982.1,  7.71391,   34.0561,     93414.6, 4.32421,     18.5192,    50639.4, 0,
    //     0,           0,        60.2175,   89.3921,     34213.9, 0,           0,          0,       0.220946,
    //     0.500617,    63.6993,  0.0359183, 0.0353868,   7.66378, 1.88185e+07, 0,          0,       0,
    //     15835.8,     23679.7,  12.4859,   35.2544,     115563,  7.70576,     19.1553,    62599.4, 0,
    //     0,           0,        105.556,   95.1216,     43452.9, 0,           0,          0,       1.2757,
    //     1.95773,     296.845,  0.277464,  0.144788,    37.1336, 2.81485e+07, 0,          0,       0,
    //     37977.7,     6679.25,  26.9146,   5.24158,     38205.9, 14.4734,     2.83312,    20641.6, 0,
    //     0,           0,        198.775,   14.6604,     14915.3, 0,           0,          0,       8.72773,
    //     0.888655,    299.078,  1.71563,   0.306414,    165.911, 1.74695e+07, 0,          0,       0,
    //     620573,      2490.5,   165.058,   0.825594,    6340.88, 91.2535,     0.45329,    3479.26, 0,
    //     0,           0,        1240.2,    2.3071,      2496.77, 0,           0,          0,       151.55,
    //     0.361321,    131.834,  110.4,     0.517923,    316.547, 6.76992e+06, 0,          0,       0};

    // transfer initial conditions to model
    for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(6); i++) {
        for (auto j = mio::Index<mio::osecirvvs::InfectionState>(0); j < mio::osecirvvs::InfectionState::Count; ++j) {
            auto flat_indx            = model.populations.get_flat_index({i, j});
            model.populations[{i, j}] = initial_conditions[flat_indx];
        }
    }

    model.apply_constraints();
    params_graph.add_node(0, model);

    // auto all_fact_variant      = std::vector<double>{1.0, 1.5, 2.0, 2.5, 3.0};
    // auto all_red_fact_exposed  = std::vector<double>{0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7};
    // auto all_t_immunity        = std::vector<double>{10, 20, 30, 40, 50, 60, 70, 80};
    // auto all_t_waning_immunity = std::vector<double>{10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130};

    // std::cout << "num runs = "
    //           << all_fact_variant.size() * all_red_fact_exposed.size() * all_t_immunity.size() *
    //                  all_t_waning_immunity.size()
    //           << std::endl;

    std::random_device rd;
    std::mt19937 gen(rd());

    // double cur_error_min = 4.67252e+7;

    // while (true) {

    // wähle ein zufälligen wert fac_variant aus dem Intervall
    std::uniform_real_distribution<> dis_rho(3.5, 3.5);
    auto fac_variant = dis_rho(gen);

    std::uniform_real_distribution<> dis_red_fact_exposed(0.4, 0.4);
    auto red_fact_exposed = dis_red_fact_exposed(gen);

    std::uniform_real_distribution<> dis_t_immunity(1, 1);
    auto t_immunity = dis_t_immunity(gen);

    std::uniform_real_distribution<> dis_t_waning_immunity(365, 365);
    auto t_waning_immunity = dis_t_waning_immunity(gen);

    auto result_dir_tmp = mio::path_join(
        result_dir.string(), "fac_variant_" + std::to_string(fac_variant) + "_red_fact_exposed_" +
                                 std::to_string(red_fact_exposed) + "_t_immunity_" + std::to_string(t_immunity) +
                                 "_t_waning_immunity_" + std::to_string(t_waning_immunity));
    BOOST_OUTCOME_TRY(set_covid_parameters(params_graph.nodes()[0].property.parameters, fac_variant, red_fact_exposed));
    auto parameter_study = mio::ParameterStudy<mio::osecirvvs::Simulation<>>{params_graph, t0, tmax, dt, num_runs};

    auto ensemble = parameter_study.run(
        [&](auto&& graph) {
            return draw_sample(graph, false);
        },
        [&](auto results_graph, auto&& run_idx) {
            auto interpolated_result = mio::interpolate_simulation_result(results_graph);
            auto params              = std::vector<mio::osecirvvs::Model<double>>();
            params.reserve(results_graph.nodes().size());
            std::transform(results_graph.nodes().begin(), results_graph.nodes().end(), std::back_inserter(params),
                           [](auto&& node) {
                               return node.property.get_simulation().get_model();
                           });

            auto num_symptomatic_infected_per_day = std::vector<double>(tmax, 0.0);
            // auto node_indx                        = 0;
            // for (auto&& node : results_graph.nodes()) {
            //     auto model_node = node.property.get_simulation().get_model();
            //     for (auto t = 0; t < tmax; t++) {
            //         auto num_symptomatic_infected = 0.0;
            //         // iterate over all age groups
            //         for (auto i = mio::AgeGroup(0); i < mio::AgeGroup(6); i++) {
            //             num_symptomatic_infected +=
            //                 interpolated_result[node_indx].get_value(t)[model_node.populations.get_flat_index(
            //                     {i, mio::osecirvvs::InfectionState::InfectedSymptomsNaive})] +
            //                 interpolated_result[node_indx].get_value(t)[model_node.populations.get_flat_index(
            //                     {i, mio::osecirvvs::InfectionState::InfectedSymptomsPartialImmunity})] +
            //                 interpolated_result[node_indx].get_value(t)[model_node.populations.get_flat_index(
            //                     {i, mio::osecirvvs::InfectionState::InfectedSymptomsImprovedImmunity})];
            //         }
            //         num_symptomatic_infected_per_day[t] += num_symptomatic_infected;
            //     }
            //     node_indx++;
            // }

            // // divide every entry from amelag by the first entry and multiply with first entry of num_symptomatic_infected
            // for (size_t i = 0; i < amelag_data.size(); i++) {
            //     amelag_data[i] = amelag_data[i] / amelag_data[0] * num_symptomatic_infected_per_day[0];
            // }

            // // calculate the error
            // auto error = 0.0;
            // for (size_t i = 0; i < amelag_data.size(); i++) {
            //     error += std::pow(amelag_data[i] - num_symptomatic_infected_per_day[i], 2);
            // }
            // error = std::sqrt(error);

            // if (error < cur_error_min) {
            //     cur_error_min = error;
            //     std::cout << "new min error = " << cur_error_min << std::endl;
            //     result_dir_tmp = result_dir_tmp + "error_" + std::to_string(cur_error_min);
            //     boost::filesystem::path dir(result_dir_tmp);
            //     bool created = boost::filesystem::create_directories(dir);
            //     auto save_status_run =
            //         save_result_with_params(interpolated_result, params, {0}, result_dir_tmp, run_idx);
            //     mio::unused(save_status_run, created);
            // }
            mio::unused(run_idx);
            return std::make_pair(interpolated_result, params);
        });

    // if (cur_error_min < curr_error_copy) {
    // curr_error_copy = cur_error_min;
    // result_dir_tmp  = result_dir_tmp + "error_" + std::to_string(cur_error_min);
    boost::filesystem::path dir(result_dir_tmp);
    bool created = boost::filesystem::create_directories(dir);
    mio::unused(created);

    auto ensemble_results = std::vector<std::vector<mio::TimeSeries<double>>>{};
    ensemble_results.reserve(ensemble.size());
    auto ensemble_params = std::vector<std::vector<mio::osecirvvs::Model<double>>>{};
    ensemble_params.reserve(ensemble.size());
    for (auto&& run : ensemble) {
        ensemble_results.emplace_back(std::move(run.first));
        ensemble_params.emplace_back(std::move(run.second));
    }
    BOOST_OUTCOME_TRY(save_results(ensemble_results, ensemble_params, {0}, result_dir_tmp, false));

    std::cout << "saved results to " << result_dir_tmp << std::endl;
    // }
    // }
    return mio::success();
}

int main()
{
    mio::mpi::init();
    const std::string memilio_dir = "/localdata1/pure/memilio";
    const std::string data_dir    = mio::path_join(memilio_dir, "data");
    std::string results_dir       = mio::path_join(memilio_dir, "results/amelag_single_vvs/");
    auto status_run               = run(data_dir, results_dir);
    mio::mpi::finalize();
    return 0;
}
