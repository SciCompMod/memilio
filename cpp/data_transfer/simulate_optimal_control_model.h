/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Maximilian Betz
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
* 
* The documentation of the Ipopt::TNLP member functions  in Secirvvs_NLP
* is extracted from the Ipopt documentation
*/

#include "ad/ad.hpp"

#include "ode_secirvvs/model.h"
#include "ode_secirvvs/parameter_space.h"
#include "ode_secirvvs/parameters_io.h"
#include "memilio/compartments/simulation.h"
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/epidemiology/populations.h"
#include "memilio/io/mobility_io.h"
#include "memilio/io/result_io.h"
#include "memilio/io/cli.h"
#include "memilio/math/eigen.h"
#include "memilio/data/analyze_result.h"
#include "memilio/utils/logging.h"
#include "memilio/utils/time_series.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/base_dir.h"
#include "memilio/utils/stl_util.h"

#include "boost/outcome/try.hpp"
#include "boost/outcome/result.hpp"
#include "boost/optional.hpp"
#include "boost/filesystem.hpp"

#include "IpTNLP.hpp"
#include "IpIpoptApplication.hpp"
#include <fstream>
#include <type_traits>

using internal_type = double;
using gt1s_type = ad::gt1s<internal_type>::type;

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

/**
 * different types of NPI, used as DampingType.
 */
enum class Intervention
{
    Home,
    SchoolClosure,
    HomeOffice,
    GatheringBanFacilitiesClosure,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Count,
};

/**
 * different level of NPI, used as DampingLevel.
 */
enum class InterventionLevel
{
    Main,
    PhysicalDistanceAndMasks,
    SeniorAwareness,
    Holidays,
    Count,
};

static const std::map<ContactLocation, std::string> contact_locations = {{ContactLocation::Home, "home"},
                                                                         {ContactLocation::School, "school_pf_eig"},
                                                                         {ContactLocation::Work, "work"},
                                                                         {ContactLocation::Other, "other"}
                                                                        };

template <typename FP, class Tag>
void set_parameters(mio::osecirvvs::Parameters<FP>& params, const std::vector<FP>& parameter_values)
{
    std::copy(parameter_values.begin(), parameter_values.end(), params.template get<Tag>());
}

template <typename FP, class Tag>
void set_parameters(mio::osecirvvs::Parameters<FP>& params, double parameter_values)
{
    std::fill(params.template get<Tag>().begin(), params.template get<Tag>().end(), FP(parameter_values));
}

template <typename FP>
using EBModel = mio::osecirvvs::Model<FP>;

template <typename FP>
mio::IOResult<void> set_covid_parameters(mio::osecirvvs::Parameters<FP>& params, std::string temp_dir, mio::Date start_date)
{
    BOOST_OUTCOME_TRY(auto&& parameter_list, mio::read_json(mio::path_join(temp_dir, "parameter_list.json")));

    std::map<std::string, std::string> id_to_name{};
    for(auto entry:  parameter_list)
    {
        id_to_name[entry["id"].asString()] = entry["name"].asString();
    }

    // only uses group 0 for all parameters
    BOOST_OUTCOME_TRY(auto&& scenario_data_run, mio::read_json(mio::path_join(temp_dir, "scenario_data_run.json")));
    std::map<std::string, double> parameter_values{};
    for(auto parameter:  scenario_data_run["modelParameters"])
    {
        std::string parameter_name = id_to_name[parameter["parameterId"].asString()];
        parameter_values[parameter_name] = 0.5 * (parameter["values"][0]["valueMin"].asDouble() + 
                                                    parameter["values"][0]["valueMax"].asDouble());
    }

    //times
    set_parameters<FP, mio::osecirvvs::TimeExposed<FP>>(params, parameter_values["TimeExposed"]);
    set_parameters<FP, mio::osecirvvs::TimeInfectedNoSymptoms<FP>>(params, parameter_values["TimeInfectedNoSymptoms"]);
    set_parameters<FP, mio::osecirvvs::TimeInfectedSymptoms<FP>>(params, parameter_values["TimeInfectedSymptoms"]);
    set_parameters<FP, mio::osecirvvs::TimeInfectedSevere<FP>>(params, parameter_values["TimeInfectedSevere"]);
    set_parameters<FP, mio::osecirvvs::TimeInfectedCritical<FP>>(params, parameter_values["TimeInfectedCritical"]);

    //probabilities
    // set_parameters<FP, mio::osecirvvs::TransmissionProbabilityOnContact<FP>>(params, 0.2); // use this value if data doesnt conlude to solution
    set_parameters<FP, mio::osecirvvs::TransmissionProbabilityOnContact<FP>>(params, parameter_values["TransmissionProbabilityOnContact"]);
    set_parameters<FP, mio::osecirvvs::RelativeTransmissionNoSymptoms<FP>>(params, parameter_values["RelativeTransmissionNoSymptoms"]);
    // The precise value between Risk* (situation under control) and MaxRisk* (situation not under control)
    // depends on incidence and test and trace capacity
    set_parameters<FP, mio::osecirvvs::RiskOfInfectionFromSymptomatic<FP>>(params, parameter_values["RiskOfInfectionFromSymptomatic"]);
    set_parameters<FP, mio::osecirvvs::MaxRiskOfInfectionFromSymptomatic<FP>>(params, parameter_values["MaxRiskOfInfectionFromSymptomatic"]);
    set_parameters<FP, mio::osecirvvs::RecoveredPerInfectedNoSymptoms<FP>>(params, parameter_values["RecoveredPerInfectedNoSymptoms"]);
    set_parameters<FP, mio::osecirvvs::SeverePerInfectedSymptoms<FP>>(params, parameter_values["SeverePerInfectedSymptoms"]);
    set_parameters<FP, mio::osecirvvs::CriticalPerSevere<FP>>(params, parameter_values["CriticalPerSevere"]);
    set_parameters<FP, mio::osecirvvs::DeathsPerCritical<FP>>(params, parameter_values["DeathsPerCritical"]);

    set_parameters<FP, mio::osecirvvs::ReducExposedPartialImmunity<FP>>(params, parameter_values["ReducedExposedPartialImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducExposedImprovedImmunity<FP>>(params, parameter_values["ReducedExposedImprovedImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducInfectedSymptomsPartialImmunity<FP>>(params, parameter_values["ReducedInfectedSymptomsPartialImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducInfectedSymptomsImprovedImmunity<FP>>(params, parameter_values["ReducedInfectedSymptomsImprovedImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducInfectedSevereCriticalDeadPartialImmunity<FP>>(params, parameter_values["ReducedInfectedSevereCriticalDeadPartialImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducInfectedSevereCriticalDeadImprovedImmunity<FP>>(params, parameter_values["ReducedInfectedSevereCriticalDeadImprovedImmunity"]);
    set_parameters<FP, mio::osecirvvs::ReducTimeInfectedMild<FP>>(params, parameter_values["ReducedTimeInfectedMild"]);

    params.template get<mio::osecirvvs::TestAndTraceCapacity<FP>>() = 100;
    params.template get<mio::osecirvvs::StartDay<FP>>() = get_day_in_year(start_date);
    params.template get<mio::osecirvvs::Seasonality<FP>>() = FP(parameter_values["Seasonality"]);

    return mio::success();
}

template <typename FP>
mio::IOResult<void> set_contact_matrices(mio::osecirvvs::Parameters<FP>& params, std::string data_dir)
{
    auto num_groups = size_t(params.get_num_groups());
    auto contact_matrices = mio::ContactMatrixGroup<FP>(contact_locations.size(), num_groups);
    for (auto&& contact_location : contact_locations) {
        BOOST_OUTCOME_TRY(auto&& baseline,
                          mio::read_mobility_plain(mio::path_join(data_dir, "Germany", "contacts", ("baseline_" + contact_location.second + ".txt"))));

        contact_matrices[size_t(contact_location.first)].get_baseline() = baseline.cast<FP>();
        contact_matrices[size_t(contact_location.first)].get_minimum()  = Eigen::Matrix<FP, Eigen::Dynamic, Eigen::Dynamic>::Zero(num_groups, num_groups);
    }
    params.template get<mio::osecirvvs::ContactPatterns<FP>>() = mio::UncertainContactMatrix<FP>(contact_matrices);

    return mio::success();
}

mio::IOResult<void> set_population_data(EBModel<internal_type>& model, std::string data_dir, double tmax, mio::Date start_date)
{  

    auto scaling_factor_infected = std::vector<double>(size_t(model.parameters.get_num_groups()), 2.5);
    auto scaling_factor_icu      = 1.0;
    auto tnt_capacity_factor     = 7.5 / 100000.;


    std::vector<EBModel<internal_type>> nodes(1, model);

    // std::vector<int> node_ids = {0};
    std::string pydata_dir = mio::path_join(data_dir, "Germany", "pydata");
    BOOST_OUTCOME_TRY(mio::osecirvvs::read_input_data_germany<EBModel<internal_type>>(nodes, start_date, scaling_factor_infected, scaling_factor_icu, pydata_dir, tmax));

    model.populations = nodes[0].populations;
    mio::unused(tnt_capacity_factor);

    return mio::success();
}

mio::IOResult<void> set_initial_values(EBModel<internal_type>& model_double, EBModel<gt1s_type>& model_ad_gt1s, 
                                       std::string data_dir, std::string temp_dir, double tmax)
{
    auto scenario_data_run = mio::read_json(mio::path_join(temp_dir, "scenario_data_run.json")).value();
    mio::Date start_date = mio::parse_date(scenario_data_run["startDate"].asString()).value();

    auto num_groups = size_t(model_double.parameters.get_num_groups());
    // populate model_double with data
    mio::osecirvvs::Parameters<internal_type> params_double(num_groups);

    BOOST_OUTCOME_TRY(set_covid_parameters<internal_type>(params_double, temp_dir, start_date));
    BOOST_OUTCOME_TRY(set_contact_matrices<internal_type>(params_double, data_dir));
    model_double.parameters = params_double;

    BOOST_OUTCOME_TRY(set_population_data(model_double, data_dir, tmax, start_date));

    // populate model_ad_gt1s with data
    mio::osecirvvs::Parameters<gt1s_type> params_ad_gt1s(num_groups);

    BOOST_OUTCOME_TRY(set_covid_parameters<gt1s_type>(params_ad_gt1s, temp_dir, start_date));
    BOOST_OUTCOME_TRY(set_contact_matrices<gt1s_type>(params_ad_gt1s, data_dir));
    model_ad_gt1s.parameters = params_ad_gt1s;
    
    model_ad_gt1s.populations = model_double.populations.template convert<gt1s_type>();

    return mio::success();
}

template <typename FP>
mio::IOResult<void> set_npis(mio::osecirvvs::Parameters<FP>& params, FP t0, FP tmax, const std::vector<FP>& x, const int numControlIntervals)
{
    // mio::unused(t0, x, numControlIntervals);
    auto damping_helper = [=](mio::SimulationTime<FP> t, FP min, FP max, mio::DampingLevel damping_level, mio::DampingType damping_type, const std::vector<size_t> location,
                Eigen::VectorX<FP> group_weights) {
        auto p = mio::UncertainValue<FP>(0.5 * (max + min));
        p.set_distribution(mio::ParameterDistributionUniform(ad::value(min), ad::value(max)));
        return mio::DampingSampling<FP>(p, damping_level, damping_type, t, location, group_weights);
    };

    auto group_weights_all = Eigen::VectorX<FP>::Constant(size_t(params.get_num_groups()), 1.0);

    auto school_closure = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(int(InterventionLevel::Main)),
                                            mio::DampingType(int(Intervention::SchoolClosure)),
                                            {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto home_office = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(int(InterventionLevel::Main)),
                                            mio::DampingType(int(Intervention::HomeOffice)),
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto physical_distancing_school = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::School)}, group_weights_all);
    };
    auto physical_distancing_work = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::Work)}, group_weights_all);
    };
    auto physical_distancing_other = [=](mio::SimulationTime<FP> t, FP v) {
        return damping_helper(t, v, v, mio::DampingLevel(int(InterventionLevel::PhysicalDistanceAndMasks)),
                                            mio::DampingType(int(Intervention::PhysicalDistanceAndMasks)),
                                            {size_t(ContactLocation::Other)}, group_weights_all);
    };

    auto& contacts         = params.template get<mio::osecirvvs::ContactPatterns<FP>>();
    auto& contact_dampings = contacts.get_dampings();

    const int step_size = ad::value(tmax - t0) / numControlIntervals;
    int t;
    for (int controlIndex = 0; controlIndex < numControlIntervals; ++controlIndex)
    {
        t = ad::value(t0) + controlIndex * step_size;
        contact_dampings.push_back(school_closure(mio::SimulationTime<FP>(t), (1. * x.at(controlIndex))));
        contact_dampings.push_back(home_office(mio::SimulationTime<FP>(t), (0.25 * x.at(controlIndex + numControlIntervals))));
        contact_dampings.push_back(physical_distancing_school(mio::SimulationTime<FP>(t), (0.25 * x.at(controlIndex + 2 * numControlIntervals))));
        contact_dampings.push_back(physical_distancing_work(mio::SimulationTime<FP>(t), (0.25 * x.at(controlIndex + 3 * numControlIntervals))));
        contact_dampings.push_back(physical_distancing_other(mio::SimulationTime<FP>(t), (0.35 * x.at(controlIndex + 4 * numControlIntervals))));
    }
    contacts.make_matrix();

    return mio::success();
}

/**
 * @brief set control values of the pandemic ODE
 * @param model an instance of the pandemic model
 */
template <typename FP>
mio::IOResult<void> set_control_values(EBModel<FP>& model, FP t0, FP tmax, const std::vector<FP>& x, const int numControlIntervals)
{
    auto& params = model.parameters;
    BOOST_OUTCOME_TRY(set_npis<FP>(params, t0, tmax, x, numControlIntervals));

    return mio::success();
}
