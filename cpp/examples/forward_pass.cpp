/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Nils Wassmuth
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
#include "abm/household.h"
#include "abm/lockdown_rules.h"
#include "abm/model.h"
#include "abm/common_abm_loggers.h"
#include "forward_pass.h"

namespace
{
const size_t NUM_AGE_GROUPS      = 4;
const auto AGE_GROUP_0_TO_4      = mio::AgeGroup(0);
const auto AGE_GROUP_5_TO_14     = mio::AgeGroup(1);
const auto AGE_GROUP_15_TO_34    = mio::AgeGroup(2);
const auto AGE_GROUP_35_TO_59    = mio::AgeGroup(3);
} // namespace

struct ABMPopulation::Impl {
    mio::abm::Model model;

    explicit Impl(int n_households)
        : model(NUM_AGE_GROUPS)
    {
        model.parameters.get<mio::abm::AgeGroupGotoSchool>() = false;
        model.parameters.get<mio::abm::AgeGroupGotoSchool>()[AGE_GROUP_5_TO_14] = true;
        model.parameters.get<mio::abm::AgeGroupGotoWork>().set_multiple(
            {AGE_GROUP_15_TO_34, AGE_GROUP_35_TO_59}, true);

        auto child = mio::abm::HouseholdMember(NUM_AGE_GROUPS);
        child.set_age_weight(AGE_GROUP_0_TO_4, 1);
        child.set_age_weight(AGE_GROUP_5_TO_14, 1);

        auto parent = mio::abm::HouseholdMember(NUM_AGE_GROUPS);
        parent.set_age_weight(AGE_GROUP_15_TO_34, 1);
        parent.set_age_weight(AGE_GROUP_35_TO_59, 1);

        auto twoPersonHousehold_group = mio::abm::HouseholdGroup();
        auto twoPersonHousehold_full  = mio::abm::Household();
        twoPersonHousehold_full.add_members(child, 1);
        twoPersonHousehold_full.add_members(parent, 1);
        twoPersonHousehold_group.add_households(twoPersonHousehold_full, n_households);
        add_household_group_to_model(model, twoPersonHousehold_group);

        auto threePersonHousehold_group = mio::abm::HouseholdGroup();
        auto threePersonHousehold_full  = mio::abm::Household();
        threePersonHousehold_full.add_members(child, 1);
        threePersonHousehold_full.add_members(parent, 2);
        threePersonHousehold_group.add_households(threePersonHousehold_full, n_households);
        add_household_group_to_model(model, threePersonHousehold_group);

        const int n_schools    = std::max(1, n_households / 200);
        const int n_workplaces = std::max(1, 3 * n_households / 200);
        const int n_shops      = std::max(1, 5 * n_households / 500);
        const int n_events     = std::max(1, 5 * n_households / 500);

        std::vector<mio::abm::LocationId> schools, workplaces, shops, events;
        for (int i = 0; i < n_schools; ++i) {
            auto loc = model.add_location(mio::abm::LocationType::School);
            model.get_location(loc).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
            schools.push_back(loc);
        }
        for (int i = 0; i < n_workplaces; ++i) {
            auto loc = model.add_location(mio::abm::LocationType::Work);
            model.get_location(loc).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
            model.get_location(loc).get_infection_parameters().get<mio::abm::ContactRates>().get_baseline()(
                AGE_GROUP_15_TO_34.get(), AGE_GROUP_15_TO_34.get()) = 10.0;
            workplaces.push_back(loc);
        }
        for (int i = 0; i < n_shops; ++i) {
            auto loc = model.add_location(mio::abm::LocationType::BasicsShop);
            model.get_location(loc).get_infection_parameters().set<mio::abm::MaximumContacts>(20);
            shops.push_back(loc);
        }
        for (int i = 0; i < n_events; ++i) {
            auto loc = model.add_location(mio::abm::LocationType::SocialEvent);
            model.get_location(loc).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
            events.push_back(loc);
        }
        auto hospital = model.add_location(mio::abm::LocationType::Hospital);
        model.get_location(hospital).get_infection_parameters().set<mio::abm::MaximumContacts>(5);
        auto icu = model.add_location(mio::abm::LocationType::ICU);
        model.get_location(icu).get_infection_parameters().set<mio::abm::MaximumContacts>(5);

        int school_ctr = 0, work_ctr = 0, shop_ctr = 0, event_ctr = 0;
        for (auto& person : model.get_persons()) {
            const auto id = person.get_id();
            model.assign_location(id, hospital);
            model.assign_location(id, icu);
            model.assign_location(id, shops[shop_ctr++ % n_shops]);
            model.assign_location(id, events[event_ctr++ % n_events]);
            if (person.get_age() == AGE_GROUP_5_TO_14) {
                model.assign_location(id, schools[school_ctr++ % n_schools]);
            }
            if (person.get_age() == AGE_GROUP_15_TO_34 || person.get_age() == AGE_GROUP_35_TO_59) {
                model.assign_location(id, workplaces[work_ctr++ % n_workplaces]);
            }
        }
    }
};

ABMPopulation::ABMPopulation(int n_households)
    : impl(std::make_shared<Impl>(n_households))
{
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> forward_pass(const ABMPopulation& population,
                                                          ScalarType beta, ScalarType kappa)
{
    mio::set_log_level(mio::LogLevel::warn);

    // Copy the pre-built population, then set inference parameters on the copy.
    auto model = population.impl->model;
    model.parameters.get<mio::abm::TimeExposedToNoSymptoms>() =
        mio::ParameterDistributionLogNormal(log(kappa), 1.);
    model.parameters.check_constraints();
    model.parameters.get<mio::abm::AerosolTransmissionRates>() = 1 / beta;

    // Randomise initial infection states on the copied model.
    auto start_date = mio::abm::TimePoint(0);
    std::vector<ScalarType> infection_distribution{0.9, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (auto& person : model.get_persons()) {
        mio::abm::InfectionState infection_state = mio::abm::InfectionState(
            mio::DiscreteDistribution<size_t>::get_instance()(mio::thread_local_rng(), infection_distribution));
        auto rng = mio::abm::PersonalRandomNumberGenerator(model.get_rng(), person);
        if (infection_state != mio::abm::InfectionState::Susceptible) {
            person.add_new_infection(mio::abm::Infection(rng, mio::abm::VirusVariant::Wildtype, person.get_age(),
                                                         model.parameters, start_date, infection_state));
        }
    }

    auto t0   = mio::abm::TimePoint(0);
    auto tmax = t0 + mio::abm::days(30);
    auto sim  = mio::abm::Simulation(t0, std::move(model));

    // clusters: 0 = school (5–14), 1 = work (15–34), 2 = work (35–59)
    mio::History<mio::DataWriterToMemory, mio::abm::LogCTCluster>
        historyTimeSeries{
            mio::abm::LogCTCluster{std::vector<mio::abm::LogCTCluster::ClusterSpec>{
                {AGE_GROUP_5_TO_14,  mio::abm::LocationType::School},
                {AGE_GROUP_15_TO_34, mio::abm::LocationType::Work},
                {AGE_GROUP_35_TO_59, mio::abm::LocationType::Work}}}};
    sim.advance(tmax, historyTimeSeries);

    const auto& hist_log = std::get<0>(historyTimeSeries.get_log());
    const int n_days     = static_cast<int>(hist_log.size());

    Eigen::MatrixXd school_result(n_days, 42); // [day, ct_0..ct_40]  school age-group
    Eigen::MatrixXd work_result(n_days, 83);   // [day, ct_0..ct_40 (15–34), ct_0..ct_40 (35–59)]

    for (int i = 0; i < n_days; ++i) {
        const auto& [time, hists] = hist_log[i];
        const double day = std::round(time.days());
        school_result(i, 0)             = day;
        school_result.row(i).tail(41)   = hists[0].cast<double>().transpose();
        work_result(i, 0)               = day;
        work_result.row(i).segment(1, 41)  = hists[1].cast<double>().transpose();
        work_result.row(i).segment(42, 41) = hists[2].cast<double>().transpose();
    }

    return {school_result, work_result};
}

std::pair<Eigen::MatrixXd, Eigen::MatrixXd> forward_pass(ScalarType beta, ScalarType kappa)
{
    return forward_pass(ABMPopulation{}, beta, kappa);
}
