/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Martin Siggel, Daniel Abele, Martin J. Kuehn, Jan Kleinert, Khoa Nguyen
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

#include "pybind_util.h"
#include "utils/custom_index_array.h"
#include "utils/parameter_set.h"
#include "utils/index.h"
//#include "memilio/utils/history.h"
#include "abm/abm.h"
#include "abm/household.h"
#include "pybind11/attr.h"
#include "pybind11/cast.h"
#include "pybind11/pybind11.h"
#include "pybind11/operators.h"
#include <type_traits>

namespace py = pybind11;

//time point logger
struct LogTimePoint : mio::LogAlways {
    using Type = double;
    static Type log(const mio::abm::Simulation& sim)
    {
        return sim.get_time().hours();
    }
};

//LocationId logger
struct LogLocationIds : mio::LogOnce {
    using Type = std::vector<std::tuple<mio::abm::LocationType, uint32_t>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<std::tuple<mio::abm::LocationType, uint32_t>> location_ids{};
        for (auto&& location : sim.get_world().get_locations()) {
            location_ids.push_back(std::make_tuple(location.get_type(), location.get_index()));
        }
        return location_ids;
    }
};

//agent logger
struct LogPersonsPerLocationAndInfectionTime : mio::LogAlways {
    using Type = std::vector<std::tuple<mio::abm::LocationId, uint32_t, mio::abm::TimeSpan, mio::abm::InfectionState>>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<std::tuple<mio::abm::LocationId, uint32_t, mio::abm::TimeSpan, mio::abm::InfectionState>>
            location_ids_person{};
        for (auto&& person : sim.get_world().get_persons()) {
            location_ids_person.push_back(std::make_tuple(person.get_location().get_id(), person.get_person_id(),
                                                          person.get_time_since_transmission(),
                                                          person.get_infection_state(sim.get_time())));
        }
        return location_ids_person;
    }
};

PYBIND11_MODULE(_simulation_abm, m)
{
    pymio::iterable_enum<mio::abm::InfectionState>(m, "InfectionState", py::module_local{})
        .value("Susceptible", mio::abm::InfectionState::Susceptible)
        .value("Exposed", mio::abm::InfectionState::Exposed)
        .value("InfectedNoSymptoms", mio::abm::InfectionState::InfectedNoSymptoms)
        .value("InfectedSymptoms", mio::abm::InfectionState::InfectedSymptoms)
        .value("InfectedSevere", mio::abm::InfectionState::InfectedSevere)
        .value("InfectedCritical", mio::abm::InfectionState::InfectedCritical)
        .value("Recovered", mio::abm::InfectionState::Recovered)
        .value("Dead", mio::abm::InfectionState::Dead)
        .value("Count", mio::abm::InfectionState::Count);

    pymio::iterable_enum<mio::abm::AgeGroup>(m, "AgeGroup")
        .value("Age0to4", mio::abm::AgeGroup::Age0to4)
        .value("Age5to14", mio::abm::AgeGroup::Age5to14)
        .value("Age15to34", mio::abm::AgeGroup::Age15to34)
        .value("Age35to59", mio::abm::AgeGroup::Age35to59)
        .value("Age60to79", mio::abm::AgeGroup::Age60to79)
        .value("Age80plus", mio::abm::AgeGroup::Age80plus);

    pymio::iterable_enum<mio::abm::VirusVariant>(m, "VirusVariant").value("Wildtype", mio::abm::VirusVariant::Wildtype);

    pymio::iterable_enum<mio::abm::VaccinationState>(m, "VaccinationState")
        .value("Unvaccinated", mio::abm::VaccinationState::Unvaccinated)
        .value("Vaccinated", mio::abm::VaccinationState::Vaccinated);

    pymio::iterable_enum<mio::abm::LocationType>(m, "LocationType")
        .value("Home", mio::abm::LocationType::Home)
        .value("School", mio::abm::LocationType::School)
        .value("Work", mio::abm::LocationType::Work)
        .value("SocialEvent", mio::abm::LocationType::SocialEvent)
        .value("BasicsShop", mio::abm::LocationType::BasicsShop)
        .value("Hospital", mio::abm::LocationType::Hospital)
        .value("ICU", mio::abm::LocationType::ICU)
        .value("Car", mio::abm::LocationType::Car)
        .value("PublicTransport", mio::abm::LocationType::PublicTransport)
        .value("TransportWithoutContact", mio::abm::LocationType::TransportWithoutContact)
        .value("Cemetery", mio::abm::LocationType::Cemetery);

    py::class_<mio::abm::TestParameters>(m, "TestParameters")
        .def(py::init<double, double>())
        .def_readwrite("sensitivity", &mio::abm::TestParameters::sensitivity)
        .def_readwrite("specificity", &mio::abm::TestParameters::specificity);

    pymio::bind_Index<mio::abm::AgeGroup>(m, "AgeIndex");
    pymio::bind_Index<mio::abm::VaccinationState>(m, "VaccinationIndex");
    pymio::bind_CustomIndexArray<mio::UncertainValue, mio::abm::VirusVariant, mio::abm::AgeGroup,
                                 mio::abm::VaccinationState>(m, "_AgeVaccinationParameterArray");
    pymio::bind_ParameterSet<mio::abm::GlobalInfectionParameters>(m, "GlobalInfectionParameters").def(py::init<>());
    pymio::bind_ParameterSet<mio::abm::LocalInfectionParameters>(m, "LocalInfectionParameters").def(py::init<>());
    pymio::bind_ParameterSet<mio::abm::MigrationParameters>(m, "MigrationParameters").def(py::init<>());

    py::class_<mio::abm::TimeSpan>(m, "TimeSpan")
        .def(py::init<int>(), py::arg("seconds") = 0)
        .def_property_readonly("seconds", &mio::abm::TimeSpan::seconds)
        .def_property_readonly("hours", &mio::abm::TimeSpan::hours)
        .def_property_readonly("days", &mio::abm::TimeSpan::days)
        .def(py::self + py::self)
        .def(py::self += py::self)
        .def(py::self - py::self)
        .def(py::self -= py::self)
        .def(py::self * int{})
        .def(py::self *= int{})
        .def(py::self / int{})
        .def(py::self /= int{})
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self < py::self)
        .def(py::self <= py::self)
        .def(py::self > py::self)
        .def(py::self <= py::self);

    m.def("seconds", &mio::abm::seconds);
    m.def("minutes", &mio::abm::minutes);
    m.def("hours", &mio::abm::hours);
    m.def("days", py::overload_cast<int>(&mio::abm::days));

    py::class_<mio::abm::TimePoint>(m, "TimePoint")
        .def(py::init<int>(), py::arg("seconds") = 0)
        .def_property_readonly("seconds", &mio::abm::TimePoint::seconds)
        .def_property_readonly("days", &mio::abm::TimePoint::days)
        .def_property_readonly("hours", &mio::abm::TimePoint::hours)
        .def_property_readonly("day_of_week", &mio::abm::TimePoint::day_of_week)
        .def_property_readonly("hour_of_day", &mio::abm::TimePoint::hour_of_day)
        .def_property_readonly("time_since_midnight", &mio::abm::TimePoint::time_since_midnight)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self < py::self)
        .def(py::self <= py::self)
        .def(py::self > py::self)
        .def(py::self >= py::self)
        .def(py::self - py::self)
        .def(py::self + mio::abm::TimeSpan{})
        .def(py::self += mio::abm::TimeSpan{})
        .def(py::self - mio::abm::TimeSpan{})
        .def(py::self -= mio::abm::TimeSpan{});

    py::class_<mio::abm::LocationId>(m, "LocationId")
        .def(py::init([](uint32_t idx, mio::abm::LocationType type) {
            return mio::abm::LocationId{idx, type};
        }))
        .def_readwrite("index", &mio::abm::LocationId::index)
        .def_readwrite("type", &mio::abm::LocationId::type)
        .def(py::self == py::self)
        .def(py::self != py::self);

    py::class_<mio::abm::Person>(m, "Person")
        .def("set_assigned_location", py::overload_cast<mio::abm::LocationId>(&mio::abm::Person::set_assigned_location))
        .def("add_new_infection",
             [](mio::abm::Person& self, mio::abm::Infection& infection, mio::abm::TimePoint t) {
                 self.add_new_infection(std::move(infection), t);
             })
        .def_property_readonly("location", py::overload_cast<>(&mio::abm::Person::get_location, py::const_))
        .def_property_readonly("age", &mio::abm::Person::get_age)
        .def_property_readonly("is_in_quarantine", &mio::abm::Person::is_in_quarantine);

    py::class_<mio::abm::HouseholdMember>(m, "HouseholdMember")
        .def(py::init<>())
        .def("set_age_weight", &mio::abm::HouseholdMember::set_age_weight);

    py::class_<mio::abm::Household>(m, "Household")
        .def(py::init<>())
        .def("add_members", &mio::abm::Household::add_members);

    m.def("add_household_group_to_world", &mio::abm::add_household_group_to_world);

    py::class_<mio::abm::HouseholdGroup>(m, "HouseholdGroup")
        .def(py::init<>())
        .def("add_households", &mio::abm::HouseholdGroup::add_households);

    py::class_<mio::abm::TestingCriteria>(m, "TestingCriteria")
        .def(py::init<const std::vector<mio::abm::AgeGroup>&, const std::vector<mio::abm::LocationType>&,
                      const std::vector<mio::abm::InfectionState>&>(),
             py::arg("age_groups"), py::arg("location_types"), py::arg("infection_states"));

    py::class_<mio::abm::GenericTest>(m, "GenericTest").def(py::init<>());
    py::class_<mio::abm::AntigenTest, mio::abm::GenericTest>(m, "AntigenTest").def(py::init<>());
    py::class_<mio::abm::PCRTest, mio::abm::GenericTest>(m, "PCRTest").def(py::init<>());

    py::class_<mio::abm::TestingScheme>(m, "TestingScheme")
        .def(py::init<const std::vector<mio::abm::TestingCriteria>&, mio::abm::TimeSpan, mio::abm::TimePoint,
                      mio::abm::TimePoint, const mio::abm::GenericTest&, double>(),
             py::arg("testing_criteria"), py::arg("testing_min_time_since_last_test"), py::arg("start_date"),
             py::arg("end_date"), py::arg("test_type"), py::arg("probability"))
        .def_property_readonly("active", &mio::abm::TestingScheme::is_active);

    py::class_<mio::abm::TestingStrategy>(m, "TestingStrategy")
        .def(py::init<const std::vector<mio::abm::TestingScheme>&>());

    py::class_<mio::abm::Infection>(m, "Infection")
        .def(py::init<mio::abm::VirusVariant, mio::abm::AgeGroup, const mio::abm::GlobalInfectionParameters&,
                      mio::abm::TimePoint, mio::abm::InfectionState, bool>());

    py::class_<mio::abm::Location>(m, "Location")
        .def("set_capacity", &mio::abm::Location::set_capacity)
        .def_property_readonly("type", &mio::abm::Location::get_type)
        .def_property_readonly("index", &mio::abm::Location::get_index)
        .def_property_readonly("population", &mio::abm::Location::get_subpopulations,
                               py::return_value_policy::reference_internal)
        .def_property("infection_parameters",
                      py::overload_cast<>(&mio::abm::Location::get_infection_parameters, py::const_),
                      [](mio::abm::Location& self, mio::abm::LocalInfectionParameters params) {
                          self.get_infection_parameters() = params;
                      });

    pymio::bind_Range<decltype(std::declval<mio::abm::World>().get_locations())>(m, "_WorldLocationsRange");
    pymio::bind_Range<decltype(std::declval<mio::abm::World>().get_persons())>(m, "_WorldPersonsRange");

    py::class_<mio::abm::Trip>(m, "Trip")
        .def(py::init<uint32_t, mio::abm::TimePoint, mio::abm::LocationId, mio::abm::LocationId,
                      std::vector<uint32_t>>(),
             py::arg("person_id"), py::arg("time"), py::arg("destination"), py::arg("origin"),
             py::arg("cells") = std::vector<uint32_t>())
        .def_readwrite("person_id", &mio::abm::Trip::person_id)
        .def_readwrite("time", &mio::abm::Trip::time)
        .def_readwrite("destination", &mio::abm::Trip::migration_destination)
        .def_readwrite("origin", &mio::abm::Trip::migration_origin)
        .def_readwrite("cells", &mio::abm::Trip::cells);

    py::class_<mio::abm::TripList>(m, "TripList")
        .def(py::init<>())
        .def("add_trip", &mio::abm::TripList::add_trip)
        .def_property_readonly("next_trip", &mio::abm::TripList::get_next_trip)
        .def_property_readonly("num_trips", &mio::abm::TripList::num_trips);

    py::class_<mio::abm::World>(m, "World")
        .def(py::init<mio::abm::GlobalInfectionParameters>(),
             py::arg("infection_parameters") = mio::abm::GlobalInfectionParameters{})
        .def("add_location", &mio::abm::World::add_location, py::arg("location_type"), py::arg("num_cells") = 1)
        .def("add_person", &mio::abm::World::add_person, py::arg("location_id"), py::arg("age_group"),
             py::return_value_policy::reference_internal)
        .def("get_individualized_location",
             py::overload_cast<mio::abm::LocationId>(&mio::abm::World::get_individualized_location, py::const_),
             py::return_value_policy::reference_internal)
        .def_property_readonly("locations", &mio::abm::World::get_locations,
                               py::keep_alive<1, 0>{}) //keep this world alive while contents are referenced in ranges
        .def_property_readonly("persons", &mio::abm::World::get_persons, py::keep_alive<1, 0>{})
        .def_property(
            "trip_list", py::overload_cast<>(&mio::abm::World::get_trip_list),
            [](mio::abm::World& self, const mio::abm::TripList& list) {
                self.get_trip_list() = list;
            },
            py::return_value_policy::reference_internal)
        .def_property("use_migration_rules", py::overload_cast<>(&mio::abm::World::use_migration_rules, py::const_),
                      py::overload_cast<bool>(&mio::abm::World::use_migration_rules))
        .def_property(
            "infection_parameters", py::overload_cast<>(&mio::abm::World::get_global_infection_parameters, py::const_),
            [](mio::abm::World& self, mio::abm::GlobalInfectionParameters params) {
                self.get_global_infection_parameters() = params;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "migration_parameters", py::overload_cast<>(&mio::abm::World::get_migration_parameters, py::const_),
            [](mio::abm::World& self, mio::abm::MigrationParameters params) {
                self.get_migration_parameters() = params;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "testing_strategy", py::overload_cast<>(&mio::abm::World::get_testing_strategy, py::const_),
            [](mio::abm::World& self, mio::abm::TestingStrategy strategy) {
                self.get_testing_strategy() = strategy;
            },
            py::return_value_policy::reference_internal);

    m.def(
        "set_viral_load_parameters",
        [](mio::abm::GlobalInfectionParameters& infection_params, mio::abm::VirusVariant variant,
           mio::abm::AgeGroup age, mio::abm::VaccinationState state, double min_peak, double max_peak,
           double min_incline, double max_incline, double min_decline, double max_decline) {
            infection_params.get<mio::abm::ViralLoadDistributions>()[{variant, age, state}] =
                mio::abm::ViralLoadDistributionsParameters{
                    {min_peak, max_peak}, {min_incline, max_incline}, {min_decline, max_decline}};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_infectivity_parameters",
        [](mio::abm::GlobalInfectionParameters& infection_params, mio::abm::VirusVariant variant,
           mio::abm::AgeGroup age, double min_alpha, double max_alpha, double min_beta, double max_beta) {
            infection_params.get<mio::abm::InfectivityDistributions>()[{variant, age}] =
                mio::abm::InfectivityDistributionsParameters{{min_alpha, max_alpha}, {min_beta, max_beta}};
        },
        py::return_value_policy::reference_internal);

    // m.def("get_instance", [](mio::abm::GlobalInfectionParameters& params) {
    //     auto inf_params = params.get<mio::abm::InfectivityDistributions>()[{mio::abm::VirusVariant::Wildtype,
    //                                                                         mio::abm::AgeGroup::Age80plus}];
    //     return inf_params.infectivity_alpha.get_distribution_instance()(inf_params.infectivity_alpha.params);
    // });

    py::class_<mio::abm::Simulation>(m, "Simulation")
        .def(py::init<mio::abm::TimePoint>())
        .def("advance",
             &mio::abm::Simulation::advance<mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                                         LogPersonsPerLocationAndInfectionTime>>)
        //.def("advance", &mio::abm::Simulation::advance)
        .def_property_readonly("result", &mio::abm::Simulation::get_result)
        .def_property_readonly("world", py::overload_cast<>(&mio::abm::Simulation::get_world),
                               py::return_value_policy::reference_internal);

    py::class_<
        mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds, LogPersonsPerLocationAndInfectionTime>>(
        m, "History")
        .def(py::init<>())
        .def_property_readonly("log", [](mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                                      LogPersonsPerLocationAndInfectionTime>& self) {
            return self.get_log();
        });
}
