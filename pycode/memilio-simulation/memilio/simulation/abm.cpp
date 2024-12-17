/* 
* Copyright (C) 2020-2024 MEmilio
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
#include "abm/vaccine.h"
#include "pybind11/attr.h"
#include "pybind11/cast.h"
#include "pybind11/pybind11.h"
#include "pybind11/operators.h"
#include "utils/logging.h"
#include <cstddef>
#include <cstdint>
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
    using Type = std::vector<mio::abm::LocationId>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<mio::abm::LocationId> location_ids{};
        for (auto&& location : sim.get_world().get_locations()) {
            location_ids.push_back(location.get_id());
        }
        return location_ids;
    }
};

//AgentId logger
struct LogAgentIds : mio::LogOnce {
    using Type = std::vector<uint32_t>;
    static Type log(const mio::abm::Simulation& sim)
    {
        std::vector<uint32_t> agent_ids{};
        for (auto&& person : sim.get_world().get_persons()) {
            agent_ids.push_back(person.get_person_id());
        }
        return agent_ids;
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

std::pair<double, double> get_my_and_sigma(double mean, double std)
{
    double my    = log(mean * mean / sqrt(mean * mean + std * std));
    double sigma = sqrt(log(1 + std * std / (mean * mean)));
    return {my, sigma};
}

PYBIND11_MODULE(_simulation_abm, m)
{
    //pymio::bind_logging(m, "LogLevel");

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

    pymio::iterable_enum<mio::abm::ExposureType>(m, "ExposureType")
        .value("NoProtection", mio::abm::ExposureType::NoProtection)
        .value("NaturalInfection", mio::abm::ExposureType::NaturalInfection)
        .value("GenericVaccine", mio::abm::ExposureType::GenericVaccine);

    pymio::iterable_enum<mio::abm::VirusVariant>(m, "VirusVariant").value("Wildtype", mio::abm::VirusVariant::Wildtype);

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

    pymio::bind_CustomIndexArray<mio::UncertainValue, mio::abm::VirusVariant, mio::AgeGroup>(m, "_AgeParameterArray");
    // pymio::bind_CustomIndexArray<LogNormalDistribution<double>::ParamType, mio::abm::VirusVariant, mio::AgeGroup>(
    //     m, "_AgeParameterDist");
    pymio::bind_Index<mio::abm::ExposureType>(m, "ExposureTypeIndex");
    pymio::bind_ParameterSet<mio::abm::ParametersBase>(m, "ParametersBase");
    py::class_<mio::abm::Parameters, mio::abm::ParametersBase>(m, "Parameters")
        .def(py::init<int>())
        .def("check_constraints", &mio::abm::Parameters::check_constraints);

    pymio::bind_ParameterSet<mio::abm::LocalInfectionParameters>(m, "LocalInfectionParameters").def(py::init<size_t>());

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
        .def("assigned_location",
             [](mio::abm::Person& self, mio::abm::LocationType type) {
                 return self.get_assigned_location_index(type);
             })
        .def("infection_state",
             [](mio::abm::Person& self, mio::abm::TimePoint t) {
                 return self.get_infection_state(t);
             })
        .def_property_readonly("infection", py::overload_cast<>(&mio::abm::Person::get_infection, py::const_))
        .def_property_readonly("location", py::overload_cast<>(&mio::abm::Person::get_location, py::const_))
        .def_property_readonly("age", &mio::abm::Person::get_age)
        .def_property_readonly("id", &mio::abm::Person::get_person_id)
        .def_property_readonly("is_in_quarantine", &mio::abm::Person::is_in_quarantine);

    py::class_<mio::abm::HouseholdMember>(m, "HouseholdMember")
        .def(py::init<int>())
        .def("set_age_weight", &mio::abm::HouseholdMember::set_age_weight);

    py::class_<mio::abm::Household>(m, "Household")
        .def(py::init<>())
        .def("add_members", &mio::abm::Household::add_members);

    m.def("add_household_group_to_world", &mio::abm::add_household_group_to_world);

    py::class_<mio::abm::HouseholdGroup>(m, "HouseholdGroup")
        .def(py::init<>())
        .def("add_households", &mio::abm::HouseholdGroup::add_households);

    py::class_<mio::abm::TestingCriteria>(m, "TestingCriteria")
        .def(py::init<const std::vector<mio::AgeGroup>&, const std::vector<mio::abm::InfectionState>&>(),
             py::arg("age_groups"), py::arg("infection_states"));

    py::class_<mio::abm::GenericTest>(m, "GenericTest").def(py::init<>());
    py::class_<mio::abm::AntigenTest, mio::abm::GenericTest>(m, "AntigenTest").def(py::init<>());
    py::class_<mio::abm::PCRTest, mio::abm::GenericTest>(m, "PCRTest").def(py::init<>());

    py::class_<mio::abm::TestingScheme>(m, "TestingScheme")
        .def(py::init<const mio::abm::TestingCriteria&, mio::abm::TimeSpan, mio::abm::TimePoint, mio::abm::TimePoint,
                      const mio::abm::GenericTest&, double>(),
             py::arg("testing_criteria"), py::arg("testing_min_time_since_last_test"), py::arg("start_date"),
             py::arg("end_date"), py::arg("test_type"), py::arg("probability"))
        .def_property_readonly("active", &mio::abm::TestingScheme::is_active);

    py::class_<mio::abm::Vaccination>(m, "Vaccination")
        .def(py::init<mio::abm::ExposureType, mio::abm::TimePoint>(), py::arg("exposure_type"), py::arg("time"))
        .def_readwrite("exposure_type", &mio::abm::Vaccination::exposure_type)
        .def_readwrite("time", &mio::abm::Vaccination::time);

    py::class_<mio::abm::TestingStrategy>(m, "TestingStrategy")
        .def(py::init<const std::unordered_map<mio::abm::LocationId, std::vector<mio::abm::TestingScheme>>&>());

    py::class_<mio::abm::Infection>(m, "Infection")
        .def(py::init([](mio::abm::World& world, mio::abm::Person& person, mio::abm::VirusVariant variant,
                         mio::abm::TimePoint start_date, mio::abm::InfectionState start_state, bool detected) {
            auto rng = mio::abm::Person::RandomNumberGenerator(world.get_rng(), person);
            return mio::abm::Infection(rng, variant, person.get_age(), world.parameters, start_date, start_state,
                                       person.get_latest_protection(), detected);
        }))
        .def("get_infection_start", &mio::abm::Infection::get_infection_start)
        .def("get_time_in_state", [](mio::abm::Infection& self, mio::abm::InfectionState state) {
            return self.get_time_in_state(state);
        });

    py::class_<mio::abm::Location>(m, "Location")
        .def("set_capacity", &mio::abm::Location::set_capacity)
        .def_property_readonly("type", &mio::abm::Location::get_type)
        .def_property_readonly("index", &mio::abm::Location::get_index)
        .def_property_readonly("population", &mio::abm::Location::get_subpopulation)
        .def_property("infection_parameters",
                      py::overload_cast<>(&mio::abm::Location::get_infection_parameters, py::const_),
                      [](mio::abm::Location& self, mio::abm::LocalInfectionParameters params) {
                          self.get_infection_parameters() = params;
                      });

    //copying and moving of ranges enabled below, see PYMIO_IGNORE_VALUE_TYPE
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
        .def("add_trip", &mio::abm::TripList::add_trip, py::arg("trip"), py::arg("weekend") = false)
        .def("next_trip", &mio::abm::TripList::get_next_trip, py::arg("weekend") = false)
        .def("num_trips", &mio::abm::TripList::num_trips, py::arg("weekend") = false);

    py::class_<mio::abm::World>(m, "World")
        .def(py::init<int32_t>())
        .def("add_location", &mio::abm::World::add_location, py::arg("location_type"), py::arg("num_cells") = 1)
        .def("add_person", &mio::abm::World::add_person, py::arg("location_id"), py::arg("age_group"),
             py::return_value_policy::reference_internal)
        .def("get_individualized_location",
             py::overload_cast<mio::abm::LocationId>(&mio::abm::World::get_individualized_location, py::const_),
             py::return_value_policy::reference_internal)
        .def("get_rng", &mio::abm::World::get_rng, py::return_value_policy::reference_internal)
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
        .def_readwrite("parameters", &mio::abm::World::parameters)
        .def_property(
            "testing_strategy", py::overload_cast<>(&mio::abm::World::get_testing_strategy, py::const_),
            [](mio::abm::World& self, mio::abm::TestingStrategy strategy) {
                self.get_testing_strategy() = strategy;
            },
            py::return_value_policy::reference_internal);

    m.def(
        "set_viral_load_parameters",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double min_peak,
           double max_peak, double min_incline, double max_incline, double min_decline, double max_decline) {
            infection_params.get<mio::abm::ViralLoadDistributions>()[{variant, age}] =
                mio::abm::ViralLoadDistributionsParameters{
                    {min_peak, max_peak}, {min_incline, max_incline}, {min_decline, max_decline}};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_incubationPeriod",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto incubation_period_params                                      = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::IncubationPeriod>()[{variant, age}] = {incubation_period_params.first,
                                                                                  incubation_period_params.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedNoSymptomsToSymptoms",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedNoSymptomsToSymptoms = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedNoSymptomsToSymptoms>()[{variant, age}] = {
                TimeInfectedNoSymptomsToSymptoms.first, TimeInfectedNoSymptomsToSymptoms.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedNoSymptomsToRecovered",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedNoSymptomsToRecovered = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedNoSymptomsToRecovered>()[{variant, age}] = {
                TimeInfectedNoSymptomsToRecovered.first, TimeInfectedNoSymptomsToRecovered.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedSymptomsToSevere",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedSymptomsToSevere = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedSymptomsToSevere>()[{variant, age}] = {
                TimeInfectedSymptomsToSevere.first, TimeInfectedSymptomsToSevere.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedSymptomsToRecovered",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedSymptomsToRecovered = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedSymptomsToRecovered>()[{variant, age}] = {
                TimeInfectedSymptomsToRecovered.first, TimeInfectedSymptomsToRecovered.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedSevereToRecovered",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedSevereToRecovered = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedSevereToRecovered>()[{variant, age}] = {
                TimeInfectedSevereToRecovered.first, TimeInfectedSevereToRecovered.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedSevereToCritical",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedSevereToCritical = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedSevereToCritical>()[{variant, age}] = {
                TimeInfectedSevereToCritical.first, TimeInfectedSevereToCritical.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedCriticalToRecovered",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedCriticalToRecovered = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedCriticalToRecovered>()[{variant, age}] = {
                TimeInfectedCriticalToRecovered.first, TimeInfectedCriticalToRecovered.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_TimeInfectedCriticalToDead",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double mean,
           double std) {
            auto TimeInfectedCriticalToDead                                              = get_my_and_sigma(mean, std);
            infection_params.get<mio::abm::TimeInfectedCriticalToDead>()[{variant, age}] = {
                TimeInfectedCriticalToDead.first, TimeInfectedCriticalToDead.second};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_infectivity_parameters",
        [](mio::abm::Parameters& infection_params, mio::abm::VirusVariant variant, mio::AgeGroup age, double min_alpha,
           double max_alpha, double min_beta, double max_beta) {
            infection_params.get<mio::abm::InfectivityDistributions>()[{variant, age}] =
                mio::abm::InfectivityDistributionsParameters{{min_alpha, max_alpha}, {min_beta, max_beta}};
        },
        py::return_value_policy::reference_internal);

    m.def(
        "set_seeds",
        [](mio::abm::World& world, int seed) {
            auto rng = mio::RandomNumberGenerator();
            rng.seed({static_cast<uint32_t>(seed)});
            world.get_rng() = rng;
        },
        py::return_value_policy::reference_internal);

    // m.def("get_instance", [](mio::abm::GlobalInfectionParameters& params) {
    //     auto inf_params = params.get<mio::abm::InfectivityDistributions>()[{mio::abm::VirusVariant::Wildtype,
    //                                                                         mio::abm::AgeGroup::Age80plus}];
    //      return inf_params.infectivity_alpha.get_distribution_instance()(inf_params.infectivity_alpha.params);
    // });

    py::class_<mio::abm::Simulation>(m, "Simulation")
        .def(py::init<mio::abm::TimePoint, size_t>())
        .def("advance",
             &mio::abm::Simulation::advance<mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                                         LogPersonsPerLocationAndInfectionTime, LogAgentIds>>)
        //.def("advance", &mio::abm::Simulation::advance)
        //.def_property_readonly("result", &mio::abm::Simulation::get_result)
        .def_property_readonly("world", py::overload_cast<>(&mio::abm::Simulation::get_world),
                               py::return_value_policy::reference_internal);

    py::class_<mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                            LogPersonsPerLocationAndInfectionTime, LogAgentIds>>(m, "History")
        .def(py::init<>())
        .def_property_readonly("log", [](mio::History<mio::DataWriterToMemory, LogTimePoint, LogLocationIds,
                                                      LogPersonsPerLocationAndInfectionTime, LogAgentIds>& self) {
            return self.get_log();
        });
}

// PYMIO_IGNORE_VALUE_TYPE(decltype(std::declval<mio::abm::World>().get_locations()))
// PYMIO_IGNORE_VALUE_TYPE(decltype(std::declval<mio::abm::World>().get_persons()))
