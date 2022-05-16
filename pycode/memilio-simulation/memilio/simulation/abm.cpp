/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Martin Siggel, Daniel Abele, Martin J. Kuehn, Jan Kleinert
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

#include "templates.h"
#include "utils/custom_index_array.h"
#include "utils/parameter_set.h"
#include "utils/index.h"
#include "abm/abm.h"
#include "pybind11/attr.h"
#include "pybind11/cast.h"
#include "pybind11/pybind11.h"
#include "pybind11/operators.h"
#include <type_traits>

namespace py = pybind11;

PYBIND11_MODULE(_simulation_abm, m)
{
    pymio::iterable_enum<mio::InfectionState>(m, "InfectionState", py::module_local{})
        .value("Susceptible", mio::InfectionState::Susceptible)
        .value("Exposed", mio::InfectionState::Exposed)
        .value("Carrier", mio::InfectionState::Carrier)
        .value("Infected", mio::InfectionState::Infected)
        .value("Infected_Severe", mio::InfectionState::Infected_Severe)
        .value("Infected_Critical", mio::InfectionState::Infected_Critical)
        .value("Recovered_Carrier", mio::InfectionState::Recovered_Carrier)
        .value("Recovered_Infected", mio::InfectionState::Recovered_Infected)
        .value("Dead", mio::InfectionState::Dead);

    pymio::iterable_enum<mio::AbmAgeGroup>(m, "AgeGroup")
        .value("Age0to4", mio::AbmAgeGroup::Age0to4)
        .value("Age5to14", mio::AbmAgeGroup::Age5to14)
        .value("Age15to34", mio::AbmAgeGroup::Age15to34)
        .value("Age35to59", mio::AbmAgeGroup::Age35to59)
        .value("Age60to79", mio::AbmAgeGroup::Age60to79)
        .value("Age80plus", mio::AbmAgeGroup::Age80plus);

    pymio::iterable_enum<mio::VaccinationState>(m, "VaccinationState")
        .value("Unvaccinated", mio::VaccinationState::Unvaccinated)
        .value("Vaccinated", mio::VaccinationState::Vaccinated);

    pymio::iterable_enum<mio::LocationType>(m, "LocationType")
        .value("Home", mio::LocationType::Home)
        .value("School", mio::LocationType::School)
        .value("Work", mio::LocationType::Work)
        .value("SocialEvent", mio::LocationType::SocialEvent)
        .value("BasicsShop", mio::LocationType::BasicsShop)
        .value("Hospital", mio::LocationType::Hospital)
        .value("ICU", mio::LocationType::ICU)
        .value("Car", mio::LocationType::Car)
        .value("PublicTransport", mio::LocationType::PublicTransport)
        .value("TransportWithoutContact", mio::LocationType::TransportWithoutContact);

    py::class_<mio::TestParameters>(m, "TestParameters")
        .def(py::init<double, double>())
        .def_readwrite("sensitivity", &mio::TestParameters::sensitivity)
        .def_readwrite("specificity", &mio::TestParameters::specificity);

    pymio::bind_Index<mio::AbmAgeGroup>(m, "AgeIndex");
    pymio::bind_Index<mio::VaccinationState>(m, "VaccinationIndex");
    pymio::bind_CustomIndexArray<double, mio::AbmAgeGroup, mio::VaccinationState>(m, "_AgeVaccinationParameterArray");
    pymio::bind_ParameterSet<mio::GlobalInfectionParameters>(m, "GlobalInfectionParameters").def(py::init<>());
    pymio::bind_ParameterSet<mio::GlobalTestingParameters>(m, "GlobalTestingParameters").def(py::init<>());
    pymio::bind_ParameterSet<mio::LocalInfectionParameters>(m, "LocalInfectionParameters").def(py::init<>());
    pymio::bind_ParameterSet<mio::AbmMigrationParameters>(m, "MigrationParameters").def(py::init<>());

    py::class_<mio::TimeSpan>(m, "TimeSpan")
        .def(py::init<int>(), py::arg("seconds") = 0)
        .def_property_readonly("seconds", &mio::TimeSpan::seconds)
        .def_property_readonly("hours", &mio::TimeSpan::hours)
        .def_property_readonly("days", &mio::TimeSpan::days)
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

    m.def("seconds", &mio::seconds);
    m.def("minutes", &mio::minutes);
    m.def("hours", &mio::hours);
    m.def("days", &mio::days);

    py::class_<mio::TimePoint>(m, "TimePoint")
        .def(py::init<int>(), py::arg("seconds") = 0)
        .def_property_readonly("seconds", &mio::TimePoint::seconds)
        .def_property_readonly("days", &mio::TimePoint::days)
        .def_property_readonly("hours", &mio::TimePoint::hours)
        .def_property_readonly("day_of_week", &mio::TimePoint::day_of_week)
        .def_property_readonly("hour_of_day", &mio::TimePoint::hour_of_day)
        .def_property_readonly("time_since_midnight", &mio::TimePoint::time_since_midnight)
        .def(py::self == py::self)
        .def(py::self != py::self)
        .def(py::self < py::self)
        .def(py::self <= py::self)
        .def(py::self > py::self)
        .def(py::self >= py::self)
        .def(py::self - py::self)
        .def(py::self + mio::TimeSpan{})
        .def(py::self += mio::TimeSpan{})
        .def(py::self - mio::TimeSpan{})
        .def(py::self -= mio::TimeSpan{});

    py::class_<mio::LocationId>(m, "LocationId")
        .def(py::init([](uint32_t idx, mio::LocationType type) {
            return mio::LocationId{idx, type};
        }))
        .def_readwrite("index", &mio::LocationId::index)
        .def_readwrite("type", &mio::LocationId::type)
        .def(py::self == py::self)
        .def(py::self != py::self);

    py::class_<mio::InfectionProperties>(m, "InfectionProperties")
        .def_readwrite("state", &mio::InfectionProperties::state)
        .def_readwrite("detected", &mio::InfectionProperties::detected)
        .def(py::self == py::self)
        .def(py::self != py::self);

    py::class_<mio::Person>(m, "Person")
        .def("set_assigned_location", py::overload_cast<mio::LocationId>(&mio::Person::set_assigned_location))
        .def_property("infection_state", &mio::Person::get_infection_state, &mio::Person::set_infection_state)
        .def_property_readonly("vaccination_state", &mio::Person::get_vaccination_state)
        .def_property_readonly("location_id", &mio::Person::get_location_id)
        .def_property_readonly("age", &mio::Person::get_age)
        .def_property_readonly("is_in_quarantine", &mio::Person::is_in_quarantine);

    py::class_<mio::TestingScheme>(m, "TestingScheme")
        .def(py::init<mio::TimeSpan, double>(), py::arg("interval"), py::arg("probability"))
        .def_property("interval", &mio::TestingScheme::get_interval, &mio::TestingScheme::set_interval)
        .def_property("probability", &mio::TestingScheme::get_probability, &mio::TestingScheme::set_probability);

    py::class_<mio::Location>(m, "Location")
        .def_property_readonly("type", &mio::Location::get_type)
        .def_property_readonly("index", &mio::Location::get_index)
        .def_property("infection_parameters", py::overload_cast<>(&mio::Location::get_infection_parameters, py::const_),
                      [](mio::Location& self, mio::LocalInfectionParameters params) {
                          self.get_infection_parameters() = params;
                      })
        .def_property("testing_scheme", &mio::Location::get_testing_scheme,
                      [](mio::Location& self, mio::TestingScheme scheme) {
                          self.set_testing_scheme(scheme.get_interval(), scheme.get_probability());
                      });

    pymio::bind_Range<decltype(std::declval<mio::World>().get_locations())>(m, "_WorldLocationsRange");
    pymio::bind_Range<decltype(std::declval<mio::World>().get_persons())>(m, "_WorldPersonsRange");

    py::class_<mio::Trip>(m, "Trip")
        .def(py::init<uint32_t, mio::TimePoint, mio::LocationId, mio::LocationId, std::vector<uint32_t>>(),
             py::arg("person_id"), py::arg("time"), py::arg("destination"), py::arg("origin"),
             py::arg("cells") = std::vector<uint32_t>())
        .def_readwrite("person_id", &mio::Trip::person_id)
        .def_readwrite("time", &mio::Trip::time)
        .def_readwrite("destination", &mio::Trip::migration_destination)
        .def_readwrite("origin", &mio::Trip::migration_origin)
        .def_readwrite("cells", &mio::Trip::cells);

    py::class_<mio::TripList>(m, "TripList")
        .def(py::init<>())
        .def("add_trip", &mio::TripList::add_trip)
        .def_property_readonly("next_trip", &mio::TripList::get_next_trip)
        .def_property_readonly("num_trips", &mio::TripList::num_trips);

    py::class_<mio::World>(m, "World")
        .def(py::init<mio::GlobalInfectionParameters>(),
             py::arg("infection_parameters") = mio::GlobalInfectionParameters{})
        .def("add_location", &mio::World::add_location, py::arg("location_type"), py::arg("num_cells") = 0,
             py::return_value_policy::reference_internal)
        .def("add_person", &mio::World::add_person, py::return_value_policy::reference_internal)
        .def_property_readonly("locations", &mio::World::get_locations,
                               py::keep_alive<1, 0>{}) //keep this world alive while contents are referenced in ranges
        .def_property_readonly("persons", &mio::World::get_persons, py::keep_alive<1, 0>{})
        .def_property(
            "trip_list", py::overload_cast<>(&mio::World::get_trip_list),
            [](mio::World& self, const mio::TripList& list) {
                self.get_trip_list() = list;
            },
            py::return_value_policy::reference_internal)
        .def_property("use_migration_rules", py::overload_cast<>(&mio::World::use_migration_rules, py::const_),
                      py::overload_cast<bool>(&mio::World::use_migration_rules))
        .def_property(
            "infection_parameters", py::overload_cast<>(&mio::World::get_global_infection_parameters, py::const_),
            [](mio::World& self, mio::GlobalInfectionParameters params) {
                self.get_global_infection_parameters() = params;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "migration_parameters", py::overload_cast<>(&mio::World::get_migration_parameters, py::const_),
            [](mio::World& self, mio::AbmMigrationParameters params) {
                self.get_migration_parameters() = params;
            },
            py::return_value_policy::reference_internal)
        .def_property(
            "testing_parameters", py::overload_cast<>(&mio::World::get_global_testing_parameters, py::const_),
            [](mio::World& self, mio::GlobalTestingParameters params) {
                self.get_global_testing_parameters() = params;
            },
            py::return_value_policy::reference_internal);

    py::class_<mio::AbmSimulation>(m, "Simulation")
        .def(py::init<mio::TimePoint>())
        .def("advance", &mio::AbmSimulation::advance)
        .def_property_readonly("result", &mio::AbmSimulation::get_result)
        .def_property_readonly("world", py::overload_cast<>(&mio::AbmSimulation::get_world));
}