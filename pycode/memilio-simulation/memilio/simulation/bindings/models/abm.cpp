/* 
* Copyright (C) 2020-2026 MEmilio
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

//Includes from pymio
#include "abm/person_id.h"
#include "pybind_util.h"
#include "utils/custom_index_array.h"
#include "utils/parameter_set.h"
#include "utils/index.h"

//Includes from MEmilio
#include "abm/pcr_surveillance.h"
#include "abm/simulation.h"
#include "forward_pass.h"

#include "pybind11/attr.h"
#include "pybind11/cast.h"
#include "pybind11/pybind11.h"
#include "pybind11/operators.h"
#include <pybind11/eigen.h>
#include <pybind11/stl.h>
#include <cstdint>
#include <type_traits>

namespace py = pybind11;

PYBIND11_MODULE(_simulation_abm, m)
{
    pymio::iterable_enum<mio::abm::InfectionState>(m, "InfectionState")
        .value("Susceptible", mio::abm::InfectionState::Susceptible)
        .value("Exposed", mio::abm::InfectionState::Exposed)
        .value("InfectedNoSymptoms", mio::abm::InfectionState::InfectedNoSymptoms)
        .value("InfectedSymptoms", mio::abm::InfectionState::InfectedSymptoms)
        .value("InfectedSevere", mio::abm::InfectionState::InfectedSevere)
        .value("InfectedCritical", mio::abm::InfectionState::InfectedCritical)
        .value("Recovered", mio::abm::InfectionState::Recovered)
        .value("Dead", mio::abm::InfectionState::Dead);

    pymio::iterable_enum<mio::abm::ProtectionType>(m, "ProtectionType")
        .value("NoProtection", mio::abm::ProtectionType::NoProtection)
        .value("NaturalInfection", mio::abm::ProtectionType::NaturalInfection)
        .value("GenericVaccine", mio::abm::ProtectionType::GenericVaccine);

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
        .value("TransportWithoutContact", mio::abm::LocationType::TransportWithoutContact);

    pymio::iterable_enum<mio::abm::TestType>(m, "TestType")
        .value("Generic", mio::abm::TestType::Generic)
        .value("Antigen", mio::abm::TestType::Antigen)
        .value("PCR", mio::abm::TestType::PCR);

    pymio::bind_class<mio::abm::TimeSpan, pymio::EnablePickling::Never>(m, "TimeSpan")
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

    pymio::bind_class<mio::abm::TimePoint, pymio::EnablePickling::Never>(m, "TimePoint")
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

    pymio::bind_class<mio::abm::TestParameters, pymio::EnablePickling::Never>(m, "TestParameters")
        .def(py::init<double, double, mio::abm::TimeSpan, mio::abm::TestType>())
        .def_readwrite("sensitivity", &mio::abm::TestParameters::sensitivity)
        .def_readwrite("specificity", &mio::abm::TestParameters::specificity)
        .def_readwrite("required_time", &mio::abm::TestParameters::required_time)
        .def_readwrite("type", &mio::abm::TestParameters::type);

    pymio::bind_CustomIndexArray<mio::UncertainValue<double>, mio::abm::VirusVariant, mio::AgeGroup>(
        m, "_AgeParameterArray");
    pymio::bind_CustomIndexArray<mio::abm::TestParameters, mio::abm::TestType>(m, "_TestData");
    pymio::bind_Index<mio::abm::ProtectionType>(m, "ProtectionTypeIndex");
    pymio::bind_ParameterSet<mio::abm::ParametersBase, pymio::EnablePickling::Never>(m, "ParametersBase");
    pymio::bind_class<mio::abm::Parameters, pymio::EnablePickling::Never, mio::abm::ParametersBase>(m, "Parameters")
        .def(py::init<int>())
        .def("check_constraints", &mio::abm::Parameters::check_constraints);

    pymio::bind_ParameterSet<mio::abm::LocalInfectionParameters, pymio::EnablePickling::Never>(
        m, "LocalInfectionParameters")
        .def(py::init<size_t>());

    pymio::bind_class<mio::abm::LocationId, pymio::EnablePickling::Never>(m, "LocationId")
        .def(py::init<uint32_t>(), py::arg("id"))
        .def("index", &mio::abm::LocationId::get);

    pymio::bind_class<mio::abm::PersonId, pymio::EnablePickling::Never>(m, "PersonId")
        .def(py::init<uint64_t>(), py::arg("id"))
        .def("index", &mio::abm::PersonId::get);

    pymio::bind_class<mio::abm::Person, pymio::EnablePickling::Never>(m, "Person")
        .def("set_assigned_location", py::overload_cast<mio::abm::LocationType, mio::abm::LocationId, int>(
                                          &mio::abm::Person::set_assigned_location))
        .def_property_readonly("location", py::overload_cast<>(&mio::abm::Person::get_location, py::const_))
        .def_property_readonly("age", &mio::abm::Person::get_age)
        .def_property_readonly("is_in_quarantine", &mio::abm::Person::is_in_quarantine);

    pymio::bind_class<mio::abm::TestingCriteria, pymio::EnablePickling::Never>(m, "TestingCriteria")
        .def(py::init<const std::vector<mio::AgeGroup>&, const std::vector<mio::abm::InfectionState>&>(),
             py::arg("age_groups"), py::arg("infection_states"));

    pymio::bind_class<mio::abm::TestingScheme, pymio::EnablePickling::Never>(m, "TestingScheme")
        .def(py::init<const mio::abm::TestingCriteria&, mio::abm::TimeSpan, mio::abm::TimePoint, mio::abm::TimePoint,
                      const mio::abm::TestParameters&, double>(),
             py::arg("testing_criteria"), py::arg("testing_validity_period"), py::arg("start_date"),
             py::arg("end_date"), py::arg("test_parameters"), py::arg("probability"));

    pymio::bind_class<mio::abm::ProtectionEvent, pymio::EnablePickling::Never>(m, "ProtectionEvent")
        .def(py::init<mio::abm::ProtectionType, mio::abm::TimePoint>(), py::arg("type"), py::arg("time"))
        .def_readwrite("type", &mio::abm::ProtectionEvent::type)
        .def_readwrite("time", &mio::abm::ProtectionEvent::time);

    pymio::bind_class<mio::abm::TestingStrategy, pymio::EnablePickling::Never>(m, "TestingStrategy")
        .def(py::init<const std::vector<mio::abm::TestingStrategy::LocalStrategy>&,
                      const std::vector<mio::abm::TestingStrategy::LocalStrategy>&>());

    pymio::bind_class<mio::abm::Location, pymio::EnablePickling::Never>(m, "Location")
        .def_property_readonly("type", &mio::abm::Location::get_type)
        .def_property_readonly("id", &mio::abm::Location::get_id)
        .def_property("infection_parameters",
                      py::overload_cast<>(&mio::abm::Location::get_infection_parameters, py::const_),
                      [](mio::abm::Location& self, mio::abm::LocalInfectionParameters params) {
                          self.get_infection_parameters() = params;
                      });

    //copying and moving of ranges enabled below, see PYMIO_IGNORE_VALUE_TYPE
    pymio::bind_Range<decltype(std::declval<const mio::abm::Model>().get_locations())>(m, "_ModelLocationsRange");
    pymio::bind_Range<decltype(std::declval<const mio::abm::Model>().get_persons())>(m, "_ModelPersonsRange");

    pymio::bind_class<mio::abm::Trip, pymio::EnablePickling::Never>(m, "Trip")
        .def(py::init<mio::abm::PersonId, mio::abm::TimePoint, mio::abm::LocationId>(), py::arg("person_id"),
             py::arg("time"), py::arg("destination"))
        .def_readwrite("person_id", &mio::abm::Trip::person_id)
        .def_readwrite("trip_time", &mio::abm::Trip::trip_time)
        .def_readwrite("destination", &mio::abm::Trip::destination);

    pymio::bind_class<mio::abm::TripList, pymio::EnablePickling::Never>(m, "TripList")
        .def(py::init<>())
        .def("add_trips", &mio::abm::TripList::add_trips, py::arg("trips") = std::vector<mio::abm::Trip>())
        .def("next_trip", &mio::abm::TripList::get_next_trip)
        .def("num_trips", &mio::abm::TripList::num_trips);

    pymio::bind_class<mio::abm::Model, pymio::EnablePickling::Never>(m, "Model")
        .def(py::init<int32_t>())
        .def("add_location", &mio::abm::Model::add_location, py::arg("location_type"), py::arg("num_cells") = 1)
        .def("add_person", py::overload_cast<mio::abm::LocationId, mio::AgeGroup>(&mio::abm::Model::add_person),
             py::arg("location_id"), py::arg("age_group"))
        .def("assign_location",
             py::overload_cast<mio::abm::PersonId, mio::abm::LocationId>(&mio::abm::Model::assign_location),
             py::arg("person_id"), py::arg("location_id"))
        .def_property_readonly("locations", py::overload_cast<>(&mio::abm::Model::get_locations, py::const_),
                               py::keep_alive<1, 0>{}) //keep this model alive while contents are referenced in ranges
        .def_property_readonly("persons", py::overload_cast<>(&mio::abm::Model::get_persons, py::const_),
                               py::keep_alive<1, 0>{})
        .def_property(
            "trip_list", py::overload_cast<>(&mio::abm::Model::get_trip_list),
            [](mio::abm::Model& self, const mio::abm::TripList& list) {
                self.get_trip_list() = list;
            },
            py::return_value_policy::reference_internal)
        .def_property("use_mobility_rules", py::overload_cast<>(&mio::abm::Model::use_mobility_rules, py::const_),
                      py::overload_cast<bool>(&mio::abm::Model::use_mobility_rules))
        .def_readwrite("parameters", &mio::abm::Model::parameters)
        .def_property(
            "testing_strategy", py::overload_cast<>(&mio::abm::Model::get_testing_strategy, py::const_),
            [](mio::abm::Model& self, mio::abm::TestingStrategy strategy) {
                self.get_testing_strategy() = strategy;
            },
            py::return_value_policy::reference_internal);

    pymio::bind_class<mio::abm::Simulation<>, pymio::EnablePickling::Never>(m, "Simulation")
        .def(py::init<mio::abm::TimePoint, size_t>())
        .def("advance",
             static_cast<void (mio::abm::Simulation<>::*)(mio::abm::TimePoint)>(&mio::abm::Simulation<>::advance),
             py::arg("tmax"))
        .def_property_readonly("model", py::overload_cast<>(&mio::abm::Simulation<>::get_model));

    m.attr("__version__") = "dev";

    py::class_<mio::abm::TestingBudget>(m, "TestingBudget",
      "Active-surveillance testing design d: event_budget_fraction is the share of the "
      "population sampled at each surveillance event. test_period_days is "
      "the cadence (sample every N days), test_hour the hour of day, reporting_delay_days the lab "
      "turnaround, min_retest_gap_days the minimum time between two PCR tests of the same Person "
      "event_budget_fraction = 0 disables surveillance. ")
        .def(py::init([](ScalarType event_budget_fraction, int test_period_days, int test_hour,
                         ScalarType reporting_delay_days, ScalarType min_retest_gap_days) {
                 mio::abm::TestingBudget design;
                 design.event_budget_fraction = event_budget_fraction;
                 design.test_period_days      = test_period_days;
                 design.test_hour             = test_hour;
                 design.reporting_delay       = mio::abm::days(reporting_delay_days);
                 design.min_retest_gap        = mio::abm::days(min_retest_gap_days);
                 return design;
             }),
             py::arg("event_budget_fraction") = 0.01, py::arg("test_period_days") = 1, py::arg("test_hour") = 8,
             py::arg("reporting_delay_days") = 1.0, py::arg("min_retest_gap_days") = 2.0)
        .def_readwrite("event_budget_fraction", &mio::abm::TestingBudget::event_budget_fraction)
        .def_readwrite("test_period_days", &mio::abm::TestingBudget::test_period_days)
        .def_readwrite("test_hour", &mio::abm::TestingBudget::test_hour)
        .def(py::pickle(
            [](const mio::abm::TestingBudget& d) { // __getstate__
                return py::make_tuple(d.event_budget_fraction, d.test_period_days, d.test_hour,
                                      d.reporting_delay.seconds(), d.min_retest_gap.seconds());
            },
            [](py::tuple t) { // __setstate__
                if (t.size() != 5) {
                    throw std::runtime_error("Invalid TestingBudget pickled state.");
                }
                mio::abm::TestingBudget d;
                d.event_budget_fraction = t[0].cast<ScalarType>();
                d.test_period_days      = t[1].cast<int>();
                d.test_hour             = t[2].cast<int>();
                d.reporting_delay       = mio::abm::TimeSpan(t[3].cast<int>());
                d.min_retest_gap        = mio::abm::TimeSpan(t[4].cast<int>());
                return d;
            }));

    py::class_<ABMPopulation>(m, "ABMPopulation",
      "Pre-built ABM population from CityBuilder: German-demographics-based households, "
      "schools, workplaces, shops, events, and hospital/ICU. Construct once, reuse across many "
      "forward_pass() calls to avoid rebuilding the structure.")
        .def(py::init<int>(), py::arg("total_population") = 500,
             py::call_guard<py::gil_scoped_release>());

    m.def("forward_pass", &forward_pass,
      py::arg("population"), py::arg("params"), py::arg("design") = mio::abm::TestingBudget{0.01},
      py::call_guard<py::gil_scoped_release>(),
      "Run the ABM forward pass using a pre-built population, with active PCR surveillance "
      "plus a fixed diagnostic care-seeking process. `params` is a dict mapping parameter "
      "names to values (e.g. {'beta': 1.0, 'kappa': 5.0}); recognised keys: "
      "'beta', 'kappa', 'time_presymptomatic', 'time_asymptomatic_recovery', 'symptom_prob', "
      "'quarantine_compliance'. "
      "Unknown keys raise ValueError; omitted keys keep the model default. Releases the GIL "
      "to allow concurrent calls from multiple Python threads sharing one ABMPopulation. "
      "Returns a dict of raw resolved tests: "
      "'positives' (n_pos, 6) [test_day, person_id, age, location, ct, source] one row per positive test "
      "(source 0 = survey / active surveillance = the analysed stream X, 1 = diagnostic); "
      "'negatives' (n_groups, 5) [test_day, source, age, location, count] negatives aggregated to "
      "counts; 'n_ever_infected' (1, 1) an epidemic-size diagnostic.");
}

PYMIO_IGNORE_VALUE_TYPE(decltype(std::declval<mio::abm::Model>().get_locations()))
PYMIO_IGNORE_VALUE_TYPE(decltype(std::declval<mio::abm::Model>().get_persons()))
