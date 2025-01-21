/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "pybind_util.h"
#include "epidemiology/age_group.h"
#include "epidemiology/damping.h"
#include "epidemiology/contact_matrix.h"
#include "epidemiology/damping_sampling.h"
#include "epidemiology/uncertain_matrix.h"
#include "epidemiology/dynamic_npis.h"
#include "epidemiology/simulation_day.h"
#include "math/integrator.h"
#include "mobility/metapopulation_mobility_instant.h"
#include "utils/date.h"
#include "utils/logging.h"
#include "utils/time_series.h"
#include "utils/parameter_distributions.h"
#include "utils/uncertain_value.h"
#include "utils/index.h"
#include "utils/custom_index_array.h"

//Includes from MEmilio
#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/date.h"
#include "memilio/geography/regions.h"
#include "memilio/epidemiology/contact_matrix.h"
#include "memilio/io/mobility_io.h"
#include "memilio/io/epi_data.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;

PYBIND11_MODULE(_simulation, m)
{
    pymio::bind_parameter_distribution(m, "ParameterDistribution");
    pymio::bind_parameter_distribution_normal(m, "ParameterDistributionNormal");
    pymio::bind_parameter_distribution_uniform(m, "ParameterDistributionUniform");
    pymio::bind_uncertain_value(m, "UncertainValue");

    pymio::bind_CustomIndexArray<mio::UncertainValue<double>, mio::AgeGroup>(m, "AgeGroupArray");
    pymio::bind_class<mio::AgeGroup, pymio::EnablePickling::Required, mio::Index<mio::AgeGroup>>(m, "AgeGroup")
        .def(py::init<size_t>());

    pymio::bind_CustomIndexArray<double, mio::AgeGroup, mio::SimulationDay>(m, "AgeGroupSimulationDayArray");
    pymio::bind_class<mio::SimulationDay, pymio::EnablePickling::Never, mio::Index<mio::SimulationDay>>(m,
                                                                                                        "SimulationDay")
        .def(py::init<size_t>());

    pymio::bind_date(m, "Date");

    auto damping_class = pymio::bind_class<mio::SquareDamping, pymio::EnablePickling::Required>(m, "Damping");
    pymio::bind_damping_members(damping_class);

    auto dampings_class = pymio::bind_class<mio::SquareDampings, pymio::EnablePickling::Required>(m, "Dampings");
    pymio::bind_dampings_members(dampings_class);

    pymio::bind_time_series(m, "TimeSeries");

    pymio::bind_Integrator_Core(m);

    auto contact_matrix_class =
        pymio::bind_class<mio::ContactMatrix, pymio::EnablePickling::Required>(m, "ContactMatrix");
    pymio::bind_damping_expression_members(contact_matrix_class);
    contact_matrix_class.def_property_readonly("num_groups", &mio::ContactMatrix::get_num_groups);

    auto contact_matrix_group_class =
        pymio::bind_class<mio::ContactMatrixGroup, pymio::EnablePickling::Required>(m, "ContactMatrixGroup");
    pymio::bind_damping_expression_group_members(contact_matrix_group_class);
    contact_matrix_group_class.def_property_readonly("num_groups", &mio::ContactMatrixGroup::get_num_groups);

    pymio::bind_damping_sampling(m, "DampingSampling");

    pymio::bind_uncertain_contact_matrix(m, "UncertainContactMatrix");

    auto mobility_damping_class =
        pymio::bind_class<mio::VectorDamping, pymio::EnablePickling::Required>(m, "MobilityDamping");
    pymio::bind_damping_members(mobility_damping_class);

    auto mobility_dampings_class =
        pymio::bind_class<mio::VectorDampings, pymio::EnablePickling::Required>(m, "MobilityDampings");
    pymio::bind_dampings_members(mobility_dampings_class);

    auto mobility_coeffs_class =
        pymio::bind_class<mio::MobilityCoefficients, pymio::EnablePickling::Required>(m, "MobilityCoefficients");
    pymio::bind_damping_expression_members(mobility_coeffs_class);

    auto mobility_coeff_group_class = pymio::bind_class<mio::MobilityCoefficientGroup, pymio::EnablePickling::Required>(
        m, "MobilityCoefficientGroup");
    pymio::bind_damping_expression_group_members(mobility_coeff_group_class);

    pymio::bind_dynamicNPI_members(m, "DynamicNPIs");

    pymio::bind_mobility_parameters(m, "MobilityParameters");
    pymio::bind_mobility(m, "Mobility");
    pymio::bind_mobility_edge(m, "MobilityEdge");
    pymio::bind_mobility_parameter_edge(m, "MobilityParameterEdge");

    m.def(
        "get_state_id_de",
        [](int county) {
            return int(mio::regions::get_state_id(int(mio::regions::CountyId(county))));
        },
        py::arg("county_id"));
    m.def(
        "get_holidays_de",
        [](int state, mio::Date start_date, mio::Date end_date) {
            auto h = mio::regions::get_holidays(mio::regions::StateId(state), start_date, end_date);
            return std::vector<std::pair<mio::Date, mio::Date>>(h.begin(), h.end());
        },
        py::arg("state_id"), py::arg("start_date") = mio::Date(std::numeric_limits<int>::min(), 1, 1),
        py::arg("end_date") = mio::Date(std::numeric_limits<int>::max(), 1, 1));

    m.def(
        "read_mobility_plain",
        [](const std::string& filename) {
            auto result = mio::read_mobility_plain(filename);
            return pymio::check_and_throw(result);
        },
        py::return_value_policy::move);

#ifdef MEMILIO_HAS_JSONCPP
    m.def(
        "get_node_ids",
        [](const std::string& path, bool is_node_for_county) {
            auto result = mio::get_node_ids(path, is_node_for_county);
            return pymio::check_and_throw(result);
        },
        py::return_value_policy::move);
#endif // MEMILIO_HAS_JSONCPP

    pymio::bind_logging(m, "LogLevel");

    m.def("seed_random_number_generator", [] {
        mio::thread_local_rng().seed(mio::RandomNumberGenerator::generate_seeds());
    });

    m.attr("__version__") = "dev";
}
