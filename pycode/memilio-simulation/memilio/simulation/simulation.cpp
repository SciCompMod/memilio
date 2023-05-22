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

#include "pybind_util.h"
#include "epidemiology/damping.h"
#include "epidemiology/contact_matrix.h"
#include "epidemiology/damping_sampling.h"
#include "epidemiology/uncertain_matrix.h"
#include "mobility/metapopulation_mobility_instant.h"
#include "utils/date.h"
#include "utils/logging.h"
#include "utils/time_series.h"
#include "utils/parameter_distributions.h"
#include "utils/uncertain_value.h"

#include "memilio/mobility/metapopulation_mobility_instant.h"
#include "memilio/utils/date.h"
#include "memilio/geography/regions.h"
#include "memilio/epidemiology/contact_matrix.h"

namespace py = pybind11;

PYBIND11_MODULE(_simulation, m)
{
    pymio::bind_date(m, "Date");

    auto damping_class = py::class_<mio::SquareDamping>(m, "Damping");
    pymio::bind_damping_members(damping_class);

    auto dampings_class = py::class_<mio::SquareDampings>(m, "Dampings");
    pymio::bind_dampings_members(dampings_class);

    pymio::bind_time_series(m, "TimeSeries");

    pymio::bind_parameter_distribution(m, "ParameterDistribution");
    pymio::bind_parameter_distribution_normal(m, "ParameterDistributionNormal");
    pymio::bind_parameter_distribution_uniform(m, "ParameterDistributionUniform");

    pymio::bind_uncertain_value(m, "UncertainValue");

    auto contact_matrix_class = py::class_<mio::ContactMatrix>(m, "ContactMatrix");
    pymio::bind_damping_expression_members(contact_matrix_class);
    contact_matrix_class.def_property_readonly("num_groups", &mio::ContactMatrix::get_num_groups);

    auto contact_matrix_group_class = py::class_<mio::ContactMatrixGroup>(m, "ContactMatrixGroup");
    pymio::bind_damping_expression_group_members(contact_matrix_group_class);
    contact_matrix_group_class.def_property_readonly("num_groups", &mio::ContactMatrixGroup::get_num_groups);

    pymio::bind_damping_sampling(m, "DampingSampling");

    pymio::bind_uncertain_contact_matrix(m, "UncertainContactMatrix");

    auto migration_damping_class = py::class_<mio::VectorDamping>(m, "MigrationDamping");
    pymio::bind_damping_members(migration_damping_class);

    auto migration_dampings_class = py::class_<mio::VectorDampings>(m, "MigrationDampings");
    pymio::bind_dampings_members(migration_dampings_class);

    auto migration_coeffs_class = py::class_<mio::MigrationCoefficients>(m, "MigrationCoefficients");
    pymio::bind_damping_expression_members(migration_coeffs_class);

    auto migration_coeff_group_class = py::class_<mio::MigrationCoefficientGroup>(m, "MigrationCoefficientGroup");
    pymio::bind_damping_expression_group_members(migration_coeff_group_class);

    pymio::bind_migration_parameters(m, "MigrationParameters");
    pymio::bind_migration_parameter_edge(m, "MigrationParameterEdge");
    pymio::bind_migration(m, "Migration");
    pymio::bind_migration_edge(m, "MigrationEdge");

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

    pymio::bind_logging(m, "LogLevel");

    m.def("seed_random_number_generator", [] {
        mio::thread_local_rng().seed(mio::RandomNumberGenerator::generate_seeds());
    });

    m.attr("__version__") = "dev";
}
