/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
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
#ifndef MIO_EPI_REGIONS_H
#define MIO_EPI_REGIONS_H

#include "memilio/utils/date.h"
#include "memilio/utils/stl_util.h"
#include "memilio/utils/type_safe.h"
#include "memilio/mobility/graph.h"
#include "memilio/mobility/mobility.h"

#include "boost/filesystem.hpp"

namespace fs = boost::filesystem;

namespace mio
{
/**
 * Contains utilities that depend on geographical regions.
 */
namespace regions
{
/**
     * Germany.
     */
namespace de
{
/**
         * Id of a state.
         * 1 = Schleswig-Holstein
         * 2 = Hamburg
         * 3 = Niedersachsen
         * 4 = Bremen
         * 5 = Nordrhein-Westfalen
         * 6 = Hessen
         * 7 = Rheinland-Pfalz
         * 8 = Baden-Württemberg
         * 9 = Bayern
         * 10 = Saarland
         * 11 = Berlin
         * 12 = Brandenburg
         * 13 = Mecklenburg-Vorpommern
         * 14 = Sachsen
         * 15 = Sachsen-Anhalt
         * 16 = Thüringen
         */
DECL_TYPESAFE(int, StateId);

/**
         * Id of a county.
         * Format ssxxx where ss is the id of the state that the county is in (first s may be 0) and xxx are other digits.
         * Ids are generally not consecutive, even within one state.
         */
DECL_TYPESAFE(int, CountyId);

DECL_TYPESAFE(int, DistrictId);

/**
         * get the id of the state that the specified county is in. 
         * @param county a county id.
         */
StateId get_state_id(CountyId county);

/**
         * get the holidays in a german state.
         * @param state id of the state.
         * @return range of pairs of start and end dates of holiday periods, sorted by start date.
         */
Range<std::pair<std::vector<std::pair<Date, Date>>::const_iterator, std::vector<std::pair<Date, Date>>::const_iterator>>
get_holidays(StateId state);

/**
         * get the holidays in a german state in a given time period.
         * The returned periods may not be completely included in the queried period,
         * they may only partially overlap.
         * @param state id of the state.
         * @param start_date start of the queried period.
         * @param end_date end of the queried period.
         * @return range of pairs of start and end dates of holiday periods, sorted by start date.
         */
Range<std::pair<std::vector<std::pair<Date, Date>>::const_iterator, std::vector<std::pair<Date, Date>>::const_iterator>>
get_holidays(StateId state, Date start_date, Date end_date);

} // namespace de

IOResult<std::vector<int>> get_node_ids(const std::string& path);

template <typename TNT, typename ContactPattern, class Model, class Parameters, class ReadFunction>
IOResult<void> set_nodes(const Parameters& params, Date start_date, Date end_date, const fs::path& data_dir,
                         Graph<Model, MigrationParameters>& params_graph, ReadFunction&& read_func,
                         const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                         double tnt_capacity_factor, int num_days = 0, bool export_time_series = false)
{
    BOOST_OUTCOME_TRY(county_ids, mio::get_county_ids((data_dir / "pydata" / "Germany").string()));
    std::vector<Model> counties(county_ids.size(), Model(int(size_t(params.get_num_groups()))));
    for (auto& county : counties) {
        county.parameters = params;
    }

    BOOST_OUTCOME_TRY(read_func(counties, start_date, county_ids, scaling_factor_inf, scaling_factor_icu,
                                (data_dir / "pydata" / "Germany").string(), num_days, export_time_series));

    for (size_t county_idx = 0; county_idx < counties.size(); ++county_idx) {

        auto tnt_capacity = counties[county_idx].populations.get_total() * tnt_capacity_factor;

        //local parameters
        auto& tnt_value = counties[county_idx].parameters.get<TNT>();
        tnt_value       = UncertainValue(0.5 * (1.2 * tnt_capacity + 0.8 * tnt_capacity));
        tnt_value.set_distribution(mio::ParameterDistributionUniform(0.8 * tnt_capacity, 1.2 * tnt_capacity));

        //holiday periods
        auto holiday_periods =
            de::get_holidays(de::get_state_id(de::CountyId(county_ids[county_idx])), start_date, end_date);
        auto& contacts = counties[county_idx].parameters.get<ContactPattern>();
        contacts.get_school_holidays() =
            std::vector<std::pair<mio::SimulationTime, mio::SimulationTime>>(holiday_periods.size());
        std::transform(
            holiday_periods.begin(), holiday_periods.end(), contacts.get_school_holidays().begin(), [=](auto& period) {
                return std::make_pair(mio::SimulationTime(mio::get_offset_in_days(period.first, start_date)),
                                      mio::SimulationTime(mio::get_offset_in_days(period.second, start_date)));
            });

        //uncertainty in populations
        for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
            for (auto j = Index<Model::Compartments>(0); j < Model::Compartments::Count; ++j) {
                auto& compartment_value = counties[county_idx].populations[{i, j}];
                compartment_value =
                    UncertainValue(0.5 * (1.1 * double(compartment_value) + 0.9 * double(compartment_value)));
                compartment_value.set_distribution(mio::ParameterDistributionUniform(0.9 * double(compartment_value),
                                                                                     1.1 * double(compartment_value)));
            }
        }

        params_graph.add_node(county_ids[county_idx], counties[county_idx]);
    }
    return success();
}

template <class ContactLocation, class Model, class InfectionState>
IOResult<void> set_edges(const fs::path& data_dir, Graph<Model, MigrationParameters>& params_graph,
                         std::initializer_list<InfectionState>& migrating_compartments, size_t contact_locations_size)
{
    // mobility between nodes
    BOOST_OUTCOME_TRY(mobility_data_commuter,
                      mio::read_mobility_plain((data_dir / "mobility" / "commuter_migration_scaled.txt").string()));
    BOOST_OUTCOME_TRY(mobility_data_twitter,
                      mio::read_mobility_plain((data_dir / "mobility" / "twitter_scaled_1252.txt").string()));
    if (mobility_data_commuter.rows() != Eigen::Index(params_graph.nodes().size()) ||
        mobility_data_commuter.cols() != Eigen::Index(params_graph.nodes().size()) ||
        mobility_data_twitter.rows() != Eigen::Index(params_graph.nodes().size()) ||
        mobility_data_twitter.cols() != Eigen::Index(params_graph.nodes().size())) {
        return mio::failure(mio::StatusCode::InvalidValue,
                            "Mobility matrices do not have the correct size. You may need to run "
                            "transformMobilitydata.py from pycode memilio epidata package.");
    }

    for (size_t county_idx_i = 0; county_idx_i < params_graph.nodes().size(); ++county_idx_i) {
        for (size_t county_idx_j = 0; county_idx_j < params_graph.nodes().size(); ++county_idx_j) {
            auto& populations = params_graph.nodes()[county_idx_i].property.populations;
            // mobility coefficients have the same number of components as the contact matrices.
            // so that the same NPIs/dampings can be used for both (e.g. more home office => fewer commuters)
            auto mobility_coeffs = MigrationCoefficientGroup(contact_locations_size, populations.numel());

            //commuters
            auto working_population = 0.0;
            auto min_commuter_age   = AgeGroup(2);
            auto max_commuter_age   = AgeGroup(4); //this group is partially retired, only partially commutes
            for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                working_population += populations.get_group_total(age) * (age == max_commuter_age ? 0.33 : 1.0);
            }
            auto commuter_coeff_ij = mobility_data_commuter(county_idx_i, county_idx_j) /
                                     working_population; //data is absolute numbers, we need relative
            for (auto age = min_commuter_age; age <= max_commuter_age; ++age) {
                for (auto compartment : migrating_compartments) {
                    auto coeff_index = populations.get_flat_index({age, compartment});
                    mobility_coeffs[size_t(ContactLocation::Work)].get_baseline()[coeff_index] =
                        commuter_coeff_ij * (age == max_commuter_age ? 0.33 : 1.0);
                }
            }
            //others
            auto total_population = populations.get_total();
            auto twitter_coeff    = mobility_data_twitter(county_idx_i, county_idx_j) /
                                 total_population; //data is absolute numbers, we need relative
            for (auto age = AgeGroup(0); age < populations.size<mio::AgeGroup>(); ++age) {
                for (auto compartment : migrating_compartments) {
                    auto coeff_idx = populations.get_flat_index({age, compartment});
                    mobility_coeffs[size_t(ContactLocation::Other)].get_baseline()[coeff_idx] = twitter_coeff;
                }
            }

            //only add edges with mobility above thresholds for performance
            //thresholds are chosen empirically so that more than 99% of mobility is covered, approx. 1/3 of the edges
            if (commuter_coeff_ij > 4e-5 || twitter_coeff > 1e-5) {
                params_graph.add_edge(county_idx_i, county_idx_j, std::move(mobility_coeffs));
            }
        }
    }

    return success();
}

} // namespace regions
} // namespace mio

#endif //MIO_EPI_REGIONS_H
