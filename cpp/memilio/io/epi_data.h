/* 
* Copyright (C) 2020-2022 German Aerospace Center (DLR-SC)
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
#ifndef MEMILIO_IO_EPI_DATA_H
#define MEMILIO_IO_EPI_DATA_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP

#include "memilio/epidemiology/age_group.h"
#include "memilio/geography/regions.h"
#include "memilio/io/io.h"
#include "memilio/io/json_serializer.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/date.h"
#include "memilio/io/mobility_io.h"

#include "json/value.h"
#include <string>
#include <vector>

namespace mio
{

/**
 * Date that serializes into a string.
 */
class StringDate : public Date
{
public:
    using Date::Date;

    StringDate(const Date& other)
        : Date(other)
    {
    }

    template <class IoContext>
    static IOResult<StringDate> deserialize(IoContext& io)
    {
        auto str = mio::deserialize(io, Tag<std::string>{});
        return apply(
            io,
            [](auto&& str_) -> IOResult<StringDate> {
                BOOST_OUTCOME_TRY(date, parse_date(str_));
                return success(date);
            },
            str);
    }
};

/**
 * Represents the entries of a confirmed cases data file, e.g., from RKI.
 * Number of confirmed, recovered and deceased in a region on a specific date.
 * Region can be a county, a state, or a country. If it is a country, both
 * state_id and county_id will be empty.
 */
class ConfirmedCasesDataEntry
{
public:
    static const std::array<const char*, 6> age_group_names;

    double num_confirmed;
    double num_recovered;
    double num_deaths;
    Date date;
    AgeGroup age_group;
    boost::optional<regions::StateId> state_id;
    boost::optional<regions::CountyId> county_id;

    template <class IOContext>
    static IOResult<ConfirmedCasesDataEntry> deserialize(IOContext& io)
    {
        auto obj           = io.expect_object("ConfirmedCasesDataEntry");
        auto num_confirmed = obj.expect_element("Confirmed", Tag<double>{});
        auto num_recovered = obj.expect_element("Recovered", Tag<double>{});
        auto num_deaths    = obj.expect_element("Deaths", Tag<double>{});
        auto date          = obj.expect_element("Date", Tag<StringDate>{});
        auto age_group_str = obj.expect_element("Age_RKI", Tag<std::string>{});
        auto state_id      = obj.expect_optional("ID_State", Tag<regions::StateId>{});
        auto county_id     = obj.expect_optional("ID_County", Tag<regions::CountyId>{});
        return apply(
            io,
            [](auto&& nc, auto&& nr, auto&& nd, auto&& d, auto&& a_str, auto&& sid,
               auto&& cid) -> IOResult<ConfirmedCasesDataEntry> {
                auto a  = AgeGroup(0);
                auto it = std::find(age_group_names.begin(), age_group_names.end(), a_str);
                if (it != age_group_names.end()) {
                    a = AgeGroup(size_t(it - age_group_names.begin()));
                }
                else if (a_str == "unknown") {
                    a = AgeGroup(age_group_names.size());
                }
                else {
                    return failure(StatusCode::InvalidValue, "Invalid confirmed cases data age group.");
                }
                return success(ConfirmedCasesDataEntry{nc, nr, nd, d, a, sid, cid});
            },
            num_confirmed, num_recovered, num_deaths, date, age_group_str, state_id, county_id);
    }
};

/**
 * Read list of ConfirmedCasesDataEntry from json.
 * @param jsvalue json value, must be an array of objects, objects must match ConfirmedCasesDataEntry.
 * @return list of entries; entries of unknown age group are omitted.
 */
inline IOResult<std::vector<ConfirmedCasesDataEntry>> deserialize_confirmed_cases_data(const Json::Value& jsvalue)
{
    BOOST_OUTCOME_TRY(cases_data, deserialize_json(jsvalue, Tag<std::vector<ConfirmedCasesDataEntry>>{}));
    //filter entries with unknown age group
    auto it = std::remove_if(cases_data.begin(), cases_data.end(), [](auto&& rki_entry) {
        return rki_entry.age_group >= AgeGroup(ConfirmedCasesDataEntry::age_group_names.size());
    });
    cases_data.erase(it, cases_data.end());
    return success(std::move(cases_data));
}

/**
 * Read list of ConfirmedCasesDataEntry from a json file.
 * @param filename name of the json file. File content must be an array of objects, objects must match ConfirmedCasesDataEntry.
 * @return list of entries; entries of unknown age group are omitted.
 */
inline IOResult<std::vector<ConfirmedCasesDataEntry>> read_confirmed_cases_data(const std::string& filename)
{
    BOOST_OUTCOME_TRY(jsvalue, read_json(filename));
    return deserialize_confirmed_cases_data(jsvalue);
}

/**
 * Represents entries in a DIVI data file.
 * Number of persons in the ICU in a region on a specific date.
 * Region can be a county, a state, or a country. If it is a country, both
 * state_id and county_id will be empty.
 */
class DiviEntry
{
public:
    double num_icu;
    Date date;
    boost::optional<regions::StateId> state_id;
    boost::optional<regions::CountyId> county_id;

    template <class IoContext>
    static IOResult<DiviEntry> deserialize(IoContext& io)
    {
        auto obj       = io.expect_object("DiviEntry");
        auto num_icu   = obj.expect_element("ICU", Tag<double>{});
        auto date      = obj.expect_element("Date", Tag<StringDate>{});
        auto state_id  = obj.expect_optional("ID_State", Tag<regions::StateId>{});
        auto county_id = obj.expect_optional("ID_County", Tag<regions::CountyId>{});
        return apply(
            io,
            [](auto&& ni, auto&& d, auto&& sid, auto&& cid) {
                return DiviEntry{ni, d, sid, cid};
            },
            num_icu, date, state_id, county_id);
    }
};

/**
 * Deserialize a list of DiviEntry from json.
 * @param jsvalue Json value that contains DIVI data.
 * @return list of DiviEntry.
 */
inline IOResult<std::vector<DiviEntry>> deserialize_divi_data(const Json::Value& jsvalue)
{
    return deserialize_json(jsvalue, Tag<std::vector<DiviEntry>>{});
}

/**
 * deserialize a list of DiviEntry from json.
 * @param filename Json file that contains DIVI data.
 * @return list of DiviEntry.
 */
inline IOResult<std::vector<DiviEntry>> read_divi_data(const std::string& filename)
{
    return read_json(filename, Tag<std::vector<DiviEntry>>{});
}

namespace details
{
//check all results in a vector and unpack each
template <class T>
IOResult<std::vector<T>> unpack_all(const std::vector<IOResult<T>>& v)
{
    std::vector<T> w;
    w.reserve(v.size());
    for (auto&& r : v) {
        BOOST_OUTCOME_TRY(t, r);
        w.push_back(t);
    }
    return success(w);
}
} // namespace details

/**
 * Represents an entry of a population data file.
 * Population per age group in a region.
 * Region can be a county, a state, or a country. If it is a country, both
 * state_id and county_id will be empty.
 */
class PopulationDataEntry
{
public:
    static const std::array<const char*, 11> age_group_names;

    CustomIndexArray<double, AgeGroup> population;
    boost::optional<regions::StateId> state_id;
    boost::optional<regions::CountyId> county_id;

    template <class IoContext>
    static IOResult<PopulationDataEntry> deserialize(IoContext& io)
    {
        auto obj       = io.expect_object("PopulationDataEntry");
        auto state_id  = obj.expect_optional("ID_State", Tag<regions::StateId>{});
        auto county_id = obj.expect_optional("ID_County", Tag<regions::CountyId>{});
        std::vector<IOResult<double>> age_groups;
        age_groups.reserve(age_group_names.size());
        std::transform(age_group_names.begin(), age_group_names.end(), std::back_inserter(age_groups),
                       [&obj](auto&& age_name) {
                           return obj.expect_element(age_name, Tag<double>{});
                       });
        return apply(
            io,
            [](auto&& ag, auto&& sid, auto&& cid) {
                return PopulationDataEntry{
                    CustomIndexArray<double, AgeGroup>(AgeGroup(ag.size()), ag.begin(), ag.end()), sid, cid};
            },
            details::unpack_all(age_groups), state_id, county_id);
    }
};

namespace details
{
inline void get_rki_age_interpolation_coefficients(const std::vector<double>& age_ranges,
                                                   std::vector<std::vector<double>>& interpolation,
                                                   std::vector<bool>& carry_over)
{
    std::array<double, 6> param_ranges = {5., 10., 20., 25., 20., 20.};
    static_assert(param_ranges.size() == ConfirmedCasesDataEntry::age_group_names.size(),
                  "Number of RKI age groups does not match number of age ranges.");

    //counter for parameter age groups
    size_t counter = 0;

    //residual of param age groups
    double res = 0.0;
    for (size_t i = 0; i < age_ranges.size(); i++) {

        // if current param age group didn't fit into previous rki age group, transfer residual to current age group
        if (res < 0) {
            interpolation[i].push_back(std::min(-res / age_ranges[i], 1.0));
        }

        if (counter < param_ranges.size() - 1) {
            res += age_ranges[i];
            if (std::abs(res) < age_ranges[i]) {
                counter++;
            }
            // iterate over param age groups while there is still room in the current rki age group
            while (res > 0) {
                res -= param_ranges[counter];
                interpolation[i].push_back((param_ranges[counter] + std::min(res, 0.0)) / age_ranges[i]);
                if (res >= 0) {
                    counter++;
                }
            }
            if (res < 0) {
                carry_over.push_back(true);
            }
            else if (res == 0) {
                carry_over.push_back(false);
            }
        }
        // if last param age group is reached
        else {
            interpolation[i].push_back((age_ranges[i] + res) / age_ranges[i]);
            if (res < 0 || counter == 0) {
                carry_over.push_back(true);
            }
            else if (res == 0) {
                carry_over.push_back(false);
            }
            res = 0;
        }
    }
}

inline std::vector<PopulationDataEntry>
interpolate_to_rki_age_groups(const std::vector<PopulationDataEntry>& population_data)
{
    std::vector<double> age_ranges = {3., 3., 9., 3., 7., 5., 10., 10., 15., 10., 25.};
    std::vector<std::vector<double>> coefficients{age_ranges.size()};
    std::vector<bool> carry_over{};
    get_rki_age_interpolation_coefficients(age_ranges, coefficients, carry_over);

    std::vector<PopulationDataEntry> interpolated{population_data};
    for (auto region_entry_idx = size_t(0); region_entry_idx < population_data.size(); ++region_entry_idx) {
        interpolated[region_entry_idx].population =
            CustomIndexArray<double, AgeGroup>(AgeGroup(ConfirmedCasesDataEntry::age_group_names.size()), 0.0);

        size_t interpolated_age_idx = 0;
        for (size_t age_idx = 0; age_idx < coefficients.size(); age_idx++) {
            for (size_t coeff_idx = 0; coeff_idx < coefficients[age_idx].size(); coeff_idx++) {
                interpolated[region_entry_idx].population[AgeGroup(interpolated_age_idx)] +=
                    coefficients[age_idx][coeff_idx] * population_data[region_entry_idx].population[AgeGroup(age_idx)];
                if (coeff_idx < coefficients[age_idx].size() - 1 || !carry_over[age_idx]) {
                    interpolated_age_idx++;
                }
            }
        }
    }

    return interpolated;
}
} // namespace details

/**
 * Deserialize population data from a JSON value.
 * Age groups are interpolated to RKI age groups.
 * @param jsvalue JSON value that contains the population data.
 * @return list of population data.
 */
inline IOResult<std::vector<PopulationDataEntry>> deserialize_population_data(const Json::Value& jsvalue)
{
    BOOST_OUTCOME_TRY(population_data, deserialize_json(jsvalue, Tag<std::vector<PopulationDataEntry>>{}));
    return success(details::interpolate_to_rki_age_groups(population_data));
}

/**
 * Deserialize population data from a JSON file.
 * Age groups are interpolated to RKI age groups.
 * @param filename JSON file that contains the population data.
 * @return list of population data.
 */
inline IOResult<std::vector<PopulationDataEntry>> read_population_data(const std::string& filename)
{
    BOOST_OUTCOME_TRY(jsvalue, read_json(filename));
    return deserialize_population_data(jsvalue);
}

/**
 * @brief returns a vector with the ids of all German counties.
 * @param path directory to population data
 * @return list of county ids.
 */
IOResult<std::vector<int>> get_county_ids(const std::string& path);

template <typename TestNTrace, typename ContactPattern, class Model, class Parameters, class ReadFunction>
IOResult<void> set_nodes(const Parameters& params, Date start_date, Date end_date, const fs::path& data_dir,
                         Graph<Model, MigrationParameters>& params_graph, ReadFunction&& read_func,
                         const std::vector<double>& scaling_factor_inf, double scaling_factor_icu,
                         double tnt_capacity_factor, int num_days = 0, bool export_time_series = false)
{
    BOOST_OUTCOME_TRY(county_ids, get_county_ids((data_dir / "pydata" / "Germany").string()));
    std::vector<Model> counties(county_ids.size(), Model(int(size_t(params.get_num_groups()))));
    for (auto& county : counties) {
        county.parameters = params;
    }

    BOOST_OUTCOME_TRY(read_func(counties, start_date, county_ids, scaling_factor_inf, scaling_factor_icu,
                                (data_dir / "pydata" / "Germany").string(), num_days, export_time_series));

    for (size_t county_idx = 0; county_idx < counties.size(); ++county_idx) {

        auto tnt_capacity = counties[county_idx].populations.get_total() * tnt_capacity_factor;

        //local parameters
        auto& tnt_value = counties[county_idx].parameters.template get<TestNTrace>();
        tnt_value       = UncertainValue(0.5 * (1.2 * tnt_capacity + 0.8 * tnt_capacity));
        tnt_value.set_distribution(mio::ParameterDistributionUniform(0.8 * tnt_capacity, 1.2 * tnt_capacity));

        //holiday periods
        auto holiday_periods = regions::get_holidays(regions::get_state_id(regions::CountyId(county_ids[county_idx])),
                                                     start_date, end_date);
        auto& contacts       = counties[county_idx].parameters.template get<ContactPattern>();
        contacts.get_school_holidays() =
            std::vector<std::pair<mio::SimulationTime, mio::SimulationTime>>(holiday_periods.size());
        std::transform(
            holiday_periods.begin(), holiday_periods.end(), contacts.get_school_holidays().begin(), [=](auto& period) {
                return std::make_pair(mio::SimulationTime(mio::get_offset_in_days(period.first, start_date)),
                                      mio::SimulationTime(mio::get_offset_in_days(period.second, start_date)));
            });

        //uncertainty in populations
        for (auto i = mio::AgeGroup(0); i < params.get_num_groups(); i++) {
            for (auto j = Index<typename Model::Compartments>(0); j < Model::Compartments::Count; ++j) {
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
            for (auto age = AgeGroup(0); age < populations.template size<mio::AgeGroup>(); ++age) {
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

template <typename TestNTrace, typename ContactPattern, class ContactLocation, class InfectionState, class Model,
          class Parameters, class ReadFunction>
IOResult<Graph<Model, MigrationParameters>>
create_graph(Graph<Model, MigrationParameters>& params_graph, const Parameters& params, Date start_date, Date end_date,
             const fs::path& data_dir, ReadFunction&& read_func, const std::vector<double>& scaling_factor_inf,
             double scaling_factor_icu, double tnt_capacity_factor,
             std::initializer_list<InfectionState>& migrating_compartments, size_t contact_locations_size,
             int num_days = 0, bool export_time_series = false)
{
    const auto& set_node_function = set_nodes<TestNTrace, ContactPattern, Model, Parameters, ReadFunction>;
    const auto& set_edge_function = set_edges<ContactLocation, Model, InfectionState>;
    BOOST_OUTCOME_TRY(set_node_function(params, start_date, end_date, data_dir, params_graph, read_func,
                                        scaling_factor_inf, scaling_factor_icu, tnt_capacity_factor, num_days,
                                        export_time_series));
    BOOST_OUTCOME_TRY(set_edge_function(data_dir, params_graph, migrating_compartments, contact_locations_size));
    return success(params_graph);
}

/**
 * Represents an entry in a vaccination data file.
 */
class VaccinationDataEntry
{
public:
    static const std::array<const char*, 6> age_group_names;

    double num_vaccinations_completed;
    Date date;
    AgeGroup age_group;
    boost::optional<regions::StateId> state_id;
    boost::optional<regions::CountyId> county_id;

    template <class IoContext>
    static IOResult<VaccinationDataEntry> deserialize(IoContext& io)
    {
        auto obj                        = io.expect_object("VaccinationDataEntry");
        auto num_vaccinations_completed = obj.expect_element("Vacc_completed", Tag<double>{});
        auto date                       = obj.expect_element("Date", Tag<StringDate>{});
        auto age_group_str              = obj.expect_element("Age_RKI", Tag<std::string>{});
        auto state_id                   = obj.expect_optional("ID_County", Tag<regions::StateId>{});
        auto county_id                  = obj.expect_optional("ID_County", Tag<regions::CountyId>{});
        return mio::apply(
            io,
            [](auto nf, auto d, auto&& a_str, auto sid, auto cid) -> IOResult<VaccinationDataEntry> {
                auto it = std::find(age_group_names.begin(), age_group_names.end(), a_str);
                auto a  = AgeGroup(0);
                if (it != age_group_names.end()) {
                    a = AgeGroup(size_t(it - age_group_names.begin()));
                }
                else {
                    return failure(StatusCode::InvalidValue, "Invalid vaccination data age group.");
                }
                return success(VaccinationDataEntry{nf, d, a, sid, cid});
            },
            num_vaccinations_completed, date, age_group_str, state_id, county_id);
    }
};

/**
 * Deserialize vaccination data from a JSON value.
 * @param jsvalue JSON value that contains the vaccination data.
 * @return list of vaccination data.
 */
inline IOResult<std::vector<VaccinationDataEntry>> deserialize_vaccination_data(const Json::Value& jsvalue)
{
    return deserialize_json(jsvalue, Tag<std::vector<VaccinationDataEntry>>{});
}

/**
 * Read vaccination data from a JSON file.
 * @param filename JSON file that contains the vaccination data.
 * @return list of vaccination data.
 */
inline IOResult<std::vector<VaccinationDataEntry>> read_vaccination_data(const std::string& filename)
{
    BOOST_OUTCOME_TRY(jsvalue, read_json(filename));
    return deserialize_vaccination_data(jsvalue);
}

} // namespace mio

#endif //MEMILIO_HAS_JSONCPP

#endif //MEMILIO_IO_EPI_DATA_H