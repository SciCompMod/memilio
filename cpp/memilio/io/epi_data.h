/*
* Copyright (C) 2020-2025 MEmilio
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
                BOOST_OUTCOME_TRY(auto&& date, parse_date(str_));
                return success(date);
            },
            str);
    }
};

/**
 * Represents the entries of a confirmed cases data file, e.g., from RKI.
 * Number of confirmed, recovered and deceased in Germany on a specific date.
 * ConfirmedCasesNoAgeEntry is a simplified version of ConfirmedCasesDataEntry without agegroups.
 */
class ConfirmedCasesNoAgeEntry
{
public:
    double num_confirmed;
    double num_recovered;
    double num_deaths;
    Date date;

    template <class IOContext>
    static IOResult<ConfirmedCasesNoAgeEntry> deserialize(IOContext& io)
    {
        auto obj           = io.expect_object("ConfirmedCasesNoAgeEntry");
        auto num_confirmed = obj.expect_element("Confirmed", Tag<double>{});
        auto num_recovered = obj.expect_element("Recovered", Tag<double>{});
        auto num_deaths    = obj.expect_element("Deaths", Tag<double>{});
        auto date          = obj.expect_element("Date", Tag<StringDate>{});
        return apply(
            io,
            [](auto&& nc, auto&& nr, auto&& nd, auto&& d) {
                return ConfirmedCasesNoAgeEntry{nc, nr, nd, d};
            },
            num_confirmed, num_recovered, num_deaths, date);
    }
};

/**
 * Deserialize a list of ConfirmedCasesNoAgeEntry from json.
 * @param jsvalue json value, must be an array of objects, objects must match ConfirmedCasesNoAgeEntry.
 *  Accordingly, there should be one entry per date with the values for Confirmed, Recovered and Deaths. In addition, no age groups should be specified.
 * @return List of ConfirmedCasesNoAgeEntry.
 */
inline IOResult<std::vector<ConfirmedCasesNoAgeEntry>> deserialize_confirmed_cases_noage(const Json::Value& jsvalue)
{
    return deserialize_json(jsvalue, Tag<std::vector<ConfirmedCasesNoAgeEntry>>{});
}

/**
 * Read list of ConfirmedCasesNoAgeEntry from a json file.
 * @param filename name of the json file. File content must be an array of objects, objects must match ConfirmedCasesNoAgeEntry.
 *  Accordingly, there should be one entry per date with the values for Confirmed, Recovered and Deaths. In addition, no age groups should be specified.
 * @return List of ConfirmedCasesNoAgeEntry.
 */
inline IOResult<std::vector<ConfirmedCasesNoAgeEntry>> read_confirmed_cases_noage(const std::string& filename)
{
    return read_json(filename, Tag<std::vector<ConfirmedCasesNoAgeEntry>>{});
}

/**
 * Represents the entries of a confirmed cases data file, e.g., from RKI.
 * Number of confirmed, recovered and deceased in a region on a specific date.
 * Region can be a county, a state, or a country. If it is a country, both
 * state_id and county_id will be empty.
 */
class ConfirmedCasesDataEntry
{
public:
    static std::vector<const char*> age_group_names;

    double num_confirmed;
    double num_recovered;
    double num_deaths;
    Date date;
    AgeGroup age_group;
    boost::optional<regions::StateId> state_id;
    boost::optional<regions::CountyId> county_id;
    boost::optional<regions::DistrictId> district_id;

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
        auto district_id   = obj.expect_optional("ID_District", Tag<regions::DistrictId>{});
        return apply(
            io,
            [](auto&& nc, auto&& nr, auto&& nd, auto&& d, auto&& a_str, auto&& sid, auto&& cid,
               auto&& did) -> IOResult<ConfirmedCasesDataEntry> {
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
                return success(ConfirmedCasesDataEntry{nc, nr, nd, d, a, sid, cid, did});
            },
            num_confirmed, num_recovered, num_deaths, date, age_group_str, state_id, county_id, district_id);
    }
};

/**
 * Read list of ConfirmedCasesDataEntry from json.
 * @param jsvalue json value, must be an array of objects, objects must match ConfirmedCasesDataEntry.
 * @return list of entries; entries of unknown age group are omitted.
 */
inline IOResult<std::vector<ConfirmedCasesDataEntry>> deserialize_confirmed_cases_data(const Json::Value& jsvalue)
{
    BOOST_OUTCOME_TRY(auto&& cases_data, deserialize_json(jsvalue, Tag<std::vector<ConfirmedCasesDataEntry>>{}));
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
    BOOST_OUTCOME_TRY(auto&& jsvalue, read_json(filename));
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
    boost::optional<regions::DistrictId> district_id;

    template <class IoContext>
    static IOResult<DiviEntry> deserialize(IoContext& io)
    {
        auto obj         = io.expect_object("DiviEntry");
        auto num_icu     = obj.expect_element("ICU", Tag<double>{});
        auto date        = obj.expect_element("Date", Tag<StringDate>{});
        auto state_id    = obj.expect_optional("ID_State", Tag<regions::StateId>{});
        auto county_id   = obj.expect_optional("ID_County", Tag<regions::CountyId>{});
        auto district_id = obj.expect_optional("ID_District", Tag<regions::DistrictId>{});
        return apply(
            io,
            [](auto&& ni, auto&& d, auto&& sid, auto&& cid, auto&& did) {
                return DiviEntry{ni, d, sid, cid, did};
            },
            num_icu, date, state_id, county_id, district_id);
    }
};

/**
 * @brief Checks if DIVI data is available for a given date.
 * @param date The date to check.
 * @return True if date is within 2020-04-23 and 2024-07-21, false otherwise.
 */
inline bool is_divi_data_available(const Date& date)
{
    static const Date divi_data_start(2020, 4, 23);
    static const Date divi_data_end(2024, 7, 21);
    return date >= divi_data_start && date <= divi_data_end;
}

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
        BOOST_OUTCOME_TRY(auto&& t, r);
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
    static std::vector<const char*> age_group_names;

    CustomIndexArray<double, AgeGroup> population;
    boost::optional<regions::StateId> state_id;
    boost::optional<regions::CountyId> county_id;
    boost::optional<regions::DistrictId> district_id;

    template <class IoContext>
    static IOResult<PopulationDataEntry> deserialize(IoContext& io)
    {
        auto obj         = io.expect_object("PopulationDataEntry");
        auto state_id    = obj.expect_optional("ID_State", Tag<regions::StateId>{});
        auto county_id   = obj.expect_optional("ID_County", Tag<regions::CountyId>{});
        auto district_id = obj.expect_optional("ID_District", Tag<regions::DistrictId>{});
        std::vector<IOResult<double>> age_groups;
        age_groups.reserve(age_group_names.size());
        std::transform(age_group_names.begin(), age_group_names.end(), std::back_inserter(age_groups),
                       [&obj](auto&& age_name) {
                           return obj.expect_element(age_name, Tag<double>{});
                       });
        return apply(
            io,
            [](auto&& ag, auto&& sid, auto&& cid, auto&& did) {
                return PopulationDataEntry{
                    CustomIndexArray<double, AgeGroup>(AgeGroup(ag.size()), ag.begin(), ag.end()), sid, cid, did};
            },
            details::unpack_all(age_groups), state_id, county_id, district_id);
    }
};

namespace details
{
inline void get_rki_age_interpolation_coefficients(const std::vector<double>& age_ranges,
                                                   std::vector<std::vector<double>>& interpolation,
                                                   std::vector<bool>& carry_over)
{
    std::array<double, 6> param_ranges = {5., 10., 20., 25., 20., 20.};
    assert(param_ranges.size() == ConfirmedCasesDataEntry::age_group_names.size() &&
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
 * @param rki_age_groups Specifies whether population data should be interpolated to rki age groups.
 * @return list of population data.
 */
inline IOResult<std::vector<PopulationDataEntry>> deserialize_population_data(const Json::Value& jsvalue,
                                                                              bool rki_age_groups = true)
{
    BOOST_OUTCOME_TRY(auto&& population_data, deserialize_json(jsvalue, Tag<std::vector<PopulationDataEntry>>{}));
    if (rki_age_groups) {
        return success(details::interpolate_to_rki_age_groups(population_data));
    }
    else {
        return success(population_data);
    }
}

/**
 * Deserialize population data from a JSON file.
 * Age groups are interpolated to RKI age groups.
 * @param filename JSON file that contains the population data.
 * @return list of population data.
 */
inline IOResult<std::vector<PopulationDataEntry>> read_population_data(const std::string& filename,
                                                                       bool rki_age_group = true)
{
    BOOST_OUTCOME_TRY(auto&& jsvalue, read_json(filename));
    return deserialize_population_data(jsvalue, rki_age_group);
}

/**
 * @brief Sets the age groups' names for the ConfirmedCasesDataEntry%s.
 * @param[in] names age group names
*/
IOResult<void> set_confirmed_cases_age_group_names(std::vector<const char*> names);

/**
 * @brief Sets the age groups' names for the PopulationDataEntry%s.
 * @param[in] names age group names
*/
IOResult<void> set_population_data_age_group_names(std::vector<const char*> names);

/**
 * @brief Sets the age groups' names for the VaccinationDataEntry%s.
 * @param[in] names age group names
*/
IOResult<void> set_vaccination_data_age_group_names(std::vector<const char*> names);

/**
 * @brief returns a vector with the ids of all nodes.
 * @param[in] path directory to population data
 * @param[in] is_node_for_county boolean specifying whether the nodes should be counties or districts
 * @param[in] rki_age_groups boolean specifying whether population data should be interpolated to rki age groups.
 * @return list of node ids.
 */
IOResult<std::vector<int>> get_node_ids(const std::string& path, bool is_node_for_county, bool rki_age_groups = true);

/**
 * Represents an entry in a vaccination data file.
 */
class VaccinationDataEntry
{
public:
    static std::vector<const char*> age_group_names;

    double num_vaccinations_partial, num_vaccinations_completed, num_vaccinations_refreshed_first,
        num_vaccinations_refreshed_additional;
    Date date;
    AgeGroup age_group;
    boost::optional<regions::StateId> state_id;
    boost::optional<regions::CountyId> county_id;
    boost::optional<regions::DistrictId> district_id;

    template <class IoContext>
    static IOResult<VaccinationDataEntry> deserialize(IoContext& io)
    {
        auto obj                                   = io.expect_object("VaccinationDataEntry");
        auto num_vaccinations_partial              = obj.expect_element("Vacc_partially", Tag<double>{});
        auto num_vaccinations_completed            = obj.expect_element("Vacc_completed", Tag<double>{});
        auto num_vaccinations_refreshed_first      = obj.expect_optional("Vacc_refreshed", Tag<double>{});
        auto num_vaccinations_refreshed_additional = obj.expect_optional("Vacc_refreshed_2", Tag<double>{});
        auto date                                  = obj.expect_element("Date", Tag<StringDate>{});
        auto age_group_str                         = obj.expect_element("Age_RKI", Tag<std::string>{});
        auto state_id                              = obj.expect_optional("ID_County", Tag<regions::StateId>{});
        auto county_id                             = obj.expect_optional("ID_County", Tag<regions::CountyId>{});
        auto district_id                           = obj.expect_optional("ID_District", Tag<regions::DistrictId>{});
        return mio::apply(
            io,
            [](auto np, auto nc, auto n_refreshed_1, auto n_refreshed_2, auto d, auto&& a_str, auto sid, auto cid,
               auto did) -> IOResult<VaccinationDataEntry> {
                auto it = std::find(age_group_names.begin(), age_group_names.end(), a_str);
                auto a  = AgeGroup(0);
                if (it != age_group_names.end()) {
                    a = AgeGroup(size_t(it - age_group_names.begin()));
                }
                else {
                    return failure(StatusCode::InvalidValue, "Invalid vaccination data age group.");
                }
                // Optional values are 0 if they do not exist.
                auto n_refreshed_1_value = n_refreshed_1.value_or(0.0);
                auto n_refreshed_2_value = n_refreshed_2.value_or(0.0);
                return success(
                    VaccinationDataEntry{np, nc, n_refreshed_1_value, n_refreshed_2_value, d, a, sid, cid, did});
            },
            num_vaccinations_partial, num_vaccinations_completed, num_vaccinations_refreshed_first,
            num_vaccinations_refreshed_additional, date, age_group_str, state_id, county_id, district_id);
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
    BOOST_OUTCOME_TRY(auto&& jsvalue, read_json(filename));
    return deserialize_vaccination_data(jsvalue);
}

} // namespace mio

#endif //MEMILIO_HAS_JSONCPP

#endif //MEMILIO_IO_EPI_DATA_H
