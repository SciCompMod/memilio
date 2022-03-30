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
#include "memilio/epidemiology/regions.h"
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
                BOOST_OUTCOME_TRY(date, parse_date(str_));
                return success(date);
            },
            str);
    }
};

/**
 * Age group that serializes into an RKI age string.
 */
class StringRkiAgeGroup : public AgeGroup
{
public:
    using AgeGroup::AgeGroup;

    StringRkiAgeGroup(const AgeGroup& other)
        : AgeGroup(other)
    {
    }

    static const std::array<const char*, 6> age_group_names;

    template <class IoContext>
    static IOResult<StringRkiAgeGroup> deserialize(IoContext& io)
    {
        auto str = mio::deserialize(io, Tag<std::string>{});
        return apply(
            io,
            [](auto&& str_) -> IOResult<StringRkiAgeGroup> {
                auto it = std::find(age_group_names.begin(), age_group_names.end(), str_);
                if (it != age_group_names.end()) {
                    return success(size_t(it - age_group_names.begin()));
                } else if (str_ == "unknown") {
                    return success(size_t(age_group_names.size()));
                }
                return failure(StatusCode::InvalidValue, "Invalid age group.");
            },
            str);
    }
};

/**
 * Represents the entries of an RKI data file.
 * Number of confirmed, recovered and deceased in a region on a specific date.
 * Region can be a county, a state, or a country. If it is a country, both
 * state_id and county_id will be empty.
 */
class RkiEntry
{
public:
    double num_confirmed;
    double num_recovered;
    double num_deceased;
    Date date;
    AgeGroup age_group;
    boost::optional<regions::de::StateId> state_id;
    boost::optional<regions::de::CountyId> county_id;

    template <class IoContext>
    static IOResult<RkiEntry> deserialize(IoContext& io)
    {
        auto obj           = io.expect_object("RkiEntry");
        auto num_confirmed = obj.expect_element("Confirmed", Tag<double>{});
        auto num_deceased  = obj.expect_element("Deaths", Tag<double>{});
        auto num_recovered = obj.expect_element("Recovered", Tag<double>{});
        auto date          = obj.expect_element("Date", Tag<StringDate>{});
        auto age_group     = obj.expect_element("Age_RKI", Tag<StringRkiAgeGroup>{});
        auto state_id      = obj.expect_optional("ID_State", Tag<regions::de::StateId>{});
        auto county_id     = obj.expect_optional("ID_County", Tag<regions::de::CountyId>{});
        return apply(
            io,
            [](auto&& nc, auto&& nr, auto&& nd, auto&& d, auto&& a, auto&& sid, auto&& cid) {
                return RkiEntry{nc, nr, nd, d, a, sid, cid};
            },
            num_confirmed, num_recovered, num_deceased, date, age_group, state_id, county_id);
    }
};

/**
 * Read list of RkiEntry from json.
 * @param jsvalue json value, must be an array of objects, objects must match RkiEntry.
 * @return list of entries; entries of unknown age group are omitted.
 */
inline IOResult<std::vector<RkiEntry>> deserialize_rki_data(const Json::Value& jsvalue)
{
    BOOST_OUTCOME_TRY(rki_data, deserialize_json(jsvalue, Tag<std::vector<RkiEntry>>{}));
    //filter unknown age group
    auto it = std::remove_if(rki_data.begin(), rki_data.end(), [](auto&& rki_entry) {
        return rki_entry.age_group >= AgeGroup(StringRkiAgeGroup::age_group_names.size());
    });
    rki_data.erase(it, rki_data.end());
    return success(std::move(rki_data));
}

/**
 * Read list of RkiEntry from a json file.
 * @param filename name of the json file. File content must be an array of objects, objects must match RkiEntry.
 * @return list of entries; entries of unknown age group are omitted.
 */
inline IOResult<std::vector<RkiEntry>> read_rki_data(const std::string& filename)
{
    BOOST_OUTCOME_TRY(jsvalue, read_json(filename));
    return deserialize_rki_data(jsvalue);
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
    boost::optional<regions::de::StateId> state_id;
    boost::optional<regions::de::CountyId> county_id;    

    template <class IoContext>
    static IOResult<DiviEntry> deserialize(IoContext& io)
    {
        auto obj       = io.expect_object("RkiEntry");
        auto num_icu   = obj.expect_element("ICU", Tag<double>{});
        auto date      = obj.expect_element("Date", Tag<StringDate>{});
        auto state_id  = obj.expect_optional("ID_State", Tag<regions::de::StateId>{});
        auto county_id = obj.expect_optional("ID_County", Tag<regions::de::CountyId>{});
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
inline IOResult<std::vector<DiviEntry>> deserialize_divi_data(const Json::Value& jsvalue) {
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

namespace details {
    //check all results in a vector and unpack each
    template<class T>
    IOResult<std::vector<T>> unpack_all(const std::vector<IOResult<T>>& v) {
        std::vector<T> w;
        w.reserve(v.size());
        for (auto&& r : v) {
            BOOST_OUTCOME_TRY(t, r);
            w.push_back(t);
        }
        return success(w);
    }
}

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
    boost::optional<regions::de::StateId> state_id;
    boost::optional<regions::de::CountyId> county_id;    
    
    template <class IoContext>
    static IOResult<PopulationDataEntry> deserialize(IoContext& io)
    {
        auto obj       = io.expect_object("RkiEntry");
        auto state_id  = obj.expect_optional("ID_State", Tag<regions::de::StateId>{});
        auto county_id = obj.expect_optional("ID_County", Tag<regions::de::CountyId>{});
        std::vector<IOResult<double>> age_groups;
        age_groups.reserve(age_group_names.size());
        std::transform(age_group_names.begin(), age_group_names.end(), std::back_inserter(age_groups), [&obj](auto&& age_name) {
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
        static_assert(param_ranges.size() == StringRkiAgeGroup::age_group_names.size(), "RKI age group mismatch.");

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
        std::vector<double> age_ranges     = {3., 3., 9., 3., 7., 5., 10., 10., 15., 10., 25.};
        std::vector<std::vector<double>> coefficients{ age_ranges.size() };
        std::vector<bool> carry_over{};
        get_rki_age_interpolation_coefficients(age_ranges, coefficients, carry_over);

        std::vector<PopulationDataEntry> interpolated{population_data};
        for (auto entry_idx = size_t(0); entry_idx < population_data.size(); ++entry_idx) {
            interpolated[entry_idx].population =
                CustomIndexArray<double, AgeGroup>(AgeGroup(StringRkiAgeGroup::age_group_names.size()), 0.0);

            size_t interpolated_age_idx = 0;
            for (size_t age_idx = 0; age_idx < coefficients.size(); age_idx++) {
                for (size_t coeff_idx = 0; coeff_idx < coefficients[age_idx].size(); coeff_idx++) {
                    interpolated[entry_idx].population[AgeGroup(interpolated_age_idx)] +=
                        coefficients[age_idx][coeff_idx] * population_data[entry_idx].population[AgeGroup(age_idx)];
                    if (coeff_idx < coefficients[age_idx].size() - 1 || !carry_over[age_idx]) {
                        interpolated_age_idx++;
                    }
                }
            }
        }

        return interpolated;
    }
}

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
 * @brief returns a vector with the ids of all german counties
 * @param path directory to population data
 * @return list of county ids.
 */
IOResult<std::vector<int>> get_county_ids(const std::string& path);

/**
 * Represents an age group string in vaccination data.
 */
class StringVaccinationDataAgeGroup : public AgeGroup
{
public:
    using AgeGroup::AgeGroup;

    static const std::array<const char*, 6> age_group_names;   

    template <class IoContext>
    static IOResult<StringVaccinationDataAgeGroup> deserialize(IoContext& io)
    {
        auto str = mio::deserialize(io, Tag<std::string>{});
        return apply(
            io,
            [](auto&& str_) -> IOResult<StringVaccinationDataAgeGroup> {
                auto it = std::find(age_group_names.begin(), age_group_names.end(), str_);
                if (it != age_group_names.end()) {
                    return success(size_t(it - age_group_names.begin()));
                }
                return failure(StatusCode::InvalidValue, "Invalid age group.");
            },
            str);
    }
};

/**
 * Represents an entry in a vaccination data file.
 */
class VaccinationDataEntry
{
public:
    double num_full;
    Date date;
    AgeGroup age_group;
    boost::optional<regions::de::StateId> state_id;
    boost::optional<regions::de::CountyId> county_id;

    template<class IoContext>
    static IOResult<VaccinationDataEntry> deserialize(IoContext& io)
    {
        auto obj = io.expect_object("VaccinationDataEntry");
        auto num_full = obj.expect_element("Vacc_completed", Tag<double>{});
        auto date = obj.expect_element("Date", Tag<StringDate>{});
        auto age_group = obj.expect_element("Age_RKI", Tag<StringVaccinationDataAgeGroup>{});
        auto state_id = obj.expect_optional("ID_County", Tag<regions::de::StateId>{});
        auto county_id = obj.expect_optional("ID_County", Tag<regions::de::CountyId>{});
        return mio::apply(io, [](auto nf, auto d, auto a, auto sid, auto cid) {
            return VaccinationDataEntry{nf, d, a, sid, cid};
        }, num_full, date, age_group, state_id, county_id);
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