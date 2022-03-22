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
#ifndef MEMILIO_IO_RKI_DATA_H
#define MEMILIO_IO_RKI_DATA_H

#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/regions.h"
#include "memilio/io/json_serializer.h"
#include "memilio/utils/date.h"

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

    static const std::array<const char*, 7> age_names;

    template <class IoContext>
    static IOResult<StringRkiAgeGroup> deserialize(IoContext& io)
    {
        auto str = mio::deserialize(io, Tag<std::string>{});
        return apply(
            io,
            [](auto&& str_) -> IOResult<StringRkiAgeGroup> {
                auto it = std::find(age_names.begin(), age_names.end(), str_);
                if (it != age_names.end()) {
                    return success(size_t(it - age_names.begin()));
                }
                return failure(StatusCode::InvalidValue, "Invalid age group.");
            },
            str);
    }
};

/**
 * Represents the entries of an RKI data file.
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
        auto num_deceased  = obj.expect_element("Dead", Tag<double>{});
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
inline IOResult<std::vector<RkiEntry>> read_rki_data(const Json::Value& jsvalue)
{
    BOOST_OUTCOME_TRY(rki_data, deserialize_json(jsvalue, Tag<std::vector<RkiEntry>>{}));
    //filter unknown age group
    auto it = std::remove_if(rki_data.begin(), rki_data.end(), [](auto&& rki_entry) {
        return rki_entry.age_group >= AgeGroup(StringRkiAgeGroup::age_names.size() - 1);
    });
    rki_data.erase(it, rki_data.end());
    return rki_data;
}

/**
 * Read list of RkiEntry from a json file.
 * @param filename name of the json file. File content must be an array of objects, objects must match RkiEntry.
 * @return list of entries; entries of unknown age group are omitted.
 */
inline IOResult<std::vector<RkiEntry>> read_rki_data(const std::string& filename)
{
    BOOST_OUTCOME_TRY(jsvalue, read_json(filename));
    return read_rki_data(jsvalue);
}

} // namespace mio

#endif //MEMILIO_IO_RKI_DATA_H