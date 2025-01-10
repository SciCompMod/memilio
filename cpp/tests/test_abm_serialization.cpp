/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: René Schmieding
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
#include "matchers.h"
#include "abm/config.h"
#include "abm/infection_state.h"
#include "abm/parameters.h"
#include "abm/test_type.h"
#include "abm/testing_strategy.h"
#include "abm/protection_event.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/io/json_serializer.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/uncertain_value.h"
#include "models/abm/location.h"
#include "models/abm/person.h"
#include "models/abm/trip_list.h"
#include "models/abm/model.h"
#include <bitset>
#include <cstddef>

#ifdef MEMILIO_HAS_JSONCPP

#include "json/value.h"

#include "gmock/gmock.h"
#include "gtest/gtest.h"

/**
 * @brief Test de- and serialization of an object without equality operator.
 *
 * Test that a json value x representing type T is equal to serialize(deserialize(x)).
 *
 * Assuming the (de)serialization functions' general behavior is independent of specific values of member variables,
 * i.e. the function does not contain conditionals (`if (t > 0)`), optionals (`add_optional`/`expect_optional`), etc.,
 * and assuming that deserialize is injective (meaning that two unequal instances of T do not have the same json
 * representation, which can happen e.g. if not all member variables are serialized),
 * this sufficiently tests that serialize and deserialize are inverse functions to each other.
 *
 * @tparam T The type to test.
 * @param reference_json A json value representing an instance of T.
 */
template <class T>
void test_json_serialization(const Json::Value& reference_json)
{
    // check that the json is deserializable (i.e. a valid representation)
    auto t_result = mio::deserialize_json(reference_json, mio::Tag<T>());
    ASSERT_THAT(print_wrap(t_result), IsSuccess());

    // check that the resulting type T is serializable
    auto json_result = mio::serialize_json(t_result.value());
    ASSERT_THAT(print_wrap(json_result), IsSuccess());

    EXPECT_THAT(json_result.value(), JsonEqual(reference_json));
}

TEST(TestAbmSerialization, Trip)
{
    // See test_json_serialization for info on this test.

    unsigned i = 1; // counter s.t. members have different values

    Json::Value reference_json;
    reference_json["person_id"]       = Json::UInt(i++);
    reference_json["time"]["seconds"] = Json::Int(i++);
    reference_json["destination"]     = Json::UInt(i++);
    reference_json["origin"]          = Json::UInt(i++);

    test_json_serialization<mio::abm::Trip>(reference_json);
}

TEST(TestAbmSerialization, ProtectionEvent)
{
    // See test_json_serialization for info on this test.

    Json::Value reference_json;
    reference_json["type"]            = Json::UInt(1);
    reference_json["time"]["seconds"] = Json::Int(2);

    test_json_serialization<mio::abm::ProtectionEvent>(reference_json);
}

TEST(TestAbmSerialization, Infection)
{
    // See test_json_serialization for info on this test.

    unsigned i = 1; // counter s.t. members have different values

    Json::Value viral_load;
    viral_load["decline"]               = Json::Value((double)i++);
    viral_load["end_date"]["seconds"]   = Json::Int(i++);
    viral_load["incline"]               = Json::Value((double)i++);
    viral_load["peak"]                  = Json::Value((double)i++);
    viral_load["start_date"]["seconds"] = Json::Int(i++);

    Json::Value reference_json;
    reference_json["infection_course"] = Json::Value(Json::arrayValue);
    reference_json["virus_variant"]    = Json::UInt(0);
    reference_json["viral_load"]       = viral_load;
    reference_json["log_norm_alpha"]   = Json::Value((double)i++);
    reference_json["log_norm_beta"]    = Json::Value((double)i++);
    reference_json["detected"]         = Json::Value((bool)0);

    test_json_serialization<mio::abm::Infection>(reference_json);
}

TEST(TestAbmSerialization, TestingScheme)
{
    // See test_json_serialization for info on this test.

    unsigned i = 1; // counter s.t. members have different values

    Json::Value testing_criteria;
    std::bitset<mio::abm::MAX_NUM_AGE_GROUPS> ages_bits{}; // initialize to false
    ages_bits[i++]           = true;
    testing_criteria["ages"] = mio::serialize_json(ages_bits).value();
    std::bitset<(size_t)mio::abm::InfectionState::Count> inf_st_bits{}; // initialize to false
    inf_st_bits[i++]                     = true;
    testing_criteria["infection_states"] = mio::serialize_json(inf_st_bits).value();

    Json::Value test_parameters;
    test_parameters["sensitivity"]              = mio::serialize_json(mio::UncertainValue<double>{(double)i++}).value();
    test_parameters["specificity"]              = mio::serialize_json(mio::UncertainValue<double>{(double)i++}).value();
    test_parameters["required_time"]["seconds"] = Json::Int(i++);
    test_parameters["test_type"]                = Json::UInt(0);

    Json::Value reference_json;
    reference_json["criteria"]                   = testing_criteria;
    reference_json["validity_period"]["seconds"] = Json::Int(i++);
    reference_json["start_date"]["seconds"]      = Json::Int(i++);
    reference_json["end_date"]["seconds"]        = Json::Int(i++);
    reference_json["test_params"]                = test_parameters;
    reference_json["probability"]                = Json::Value((double)i++);
    reference_json["is_active"]                  = Json::Value((bool)0);

    test_json_serialization<mio::abm::TestingScheme>(reference_json);
}

TEST(TestAbmSerialization, TestingStrategy)
{
    // See test_json_serialization for info on this test.

    unsigned i = 1; // counter s.t. members have different values

    Json::Value local_strategy;
    local_strategy["id"]      = Json::UInt(i++);
    local_strategy["schemes"] = Json::Value(Json::arrayValue);
    local_strategy["type"]    = Json::UInt(i++);

    Json::Value reference_json;
    reference_json["schemes"][0] = local_strategy;

    test_json_serialization<mio::abm::TestingStrategy>(reference_json);
}

TEST(TestAbmSerialization, TestResult)
{
    // See test_json_serialization for info on this test.

    Json::Value reference_json;
    reference_json["result"]                     = Json::Value(false);
    reference_json["time_of_testing"]["seconds"] = Json::Int(1);

    test_json_serialization<mio::abm::TestResult>(reference_json);
}

TEST(TestAbmSerialization, Person)
{
    // See test_json_serialization for info on this test.

    auto json_uint_array = [](std::vector<uint32_t> values) {
        return mio::serialize_json(values).value();
    };
    auto json_double_array = [](std::vector<double> values) {
        return mio::serialize_json(values).value();
    };

    unsigned i = 1; // counter s.t. members have different values

    Json::Value reference_json;
    reference_json["age_group"]           = Json::UInt(i++);
    reference_json["assigned_locations"]  = json_uint_array({i++, i++, i++, i++, i++, i++, i++, i++, i++, i++, i++});
    reference_json["cells"]               = json_uint_array({i++});
    reference_json["compliance"]          = json_double_array({(double)i++, (double)i++, (double)i++});
    reference_json["id"]                  = Json::UInt(i++);
    reference_json["infections"]          = Json::Value(Json::arrayValue);
    reference_json["last_transport_mode"] = Json::UInt(i++);
    reference_json["location"]            = Json::UInt(i++);
    reference_json["location_type"]       = Json::UInt(0);
    reference_json["mask"]["mask_type"]   = Json::UInt(0);
    reference_json["mask"]["time_first_used"]["seconds"] = Json::Int(i++);
    reference_json["home_isolation_start"]["seconds"]    = Json::Int(i++);
    reference_json["rnd_go_to_school_hour"]              = Json::Value((double)i++);
    reference_json["rnd_go_to_work_hour"]                = Json::Value((double)i++);
    reference_json["rnd_schoolgroup"]                    = Json::Value((double)i++);
    reference_json["rnd_workgroup"]                      = Json::Value((double)i++);
    reference_json["rng_counter"]                        = Json::UInt(i++);
    reference_json["test_results"] =
        mio::serialize_json(mio::CustomIndexArray<mio::abm::TestResult, mio::abm::TestType>{}).value();
    reference_json["time_at_location"]["seconds"] = Json::Int(i++);
    reference_json["vaccinations"]                = Json::Value(Json::arrayValue);

    test_json_serialization<mio::abm::Person>(reference_json);
}

TEST(TestAbmSerialization, Location)
{
    // See test_json_serialization for info on this test.

    unsigned i = 1; // counter s.t. members have different values

    Json::Value contact_rates = mio::serialize_json(mio::abm::ContactRates::get_default(i++)).value();

    Json::Value reference_json;
    reference_json["cells"][0]["capacity"]["persons"]                   = Json::UInt(i++);
    reference_json["cells"][0]["capacity"]["volume"]                    = Json::UInt(i++);
    reference_json["geographical_location"]["latitude"]                 = Json::Value((double)i++);
    reference_json["geographical_location"]["longitude"]                = Json::Value((double)i++);
    reference_json["id"]                                                = Json::UInt(i++);
    reference_json["parameters"]["ContactRates"]                        = contact_rates;
    reference_json["parameters"]["MaximumContacts"]                     = Json::Value((double)i++);
    reference_json["parameters"]["UseLocationCapacityForTransmissions"] = Json::Value(false);
    reference_json["required_mask"]                                     = Json::UInt(0);
    reference_json["type"]                                              = Json::UInt(0);

    test_json_serialization<mio::abm::Location>(reference_json);
}

// TEST(TestAbmSerialization, Model)
// {
//     // See test_json_serialization for info on this test.

//     auto json_uint_array = [](std::vector<uint32_t> values) {
//         return mio::serialize_json(values).value();
//     };

//     unsigned i = 1; // counter s.t. members have different values

//     Json::Value abm_parameters = mio::serialize_json(mio::abm::Parameters(i++)).value();

//     Json::Value reference_json;
//     reference_json["cemetery_id"]                 = Json::UInt(i++);
//     reference_json["location_types"]              = Json::UInt(i++);
//     reference_json["locations"]                   = Json::Value(Json::arrayValue);
//     reference_json["parameters"]                  = abm_parameters;
//     reference_json["persons"]                     = Json::Value(Json::arrayValue);
//     reference_json["rng"]["counter"]              = Json::UInt(i++);
//     reference_json["rng"]["key"]                  = Json::UInt(i++);
//     reference_json["rng"]["seeds"]                = json_uint_array({i++, i++, i++, i++, i++, i++});
//     reference_json["testing_strategy"]["schemes"] = Json::Value(Json::arrayValue);
//     reference_json["trip_list"]["index"]          = Json::UInt(i++);
//     reference_json["trip_list"]["trips_weekday"]  = Json::Value(Json::arrayValue);
//     reference_json["trip_list"]["trips_weekend"]  = Json::Value(Json::arrayValue);
//     reference_json["use_mobility_rules"]          = Json::Value(false);

//     test_json_serialization<mio::abm::Model>(reference_json);
// }

#endif
