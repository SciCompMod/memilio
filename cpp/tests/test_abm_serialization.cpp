#include "abm/parameters.h"
#include "abm/testing_strategy.h"
#include "abm/vaccine.h"
#include "matchers.h"
#include "memilio/io/json_serializer.h"
#include "memilio/utils/uncertain_value.h"
#include "models/abm/location.h"
#include "models/abm/person.h"
#include "models/abm/time.h"
#include "models/abm/trip_list.h"
#include "models/abm/world.h"

#ifdef MEMILIO_HAS_JSONCPP

void test_equal_json_representation(const Json::Value& test_json, const Json::Value& reference_json)
{
    // write the resulting json value and the reference value to string to compare their representations.
    Json::StreamWriterBuilder swb;
    swb["indentation"] = " ";
    auto js_writer     = std::unique_ptr<Json::StreamWriter>(swb.newStreamWriter());
    std::stringstream test_str, reference_str;
    js_writer->write(reference_json, &reference_str);
    js_writer->write(test_json, &test_str);
    // we compare strings here, as e.g. Json::Int(5) != Json::Uint(5), but their json representation is the same
    EXPECT_EQ(test_str.str(), reference_str.str());
}

/**
 * @brief Test de- and serialization of an object by comparing its json representation.
 *
 * Test that a json value x representing type T is equal to serialize(deserialize(x)) w.r.t json representation.
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
void test_json_serialization_by_representation(const Json::Value& reference_json)
{
    // check that the json is deserializable (i.e. a valid representation)
    auto t_result = mio::deserialize_json(reference_json, mio::Tag<T>());
    ASSERT_THAT(print_wrap(t_result), IsSuccess());

    // check that the resulting type T is serializable
    auto json_result = mio::serialize_json(t_result.value());
    ASSERT_TRUE(json_result);

    test_equal_json_representation(json_result.value(), reference_json);
}

/**
 * @brief Test de- and serialization of an object by comparing its json representation and using its equality operator.
 *
 * First, test that serializing the reference_object is equal to the reference_json (w.r.t. their representation),
 * and that deserializing the reference_json results in an object equal to the reference_object.
 * Then, repeat this step using its own results as arguments to (de)serialize, to check that serialization and
 * deserialization are inverse functions to each other.
 *
 * @tparam T The type to test.
 * @param reference_object An instance of T.
 * @param reference_json A json value representing reference_object.
 */
template <class T>
void test_json_serialization_full(const T& reference_object, const Json::Value& reference_json)
{
    // check that the reference type T is serializable
    auto json_result = mio::serialize_json(reference_object);
    ASSERT_TRUE(json_result);

    // check that the reference json is deserializable
    auto t_result = mio::deserialize_json(reference_json, mio::Tag<T>());
    ASSERT_THAT(print_wrap(t_result), IsSuccess());

    // compare both results with other reference values
    EXPECT_EQ(t_result.value(), reference_object);
    test_equal_json_representation(json_result.value(), reference_json);

    // do the same once more using the results from above
    auto json_result_2 = mio::serialize_json(t_result.value());
    ASSERT_TRUE(json_result_2);
    auto t_result_2 = mio::deserialize_json(json_result.value(), mio::Tag<T>());
    ASSERT_THAT(print_wrap(t_result_2), IsSuccess());

    EXPECT_EQ(t_result_2.value(), reference_object);
    test_equal_json_representation(json_result_2.value(), reference_json);
}

TEST(TestAbmSerialization, Trip)
{
    // Test that a json value x is equal to serialize(deserialize(x)) w.r.t json representation.
    // See test_json_serialization_by_representation for more detail.

    mio::abm::Trip trip(1, mio::abm::TimePoint(0) + mio::abm::hours(2), 3, 4);

    Json::Value reference_json; // aka x
    reference_json["person_id"]   = Json::UInt(1);
    reference_json["time"]        = Json::Int(mio::abm::hours(2).seconds());
    reference_json["destination"] = Json::UInt(3);
    reference_json["origin"]      = Json::UInt(4);

    test_json_serialization_full(trip, reference_json);
}

TEST(TestAbmSerialization, Vaccination)
{
    // Test that a json value x is equal to serialize(deserialize(x)) w.r.t json representation.
    // See test_json_serialization_by_representation for more detail.

    Json::Value reference_json; // aka x
    reference_json["exposure_type"]   = Json::Int(1);
    reference_json["time"]["seconds"] = Json::UInt(2);

    test_json_serialization_by_representation<mio::abm::Vaccination>(reference_json);
}

TEST(TestAbmSerialization, Infection)
{
    // Test that a json value x is equal to serialize(deserialize(x)) w.r.t json representation.
    // See test_json_serialization_by_representation for more detail.

    unsigned i = 1; // counter s.t. members have different values

    Json::Value viral_load;
    viral_load["decline"]               = Json::Value((double)i++);
    viral_load["end_date"]["seconds"]   = Json::UInt(i++);
    viral_load["incline"]               = Json::Value((double)i++);
    viral_load["peak"]                  = Json::Value((double)i++);
    viral_load["start_date"]["seconds"] = Json::UInt(i++);

    Json::Value reference_json; // aka x
    reference_json["infection_course"] = Json::Value(Json::arrayValue);
    reference_json["virus_variant"]    = Json::UInt(0);
    reference_json["viral_load"]       = viral_load;
    reference_json["log_norm_alpha"]   = Json::Value((double)i++);
    reference_json["log_norm_beta"]    = Json::Value((double)i++);
    reference_json["detected"]         = Json::Value((bool)0);

    test_json_serialization_by_representation<mio::abm::Infection>(reference_json);
}

TEST(TestAbmSerialization, TestingScheme)
{
    // Test that a json value x is equal to serialize(deserialize(x)) w.r.t json representation.
    // See test_json_serialization_by_representation for more detail.

    mio::abm::TestingScheme testing_scheme(mio::abm::TestingCriteria({}, {}), mio::abm::TimeSpan(1),
                                           mio::abm::TimePoint(2), mio::abm::TimePoint(3),
                                           mio::abm::TestParameters{{4.0}, {5.0}}, 6.0);

    Json::Value test_parameters;
    test_parameters["sensitivity"] = mio::serialize_json(mio::UncertainValue<double>{4.0}).value();
    test_parameters["specitivity"] = mio::serialize_json(mio::UncertainValue<double>{5.0}).value();

    Json::Value reference_json; // aka x
    reference_json["criteria"]                            = Json::Value(Json::arrayValue);
    reference_json["min_time_since_last_test"]["seconds"] = Json::UInt(1);
    reference_json["start_date"]["seconds"]               = Json::UInt(2);
    reference_json["end_date"]["seconds"]                 = Json::UInt(3);
    reference_json["test_params"]                         = test_parameters;
    reference_json["probability"]                         = Json::Value((double)6);
    reference_json["is_active"]                           = Json::Value((bool)0);

    test_json_serialization_full(testing_scheme, reference_json);
}

TEST(TestAbmSerialization, TestingStrategy)
{
    // Test that a json value x is equal to serialize(deserialize(x)) w.r.t json representation.
    // See test_json_serialization_by_representation for more detail.

    unsigned i = 1; // counter s.t. members have different values

    Json::Value local_strategy;
    local_strategy["id"]      = Json::UInt(i++);
    local_strategy["schemes"] = Json::Value(Json::arrayValue);
    local_strategy["type"]    = Json::UInt(i++);

    Json::Value reference_json; // aka x
    reference_json["schemes"][0] = local_strategy;

    test_json_serialization_by_representation<mio::abm::TestingScheme>(reference_json);
}

TEST(TestAbmSerialization, Person)
{
    // Test that a json value x is equal to serialize(deserialize(x)) w.r.t json representation.
    // See test_json_serialization_by_representation for more detail.

    auto json_uint_array = [](std::vector<uint32_t> values) {
        return mio::serialize_json(values).value();
    };
    auto json_double_array = [](std::vector<double> values) {
        return mio::serialize_json(values).value();
    };

    unsigned i = 1; // counter s.t. members have different values

    Json::Value reference_json; // aka x
    reference_json["age_group"]           = Json::UInt(i++);
    reference_json["assigned_locations"]  = json_uint_array({i++, i++, i++, i++, i++, i++, i++, i++, i++, i++, i++});
    reference_json["cells"]               = json_uint_array({i++});
    reference_json["id"]                  = Json::UInt(i++);
    reference_json["infections"]          = Json::Value(Json::arrayValue);
    reference_json["last_transport_mode"] = Json::UInt(i++);
    reference_json["location"]            = Json::UInt(i++);
    reference_json["mask"]["mask_type"]   = Json::UInt(0);
    reference_json["mask"]["time_used"]["seconds"] = Json::UInt(i++);
    reference_json["mask_compliance"] =
        json_double_array({(double)i++, (double)i++, (double)i++, (double)i++, (double)i++, (double)i++, (double)i++,
                           (double)i++, (double)i++, (double)i++, (double)i++});
    reference_json["quarantine_start"]["seconds"]  = Json::UInt(i++);
    reference_json["rnd_go_to_school_hour"]        = Json::Value((double)i++);
    reference_json["rnd_go_to_work_hour"]          = Json::Value((double)i++);
    reference_json["rnd_schoolgroup"]              = Json::Value((double)i++);
    reference_json["rnd_workgroup"]                = Json::Value((double)i++);
    reference_json["rng_counter"]                  = Json::UInt(i++);
    reference_json["time_at_location"]["seconds"]  = Json::UInt(i++);
    reference_json["time_of_last_test"]["seconds"] = Json::UInt(i++);
    reference_json["vaccinations"]                 = Json::Value(Json::arrayValue);
    reference_json["wears_mask"]                   = Json::Value(false);

    test_json_serialization_by_representation<mio::abm::Person>(reference_json);
}

TEST(TestAbmSerialization, Location)
{
    // Test that a json value x is equal to serialize(deserialize(x)) w.r.t json representation.
    // See test_json_serialization_by_representation for more detail.

    unsigned i = 1; // counter s.t. members have different values

    Json::Value reference_json; // aka x
    reference_json["cells"][0]["capacity"]["persons"]    = Json::UInt(i++);
    reference_json["cells"][0]["capacity"]["volume"]     = Json::UInt(i++);
    reference_json["geographical_location"]["latitude"]  = Json::Value((double)i++);
    reference_json["geographical_location"]["longitude"] = Json::Value((double)i++);
    reference_json["id"]                                 = Json::UInt(i++);
    reference_json["npi_active"]                         = Json::Value(false);
    reference_json["parameters"]["ContactRates"] =
        mio::serialize_json(mio::abm::ContactRates::get_default(i++)).value();
    reference_json["parameters"]["MaximumContacts"]                     = Json::Value((double)i++);
    reference_json["parameters"]["UseLocationCapacityForTransmissions"] = Json::Value(false);
    reference_json["required_mask"]                                     = Json::UInt(0);

    test_json_serialization_by_representation<mio::abm::Location>(reference_json);
}

TEST(TestAbmSerialization, World)
{
    // Test that a json value x is equal to serialize(deserialize(x)) w.r.t json representation.
    // See test_json_serialization_by_representation for more detail.

    auto json_uint_array = [](std::vector<uint32_t> values) {
        return mio::serialize_json(values).value();
    };

    unsigned i = 1; // counter s.t. members have different values

    Json::Value reference_json; // aka x
    reference_json["cemetery_id"]                 = Json::UInt(i++);
    reference_json["location_types"]              = Json::UInt(i++);
    reference_json["locations"]                   = Json::Value(Json::arrayValue);
    reference_json["parameters"]                  = mio::serialize_json(mio::abm::Parameters(i++)).value();
    reference_json["persons"]                     = Json::Value(Json::arrayValue);
    reference_json["rng"]["counter"]              = Json::UInt(i++);
    reference_json["rng"]["key"]                  = Json::UInt(i++);
    reference_json["rng"]["seeds"]                = json_uint_array({i++, i++, i++, i++, i++, i++});
    reference_json["testing_strategy"]["schemes"] = Json::Value(Json::arrayValue);
    reference_json["trip_list"]["index"]          = Json::UInt(i++);
    reference_json["trip_list"]["trips_weekday"]  = Json::Value(Json::arrayValue);
    reference_json["trip_list"]["trips_weekend"]  = Json::Value(Json::arrayValue);
    reference_json["use_migration_rules"]         = Json::Value(true);

    test_json_serialization_by_representation<mio::abm::World>(reference_json);
}

#endif
