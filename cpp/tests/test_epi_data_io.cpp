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

#include "abm/infection_state.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/geography/regions.h"
#include "memilio/io/epi_data.h"
#include "matchers.h"
#include "memilio/io/io.h"
#include "memilio/io/mobility_io.h"
#include "test_data_dir.h"
#include "gtest/gtest.h"
#include "json/value.h"
#include "ode_secirvvs/model.h"
#include "ode_secirvvs/parameters_io.h"
#include "memilio/utils/stl_util.h"
#include "boost/optional/optional_io.hpp"
#include <gmock/gmock-matchers.h>

TEST(TestEpiDataIo, read_rki)
{
    Json::Value js(Json::arrayValue);
    js[0]["ID_County"] = 1001;
    js[0]["Date"]      = "2021-12-01";
    js[0]["Confirmed"] = 1;
    js[0]["Deaths"]    = 2;
    js[0]["Recovered"] = 3;
    js[0]["Age_RKI"]   = "A80+";

    js[1]["ID_County"] = 1001;
    js[1]["Date"]      = "2021-12-02";
    js[1]["Confirmed"] = 3;
    js[1]["Deaths"]    = 4;
    js[1]["Recovered"] = 5;
    js[1]["Age_RKI"]   = "A00-A04";

    js[2]["ID_County"] = 1002;
    js[2]["Date"]      = "2021-12-02";
    js[2]["Confirmed"] = 3;
    js[2]["Deaths"]    = 4;
    js[2]["Recovered"] = 5;
    js[2]["Age_RKI"]   = "unknown";

    auto result = mio::deserialize_confirmed_cases_data(js);
    ASSERT_THAT(print_wrap(result), IsSuccess());
    auto rki_data = result.value();
    ASSERT_EQ(rki_data.size(), 2);

    ASSERT_EQ(rki_data[0].age_group, mio::AgeGroup(5));
    ASSERT_EQ(rki_data[0].date, mio::Date(2021, 12, 1));
    ASSERT_EQ(rki_data[0].num_confirmed, 1);
    ASSERT_EQ(rki_data[0].num_deaths, 2);
    ASSERT_EQ(rki_data[0].num_recovered, 3);
    ASSERT_EQ(rki_data[0].county_id, mio::regions::CountyId(1001));
    ASSERT_EQ(rki_data[0].state_id, boost::none);

    ASSERT_EQ(rki_data[1].age_group, mio::AgeGroup(0));
    ASSERT_EQ(rki_data[1].date, mio::Date(2021, 12, 2));
    ASSERT_EQ(rki_data[1].num_confirmed, 3);
    ASSERT_EQ(rki_data[1].num_deaths, 4);
    ASSERT_EQ(rki_data[1].num_recovered, 5);
    ASSERT_EQ(rki_data[1].county_id, mio::regions::CountyId(1001));
    ASSERT_EQ(rki_data[1].state_id, boost::none);
}

TEST(TestEpiDataIo, read_rki_error_age)
{
    Json::Value js(Json::arrayValue);
    js[0]["ID_County"] = 1001;
    js[0]["Date"]      = "2021-12-01";
    js[0]["Confirmed"] = 1;
    js[0]["Deaths"]    = 2;
    js[0]["Recovered"] = 3;
    js[0]["Age_RKI"]   = "A01-A05"; //error

    auto result = mio::deserialize_confirmed_cases_data(js);
    ASSERT_THAT(print_wrap(result), IsFailure(mio::StatusCode::InvalidValue));
}

TEST(TestEpiDataIo, read_divi)
{
    Json::Value js(Json::arrayValue);
    js[0]["ID_County"] = 1001;
    js[0]["ICU"]       = 10.0;
    js[0]["Date"]      = "2022-10-05";

    js[1]["ID_County"] = 1002;
    js[1]["ICU"]       = 20.0;
    js[1]["Date"]      = "2022-10-07";

    auto r = mio::deserialize_divi_data(js);
    ASSERT_THAT(print_wrap(r), IsSuccess());

    auto& divi_data = r.value();
    ASSERT_EQ(divi_data.size(), 2);

    ASSERT_EQ(divi_data[0].county_id, mio::regions::CountyId(1001));
    ASSERT_EQ(divi_data[0].date, mio::Date(2022, 10, 5));
    ASSERT_EQ(divi_data[0].num_icu, 10.0);

    ASSERT_EQ(divi_data[1].county_id, mio::regions::CountyId(1002));
    ASSERT_EQ(divi_data[1].date, mio::Date(2022, 10, 7));
    ASSERT_EQ(divi_data[1].num_icu, 20.0);
}

TEST(TestEpiDataIo, read_population)
{
    Json::Value js(Json::arrayValue);
    js[0]["ID_County"]                                   = 1001;
    js[0][mio::PopulationDataEntry::age_group_names[0]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[1]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[2]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[3]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[4]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[5]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[6]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[7]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[8]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[9]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[10]] = 10;

    js[1]["ID_County"]                                   = 1002;
    js[1][mio::PopulationDataEntry::age_group_names[0]]  = 10;
    js[1][mio::PopulationDataEntry::age_group_names[1]]  = 20;
    js[1][mio::PopulationDataEntry::age_group_names[2]]  = 30;
    js[1][mio::PopulationDataEntry::age_group_names[3]]  = 40;
    js[1][mio::PopulationDataEntry::age_group_names[4]]  = 50;
    js[1][mio::PopulationDataEntry::age_group_names[5]]  = 60;
    js[1][mio::PopulationDataEntry::age_group_names[6]]  = 70;
    js[1][mio::PopulationDataEntry::age_group_names[7]]  = 80;
    js[1][mio::PopulationDataEntry::age_group_names[8]]  = 90;
    js[1][mio::PopulationDataEntry::age_group_names[9]]  = 100;
    js[1][mio::PopulationDataEntry::age_group_names[10]] = 110;

    auto r = mio::deserialize_population_data(js);
    ASSERT_THAT(print_wrap(r), IsSuccess());

    auto& population_data = r.value();
    ASSERT_EQ(population_data.size(), 2);

    ASSERT_EQ(population_data[0].county_id, mio::regions::CountyId(1001));
    ASSERT_THAT(population_data[0].population,
                testing::ElementsAre(testing::DoubleEq(10.0 + 2 * 10.0 / 3), testing::DoubleEq(10.0 / 3 + 10.0),
                                     testing::DoubleEq(10.0 + 10.0 + 10.0 + 0.5 * 10.0),
                                     testing::DoubleEq(0.5 * 10.0 + 10.0 + 2 * 10.0 / 3),
                                     testing::DoubleEq(10.0 / 3 + 10.0 + 0.2 * 10.0), testing::DoubleEq(0.8 * 10.0)));

    ASSERT_EQ(population_data[1].county_id, mio::regions::CountyId(1002));
    ASSERT_THAT(population_data[1].population,
                testing::ElementsAre(testing::DoubleEq(10.0 + 2 * 20.0 / 3), testing::DoubleEq(20.0 / 3 + 30.0),
                                     testing::DoubleEq(40.0 + 50.0 + 60.0 + 0.5 * 70.0),
                                     testing::DoubleEq(0.5 * 70.0 + 80.0 + 2 * 90.0 / 3),
                                     testing::DoubleEq(90.0 / 3 + 100.0 + 0.2 * 110.0),
                                     testing::DoubleEq(0.8 * 110.0)));
}

TEST(TestEpiDataIo, read_population_error_age)
{
    Json::Value js(Json::arrayValue);
    js[0]["ID_County"]                                   = 1001;
    js[0]["< 4 years"]                                   = 10; //error
    js[0][mio::PopulationDataEntry::age_group_names[1]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[2]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[3]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[4]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[5]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[6]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[7]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[8]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[9]]  = 10;
    js[0][mio::PopulationDataEntry::age_group_names[10]] = 10;

    auto r = mio::deserialize_population_data(js);
    ASSERT_THAT(print_wrap(r), IsFailure(mio::StatusCode::KeyNotFound));
}

TEST(TestEpiDataIo, read_county_ids)
{
    std::vector<int> true_ids = {
        1001,  1002,  1003,  1004,  1051,  1053,  1054,  1055,  1056,  1057,  1058,  1059,  1060,  1061,  1062,  2000,
        3101,  3102,  3103,  3151,  3153,  3154,  3155,  3157,  3158,  3159,  3241,  3251,  3252,  3254,  3255,  3256,
        3257,  3351,  3352,  3353,  3354,  3355,  3356,  3357,  3358,  3359,  3360,  3361,  3401,  3402,  3403,  3404,
        3405,  3451,  3452,  3453,  3454,  3455,  3456,  3457,  3458,  3459,  3460,  3461,  3462,  4011,  4012,  5111,
        5112,  5113,  5114,  5116,  5117,  5119,  5120,  5122,  5124,  5154,  5158,  5162,  5166,  5170,  5314,  5315,
        5316,  5334,  5358,  5362,  5366,  5370,  5374,  5378,  5382,  5512,  5513,  5515,  5554,  5558,  5562,  5566,
        5570,  5711,  5754,  5758,  5762,  5766,  5770,  5774,  5911,  5913,  5914,  5915,  5916,  5954,  5958,  5962,
        5966,  5970,  5974,  5978,  6411,  6412,  6413,  6414,  6431,  6432,  6433,  6434,  6435,  6436,  6437,  6438,
        6439,  6440,  6531,  6532,  6533,  6534,  6535,  6611,  6631,  6632,  6633,  6634,  6635,  6636,  7111,  7131,
        7132,  7133,  7134,  7135,  7137,  7138,  7140,  7141,  7143,  7211,  7231,  7232,  7233,  7235,  7311,  7312,
        7313,  7314,  7315,  7316,  7317,  7318,  7319,  7320,  7331,  7332,  7333,  7334,  7335,  7336,  7337,  7338,
        7339,  7340,  8111,  8115,  8116,  8117,  8118,  8119,  8121,  8125,  8126,  8127,  8128,  8135,  8136,  8211,
        8212,  8215,  8216,  8221,  8222,  8225,  8226,  8231,  8235,  8236,  8237,  8311,  8315,  8316,  8317,  8325,
        8326,  8327,  8335,  8336,  8337,  8415,  8416,  8417,  8421,  8425,  8426,  8435,  8436,  8437,  9161,  9162,
        9163,  9171,  9172,  9173,  9174,  9175,  9176,  9177,  9178,  9179,  9180,  9181,  9182,  9183,  9184,  9185,
        9186,  9187,  9188,  9189,  9190,  9261,  9262,  9263,  9271,  9272,  9273,  9274,  9275,  9276,  9277,  9278,
        9279,  9361,  9362,  9363,  9371,  9372,  9373,  9374,  9375,  9376,  9377,  9461,  9462,  9463,  9464,  9471,
        9472,  9473,  9474,  9475,  9476,  9477,  9478,  9479,  9561,  9562,  9563,  9564,  9565,  9571,  9572,  9573,
        9574,  9575,  9576,  9577,  9661,  9662,  9663,  9671,  9672,  9673,  9674,  9675,  9676,  9677,  9678,  9679,
        9761,  9762,  9763,  9764,  9771,  9772,  9773,  9774,  9775,  9776,  9777,  9778,  9779,  9780,  10041, 10042,
        10043, 10044, 10045, 10046, 11000, 12051, 12052, 12053, 12054, 12060, 12061, 12062, 12063, 12064, 12065, 12066,
        12067, 12068, 12069, 12070, 12071, 12072, 12073, 13003, 13004, 13071, 13072, 13073, 13074, 13075, 13076, 14511,
        14521, 14522, 14523, 14524, 14612, 14625, 14626, 14627, 14628, 14713, 14729, 14730, 15001, 15002, 15003, 15081,
        15082, 15083, 15084, 15085, 15086, 15087, 15088, 15089, 15090, 15091, 16051, 16052, 16053, 16054, 16055, 16056,
        16061, 16062, 16063, 16064, 16065, 16066, 16067, 16068, 16069, 16070, 16071, 16072, 16073, 16074, 16075, 16076,
        16077};

    std::string path = mio::path_join(TEST_DATA_DIR, "county_current_population.json");
    auto read_ids    = mio::get_node_ids(path, true);
    ASSERT_THAT(print_wrap(read_ids), IsSuccess());

    EXPECT_THAT(read_ids.value(), testing::ElementsAreArray(true_ids));
}

TEST(TestEpiDataIo, get_node_ids)
{
    std::vector<int> true_ids_district = {1234, 1235};

    std::vector<int> true_ids_county = {1001};

    std::string path       = mio::path_join(TEST_DATA_DIR, "test_current_population.json");
    auto read_ids_district = mio::get_node_ids(path, false);
    auto read_ids_county   = mio::get_node_ids(path, true);
    ASSERT_THAT(print_wrap(read_ids_district), IsSuccess());
    ASSERT_THAT(print_wrap(read_ids_county), IsSuccess());

    EXPECT_THAT(read_ids_district.value(), testing::ElementsAreArray(true_ids_district));
    EXPECT_THAT(read_ids_county.value(), testing::ElementsAreArray(true_ids_county));
}

TEST(TestEpiDataIo, read_divi_data)
{
    auto divi_data = mio::read_divi_data(mio::path_join(TEST_DATA_DIR, "test_county_divi.json")).value();

    ASSERT_EQ(divi_data.size(), 4);

    ASSERT_EQ(divi_data[0].district_id, mio::regions::DistrictId(1234));
    ASSERT_EQ(divi_data[0].date, mio::Date(2022, 04, 24));
    ASSERT_EQ(divi_data[0].num_icu, 0.5437);

    ASSERT_EQ(divi_data[1].district_id, mio::regions::DistrictId(1234));
    ASSERT_EQ(divi_data[1].date, mio::Date(2022, 04, 25));
    ASSERT_EQ(divi_data[1].num_icu, 0.64532);

    ASSERT_EQ(divi_data[2].district_id, mio::regions::DistrictId(1235));
    ASSERT_EQ(divi_data[2].date, mio::Date(2022, 04, 24));
    ASSERT_EQ(divi_data[2].num_icu, 0.3574);

    ASSERT_EQ(divi_data[3].district_id, mio::regions::DistrictId(1235));
    ASSERT_EQ(divi_data[3].date, mio::Date(2022, 04, 25));
    ASSERT_EQ(divi_data[3].num_icu, 0.35437);
}

TEST(TestEpiDataIo, read_confirmed_cases_data)
{
    auto case_data = mio::read_confirmed_cases_data(mio::path_join(TEST_DATA_DIR, "test_cases_all_age.json")).value();

    ASSERT_EQ(case_data.size(), 3);

    ASSERT_EQ(case_data[0].age_group, mio::AgeGroup(0));
    ASSERT_EQ(case_data[0].date, mio::Date(2022, 04, 24));
    ASSERT_EQ(case_data[0].num_confirmed, 1);
    ASSERT_EQ(case_data[0].num_deaths, 0);
    ASSERT_EQ(case_data[0].num_recovered, 0);
    ASSERT_EQ(case_data[0].district_id, mio::regions::DistrictId(1234));
    ASSERT_EQ(case_data[0].county_id, boost::none);
    ASSERT_EQ(case_data[0].state_id, boost::none);

    ASSERT_EQ(case_data[1].age_group, mio::AgeGroup(2));
    ASSERT_EQ(case_data[1].date, mio::Date(2022, 04, 25));
    ASSERT_EQ(case_data[1].num_confirmed, 20);
    ASSERT_EQ(case_data[1].num_deaths, 1);
    ASSERT_EQ(case_data[1].num_recovered, 5);
    ASSERT_EQ(case_data[1].district_id, mio::regions::DistrictId(1234));
    ASSERT_EQ(case_data[1].county_id, boost::none);
    ASSERT_EQ(case_data[1].state_id, boost::none);

    ASSERT_EQ(case_data[2].age_group, mio::AgeGroup(5));
    ASSERT_EQ(case_data[2].date, mio::Date(2022, 04, 24));
    ASSERT_EQ(case_data[2].num_confirmed, 15);
    ASSERT_EQ(case_data[2].num_deaths, 3);
    ASSERT_EQ(case_data[2].num_recovered, 2);
    ASSERT_EQ(case_data[2].district_id, mio::regions::DistrictId(1235));
    ASSERT_EQ(case_data[2].county_id, boost::none);
    ASSERT_EQ(case_data[2].state_id, boost::none);
}

TEST(TestEpiDataIO, read_vaccination_data)
{
    auto vacc_data = mio::read_vaccination_data(mio::path_join(TEST_DATA_DIR, "test_all_ageinf_vacc.json")).value();

    ASSERT_EQ(vacc_data.size(), 2);

    ASSERT_EQ(vacc_data[0].date, mio::Date(2022, 4, 12));
    ASSERT_EQ(vacc_data[0].age_group, mio::AgeGroup(0));
    ASSERT_EQ(vacc_data[0].county_id, mio::regions::CountyId(1001));
    ASSERT_EQ(vacc_data[0].district_id, mio::regions::DistrictId(1234));
    ASSERT_EQ(vacc_data[0].num_vaccinations_completed, 5.0);

    ASSERT_EQ(vacc_data[1].date, mio::Date(2022, 4, 15));
    ASSERT_EQ(vacc_data[1].age_group, mio::AgeGroup(2));
    ASSERT_EQ(vacc_data[1].county_id, boost::none);
    ASSERT_EQ(vacc_data[1].district_id, mio::regions::DistrictId(1235));
    ASSERT_EQ(vacc_data[1].num_vaccinations_completed, 1.0);
}

TEST(TestEpiData, set_vaccination_data)
{
    auto num_age_groups = 1;
    auto num_days       = 10;

    std::vector<int> county_ids = {1001};
    mio::osecirvvs::Model model(num_age_groups);
    model.parameters.set<mio::osecirvvs::VaccinationGap>(3);
    model.parameters.set<mio::osecirvvs::DaysUntilEffectivePartialImmunity>(1);
    model.parameters.set<mio::osecirvvs::DaysUntilEffectiveImprovedImmunity>(2);
    std::vector<mio::osecirvvs::Model> model_vector{model};

    auto f = mio::osecirvvs::details::set_vaccination_data(model_vector,
                                                           mio::path_join(TEST_DATA_DIR, "vaccination_test.json"),
                                                           mio::Date(2022, 4, 15), county_ids, num_days);

    auto expected_values_PV =
        (Eigen::ArrayXd(num_age_groups * (num_days + 1)) << 7, 8, 9, 9, 10, 12, 14, 16, 18, 20, 22).finished();

    auto expected_values_FV =
        (Eigen::ArrayXd(num_age_groups * (num_days + 1)) << 2, 4, 5, 5, 7, 8, 9, 9, 10, 12, 14).finished();

    ASSERT_THAT(print_wrap(model_vector[0].parameters.template get<mio::osecirvvs::DailyFullVaccination>().array()),
                MatrixNear(print_wrap(expected_values_FV), 1e-8, 1e-8));
    ASSERT_THAT(print_wrap(model_vector[0].parameters.template get<mio::osecirvvs::DailyFirstVaccination>().array()),
                MatrixNear(print_wrap(expected_values_PV), 1e-8, 1e-8));
}

TEST(TestEpiData, vaccination_data)
{
    auto js                 = Json::Value(Json::arrayValue);
    js[0]["Date"]           = "2021-12-01";
    js[0]["ID_County"]      = 1011;
    js[0]["Vacc_completed"] = 23.05;
    js[0]["Age_RKI"]        = "5-14";

    js[1]["Date"]           = "2021-12-02";
    js[1]["ID_County"]      = 1012;
    js[1]["Vacc_completed"] = 12.0;
    js[1]["Age_RKI"]        = "80-99";

    auto r = mio::deserialize_vaccination_data(js);
    ASSERT_THAT(print_wrap(r), IsSuccess());

    auto&& vacc_data = r.value();
    ASSERT_EQ(vacc_data.size(), 2);

    ASSERT_EQ(vacc_data[0].date, mio::Date(2021, 12, 1));
    ASSERT_EQ(vacc_data[0].age_group, mio::AgeGroup(1));
    ASSERT_EQ(vacc_data[0].county_id, mio::regions::CountyId(1011));
    ASSERT_EQ(vacc_data[0].num_vaccinations_completed, 23.05);

    ASSERT_EQ(vacc_data[1].date, mio::Date(2021, 12, 2));
    ASSERT_EQ(vacc_data[1].age_group, mio::AgeGroup(5));
    ASSERT_EQ(vacc_data[1].county_id, mio::regions::CountyId(1012));
    ASSERT_EQ(vacc_data[1].num_vaccinations_completed, 12.0);
}

TEST(TestEpiData, vaccination_data_error_age)
{
    auto js                 = Json::Value(Json::arrayValue);
    js[0]["Date"]           = "2021-12-01";
    js[0]["ID_County"]      = 1011;
    js[0]["Vacc_completed"] = 23.05;
    js[0]["Age_RKI"]        = "5-15"; //error

    js[1]["Date"]           = "2021-12-02";
    js[1]["ID_County"]      = 1012;
    js[1]["Vacc_completed"] = 12.0;
    js[1]["Age_RKI"]        = "80-99";

    auto r = mio::deserialize_vaccination_data(js);
    ASSERT_THAT(print_wrap(r), IsFailure(mio::StatusCode::InvalidValue));
}
