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

#include "memilio/io/epi_data.h"
#include "matchers.h"
#include "gtest/gtest.h"

TEST(TestRkiDataIo, read) {
    Json::Value js(Json::arrayValue);
    js[0]["ID_County"] = 1001;
    js[0]["Date"] = "2021-12-01";
    js[0]["Confirmed"] = 1;
    js[0]["Deaths"] = 2;
    js[0]["Recovered"] = 3;
    js[0]["Age_RKI"] = "A80+";

    js[1]["ID_County"] = 1001;
    js[1]["Date"] = "2021-12-02";
    js[1]["Confirmed"] = 3;
    js[1]["Deaths"] = 4;
    js[1]["Recovered"] = 5;
    js[1]["Age_RKI"] = "A00-A04";

    js[2]["ID_County"] = 1002;
    js[2]["Date"] = "2021-12-02";
    js[2]["Confirmed"] = 3;
    js[2]["Deaths"] = 4;
    js[2]["Recovered"] = 5;
    js[2]["Age_RKI"] = "unknown";

    auto result = mio::deserialize_rki_data(js);
    ASSERT_THAT(print_wrap(result), IsSuccess());
    auto rki_data = result.value();
    ASSERT_EQ(rki_data.size(), 2);

    ASSERT_EQ(rki_data[0].age_group, mio::AgeGroup(5));
    ASSERT_EQ(rki_data[0].date, mio::Date(2021, 12, 1));
    ASSERT_EQ(rki_data[0].num_confirmed, 1);
    ASSERT_EQ(rki_data[0].num_deceased, 2);
    ASSERT_EQ(rki_data[0].num_recovered, 3);
    ASSERT_EQ(rki_data[0].county_id, mio::regions::de::CountyId(1001));
    ASSERT_EQ(rki_data[0].state_id, boost::none);

    ASSERT_EQ(rki_data[1].age_group, mio::AgeGroup(0));
    ASSERT_EQ(rki_data[1].date, mio::Date(2021, 12, 2));
    ASSERT_EQ(rki_data[1].num_confirmed, 3);
    ASSERT_EQ(rki_data[1].num_deceased, 4);
    ASSERT_EQ(rki_data[1].num_recovered, 5);
    ASSERT_EQ(rki_data[1].county_id, mio::regions::de::CountyId(1001));
    ASSERT_EQ(rki_data[1].state_id, boost::none);
}