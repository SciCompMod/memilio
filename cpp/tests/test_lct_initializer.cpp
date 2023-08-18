/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke
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

#include "lct_secir/infection_state.h"
#include "lct_secir/initialization.h"
#include "memilio/config.h"
#include "memilio/epidemiology/state_age_function.h"

#include <gtest/gtest.h>

TEST(TestInitializer, testErlang)
{
    ScalarType rate[]  = {0.5, 0.8, 1, 3};
    int shape[]        = {1, 3, 5, 10, 20};
    ScalarType times[] = {0, 0.5, 3, 5.555, 20, 70.34};
    for (int r = 0; r < 4; r++) {
        for (int s = 0; s < 5; s++) {
            mio::ErlangSurvivalFunction survival(rate[r], shape[s]);
            EXPECT_EQ(survival.eval(0), 1.0);
            mio::lsecir::ErlangDensity density(rate[r], 1);

            for (int tau = 0; tau < 6; tau++) {
                ScalarType f = 0;
                for (int k = 1; k < shape[s] + 1; k++) {
                    density.set_parameter(k);
                    f += density.eval(times[tau]);
                }
                EXPECT_NEAR(f * (1 / rate[r]), survival.eval(times[tau]), 1e-10);
            }
        }
    }
}