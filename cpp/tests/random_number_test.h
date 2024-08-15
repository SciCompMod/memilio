/* 
* Copyright (C) 2020-2024 MEmilio
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
#include "memilio/utils/random_number_generator.h"

#include <gtest/gtest.h>

class RandomNumberTest : public ::testing::Test
{
public:
    /**
     * @brief Draws a uniformly distributed random number.
     * @tparam FP A floating point type, defaults to double.
     * @param[in] min, max Lower and upper bound to the uniform distribution.
     * @return A random value between min and max. 
     */
    template <class FP = double>
    double random_number(FP min = FP{-1e+3}, FP max = {1e+3})
    {
        return mio::UniformDistribution<FP>::get_instance()(m_rng, min, max);
    }

    /// @brief Access the random number generator. Should only be used for debugging.
    mio::RandomNumberGenerator get_rng()
    {
        return m_rng;
    }

protected:
    void SetUp() override
    {
        log_rng_seeds(m_rng, mio::LogLevel::warn);
    }

private:
    mio::RandomNumberGenerator m_rng{}; ///< Seeded rng used by this test fixture.
};
