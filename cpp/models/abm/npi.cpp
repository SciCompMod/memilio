/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: David Kerkman
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

#include "abm/npi.h"
#include "abm/person.h"
#include "abm/location.h"

namespace mio
{
namespace abm
{

/**
 * Base NPI, but equipped with a stringency which is often useful for all kind of NPIs.
 * Alternatively, this can also be implemented in each NPI directly.
*/
struct AcceptBaseStringency : AcceptBase {
    AcceptBaseStringency() = default;
    AcceptBaseStringency(double strin)
        : stringency(strin)
    {
    }
    double stringency = 1.;
};

/**
 * @brief NPI to require the usage of masks for entering the Location.
*/
struct MaskRequired : AcceptBaseStringency {

    MaskRequired(MaskType mt)
        : mask_type(mt)
    {
    }

    MaskRequired(MaskType mt, double strin)
        : AcceptBaseStringency(strin)
        , mask_type(mt)
    {
    }

    bool operator()(Person::RandomNumberGenerator& rng, const Person& p, const Location& loc, TimePoint /*t*/) override
    {
        return p.apply_mask_intervention(rng, loc);
    }

    MaskType mask_type;
};

/**
 * @brief NPI to set the maximum capacity for entering for the Location. Use 0 for a complete shutdown.
*/
struct CapacityAbsolute : AcceptBaseStringency {

    CapacityAbsolute(uint32_t cap)
        : capacity(cap)
    {
    }

    CapacityAbsolute(uint32_t cap, double strin)
        : AcceptBaseStringency(strin)
        , capacity(cap)
    {
    }

    bool operator()(Person::RandomNumberGenerator& /*rng*/, const Person& /*p*/, const Location& loc,
                    TimePoint /*t*/) override
    {
        return loc.get_number_persons() < capacity;
    }

    uint32_t capacity;
};

/**
 * @brief NPI to set the relative capacity change for entering the Location. Use 0 for a complete shutdown.
*/
struct CapacityRelative : AcceptBaseStringency {

    CapacityRelative(double cap_fac)
        : capacity_factor(cap_fac)
    {
    }

    CapacityRelative(double cap_fac, double strin)
        : AcceptBaseStringency(strin)
        , capacity_factor(cap_fac)
    {
    }

    bool operator()(Person::RandomNumberGenerator& /*rng*/, const Person& /*p*/, const Location& loc,
                    TimePoint /*t*/) override
    {
        return loc.get_number_persons() < capacity_factor * loc.get_capacity().persons;
    }

    double capacity_factor;
};

/**
 * @brief NPI to set the required tests for entering the Location.
*/
struct TestRequired : AcceptBaseStringency {

    TestRequired(TestingScheme& ts)
        : testing_scheme(&ts)
    {
    }

    TestRequired(TestingScheme& ts, double strin)
        : AcceptBaseStringency(strin)
        , testing_scheme(&ts)
    {
    }

    bool operator()(Person::RandomNumberGenerator& rng, const Person& p, const Location& loc, TimePoint t) override
    {
        return testing_scheme->run_scheme(rng, p, loc, t);
    }

    TestingScheme* testing_scheme;
};

struct NPISet : AcceptBase {
    LocationNPI mask{AcceptBase{}};
    LocationNPI capacity{AcceptBase{}};
    LocationNPI test{AcceptBase{}};

    bool operator()(Person::RandomNumberGenerator& rng, const Person& p, const Location& loc, TimePoint t) override
    {
        return mask(rng, p, loc, t) && capacity(rng, p, loc, t) && test(rng, p, loc, t);
    }
};

} // namespace abm
} // namespace mio
