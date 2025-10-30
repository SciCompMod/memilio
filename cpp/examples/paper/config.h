/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julia Bicker, Anna Wendler
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

// Define parameters for the simulation.
#include "abm/location.h"
#include "abm/model.h"
#include "abm/simulation.h"
#include "abm/time.h"
#include "abm/virus_variant.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/io/history.h"
#include "memilio/utils/compiler_diagnostics.h"
#include <vector>
namespace params
{
const size_t num_age_groups               = 6;
const size_t scaling_factor               = 1000;
const std::vector<double> age_group_sizes = {3969138.0, 7508662, 18921292, 28666166, 18153339, 5936434};
const ScalarType total_population         = 83155031. / scaling_factor;
// const ScalarType total_confirmed_cases = 341223.;
// const ScalarType deaths = 0.;

// Define transition probabilities per age group.
const ScalarType infectedSymptomsPerInfectedNoSymptoms[] = {0.75, 0.75, 0.8, 0.8, 0.8, 0.8};
const ScalarType severePerInfectedSymptoms[]             = {0.0075, 0.0075, 0.019, 0.0615, 0.165, 0.225};
const ScalarType criticalPerSevere[]                     = {0.075, 0.075, 0.075, 0.15, 0.3, 0.4};
const ScalarType deathsPerCritical[]                     = {0.05, 0.05, 0.14, 0.14, 0.4, 0.6};

// Define lognormal parameters. For each transition, we need shape and scale.
// These are given in this order below. The transition distributions are the same for all age groups.
const ScalarType lognorm_EtINS[]     = {0.32459285, 4.26907484};
const ScalarType lognorm_INStISy[]   = {0.7158751, 0.85135303};
const ScalarType lognorm_INStR[]     = {0.24622068, 7.76114};
const ScalarType lognorm_ISytISev[]  = {0.66258947, 5.29920733};
const ScalarType lognorm_ISytR[]     = {0.24622068, 7.76114};
const ScalarType lognorm_ISevtICri[] = {1.01076765, 0.9};
const ScalarType lognorm_ISevtR[]    = {0.33816427, 17.09411753};
const ScalarType lognorm_ICritD[]    = {0.42819924, 9.76267505};
const ScalarType lognorm_ICritR[]    = {0.33816427, 17.09411753};

// Define mean stay times per age group. Note that these are different per age group as the transition probabilities
// differ between age groups.
const ScalarType timeExposed[]            = {4.5, 4.5, 4.5, 4.5, 4.5, 4.5};
const ScalarType timeInfectedNoSymptoms[] = {2.825, 2.825, 2.48, 2.48, 2.48, 2.48};
const ScalarType timeInfectedSymptoms[]   = {7.9895, 7.9895, 7.9734, 7.9139, 7.9139, 7.685};
const ScalarType timeInfectedSevere[]     = {16.855, 16.855, 16.855, 15.61, 13.12, 11.46};
const ScalarType timeInfectedCritical[]   = {17.73, 17.73, 17.064, 17.064, 15.14, 13.66};

// Define number of subcompartments for each compartment per age group. Note that these are different per age group as
// the transition probabilities differ between age groups.
// These values are used as a template argument and thus have to be constexpr.
constexpr size_t n_subcomps_E[]    = {9, 9, 9, 9, 9, 9};
constexpr size_t n_subcomps_INS[]  = {5, 5, 4, 4, 4, 4};
constexpr size_t n_subcomps_ISy[]  = {16, 16, 16, 15, 15, 13};
constexpr size_t n_subcomps_ISev[] = {8, 8, 8, 7, 6, 5};
constexpr size_t n_subcomps_ICri[] = {8, 8, 8, 8, 7, 6};

// Define epidemiological parameters.
const ScalarType transmissionProbabilityOnContact[] = {0.03, 0.06, 0.06, 0.06, 0.09, 0.175};
const ScalarType relativeTransmissionNoSymptoms     = 1;
const ScalarType riskOfInfectionFromSymptomatic     = 1.;
const ScalarType seasonality                        = 0.;
const ScalarType scale_confirmed_cases              = 1.;

// Define simulation parameters.
ScalarType t0   = 0.;
ScalarType tmax = 20;
ScalarType dt   = 0.01;

} // namespace params

namespace ABMparams
{

// Distribution of household sizes from 1 to 5
const ScalarType household_size_distribution[] = {0.415, 0.342, 0.118, 0.091, 0.034};
// Location parameters. All location sizes are sampled from a normal distribution which is capped at the given minimum value.
// Mean size and stddev of workplaces
const size_t workplace_size[] = {15, 25};
// Minimum workplace size
const size_t min_workplace_size = 5;
// Mean size and stddev of schools
const size_t school_size[] = {42, 3};
// Minimum school size
const size_t min_school_size = 15;
// Mean size and stddev of events
const size_t event_size[] = {42, 941};
// Minimum event size
const size_t min_event_size = 10;
// Mean size and stddev of shops
const size_t shop_size[] = {90, 1000};
// Minimum shop size
const size_t min_shop_size = 30;

} // namespace ABMparams

namespace ABMLoggers
{
struct LogTimePoint : mio::LogAlways {
    using Type = mio::abm::TimePoint;
    static Type log(const mio::abm::Simulation<mio::abm::Model>& sim)
    {
        return sim.get_time();
    }
};

struct LogExposureRate : mio::LogAlways {
    using Type = std::vector<double>;
    static Type log(const mio::abm::Simulation<mio::abm::Model>& sim)
    {
        std::vector<double> cum_rate(params::num_age_groups);
        auto& rates = sim.get_model().get_contact_exposure_rates();
        if (rates.size() > 0) {
            size_t num_locs = 0;
            for (auto& loc : sim.get_model().get_locations()) {
                num_locs += 1;
                for (size_t age = 0; age < params::num_age_groups; ++age) {
                    cum_rate[age] += rates[loc.get_id().get()][{mio::abm::CellIndex(0),
                                                                mio::abm::VirusVariant::Wildtype, mio::AgeGroup(age)}];
                }
            }
            // Divide rate by number of locations
            std::transform(cum_rate.begin(), cum_rate.end(), cum_rate.begin(), [num_locs](double x) {
                return x / num_locs;
            });
        }
        return cum_rate;
    }
};
} // namespace ABMLoggers
