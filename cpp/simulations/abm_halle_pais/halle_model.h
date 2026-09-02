/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Sascha Korf
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
#ifndef MIO_SIMULATIONS_ABM_HALLE_PAIS_HALLE_MODEL_H
#define MIO_SIMULATIONS_ABM_HALLE_PAIS_HALLE_MODEL_H

#include "abm/model.h"
#include "abm/time.h"
#include "memilio/io/io.h"
#include "memilio/utils/date.h"
#include "memilio/utils/random_number_generator.h"

#include <string>
#include <vector>

namespace mio
{
namespace halle
{

/// @brief Number of age groups used throughout this simulation.
constexpr size_t num_age_groups = 6;

/// @brief County id of Halle (Saale), Stadt. The surrounding Saalekreis is 15088.
constexpr int halle_county_id = 15002;

/**
 * @brief One parameter estimated by the fit, together with the bounds of its uniform prior.
 *
 * The bounds of the acute-infection parameters are taken from the grid search of the ABM paper
 * (cpp/simulations/paper_abm_bs_testing.cpp on branch abm_paper_test_bs).
 */
struct FitParameter {
    std::string name; ///< Identifier, also used as column name in the design matrix written for the fit.
    double lower; ///< Lower bound of the uniform prior.
    double upper; ///< Upper bound of the uniform prior.
};

/**
 * @brief The parameters estimated by the fit, in the order used by every parameter vector.
 *
 * This table is the single place where the fitted parameters are defined: adding or removing an entry
 * changes the dimension of the fit without any other code change. No intervention is modelled at present,
 * so neither a contact reduction nor a testing frequency appears here.
 */
const std::vector<FitParameter>& fit_parameters();

/**
 * @brief The channels of the summary statistic that the fit compares against real data.
 * @return Channel names, in the order used by every observable vector.
 */
const std::vector<std::string>& observable_channels();

/**
 * @brief Draw a parameter vector from the (uniform) prior.
 * @param[in,out] rng Random number generator.
 * @return A parameter vector with one entry per entry of fit_parameters().
 */
std::vector<double> sample_prior(RandomNumberGenerator& rng);

/**
 * @brief Check whether a parameter vector lies inside the prior support.
 * @param[in] theta A parameter vector with one entry per entry of fit_parameters().
 */
bool is_in_prior_support(const std::vector<double>& theta);

/**
 * @brief Everything needed to build the Halle model that is not being fitted.
 *
 * Every file path is required unless stated otherwise. A path that cannot be read is a hard error and
 * never a silent fallback: a run against missing data would otherwise look successful while producing
 * all-zero output.
 */
struct ModelSetup {
    std::string person_file{}; ///< CSV with columns age,home_id,school_id,work_id,shopping_id,event_id.
    std::string cases_file{}; ///< CSV with columns date,age_group,new_cases (reported), for the infection history.
    std::string vaccinations_file{}; ///< CSV with columns date,age_group,new_doses, for the vaccination history.
    std::string contact_dir{}; ///< Directory with the German baseline contact matrices.
    /**
     * @brief Length of the window before t0 from which the histories are seeded.
     *
     * 90 days by default. Longer windows overcount immunity, because a Person is seeded at most once here
     * while the reported cases they are drawn from include reinfections: over a full year around the
     * Omicron wave, reported cases times a dark figure of 4 exceed the population of Halle, which seeds
     * every agent as recovered and leaves no susceptible for the epidemic to run in. 90 days also matches
     * the window over which SeverityProtectionFactor is defined.
     */
    int history_lookback_days = 90;
    /**
     * @brief Allow building a model without infection and vaccination history.
     *
     * Only intended for smoke tests. A fit must not use this: without the histories, the population has no
     * immunity and no pre-existing PAIS at t0, which silently changes the epidemic that is being fitted.
     */
    bool allow_missing_history = false;
};

/**
 * @brief Build the Halle ABM.
 *
 * Builds the population and its locations from the Halle population file, then seeds the infection and
 * vaccination histories so that immunity and pre-existing PAIS are in their correct state at @p t0.
 *
 * @param[in] setup The non-fitted parts of the model.
 * @param[in] theta The fitted parameters, see fit_parameters().
 * @param[in] start_date Calendar date that @p t0 corresponds to, used to align the input data.
 * @param[in] t0 Start of the simulation.
 * @param[in] tmax End of the simulation. Currently unused; kept because a testing scheme would need it as
 *                 its validity period once testing is reinstated.
 * @param[in] rng Random number generator used for the whole model, so a run is reproducible from its seed.
 * @return The model, or an error if an input file could not be read.
 */
IOResult<abm::Model> make_model(const ModelSetup& setup, const std::vector<double>& theta, Date start_date,
                                abm::TimePoint t0, abm::TimePoint tmax, const RandomNumberGenerator& rng);

} // namespace halle
} // namespace mio

#endif
