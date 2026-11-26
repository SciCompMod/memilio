#ifndef LOCATION_SPLIT_PARAMETER_SETTER_H
#define LOCATION_SPLIT_PARAMETER_SETTER_H

#include "city_parameters.h"

#include "abm/model.h"
#include "abm/parameters.h"

#include <string>
#include <utility>

/**
 * @brief Convert the mean and standard deviation of a distribution to the log mean and log standard
 *        deviation of the corresponding log normal distribution.
 */
std::pair<double, double> get_mu_and_sigma(std::pair<double, double> mean_and_std);

/// @brief Set the global infection and mobility parameters of the ABM.
void set_parameters(mio::abm::Parameters& params);

/**
 * @brief Mean hours a Person spends at a location type on a day it goes there, measured from the
 * mobility rules of the ABM over a whole week (10,000 persons, region Germany). The contact matrices
 * count contacts per day, while the ABM applies the contact rates only while a Person is at the
 * Location, so the matrices are scaled with 24 h / these hours.
 */
struct HoursPerVisitDay {
    static constexpr double home         = 19.3;
    static constexpr double school       = 8.0;
    static constexpr double work         = 10.0;
    static constexpr double social_event = 4.0;
    static constexpr double shop         = 1.2;
};

/**
 * @brief Transmission weight per location type: relative infectiousness of one contact hour there.
 * The contact rate of a location type is matrix * (24 / HoursPerVisitDay) * weight.
 *
 * CALIBRATED DEFAULTS. Fitted with calibration/fit_contact_weights.py to the sources of infection of
 * adults in the French ComCor study (Home 42.5%, SocialEvent 36.7%, Work 20.7%). They belong together
 * with Config::DEFAULT_INFECTION_K = 2.11. School and shop have no target data and are not fitted.
 */
struct ContactWeights {
    double home         = 21.616;
    double school       = 1.6;
    double work         = 0.625;
    double social_event = 4.8;
    double shop         = 0.1584;
};

/**
 * @brief Parse "home=1.2,work=0.8,..." into @p weights. Keys not given keep their value.
 * @return False if a key is unknown or a value is not a non-negative number.
 */
bool parse_contact_weights(const std::string& text, ContactWeights& weights);

/// @brief The weights as "home=...,school=...,work=...,social_event=...,shop=...".
std::string contact_weights_to_string(const ContactWeights& weights);

/**
 * @brief Set the location specific contact rates for all Location%s of the Model.
 * @param[in,out] model The Model whose Location%s are parameterized.
 * @param[in] region Region whose contact matrices are used.
 * @param[in] weights Transmission weight per location type.
 */
void set_local_parameters(mio::abm::Model& model, CityParameters::Region region,
                          const ContactWeights& weights = ContactWeights{});

#endif // LOCATION_SPLIT_PARAMETER_SETTER_H
