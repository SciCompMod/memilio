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
 * @brief Set the location specific contact rates for all Location%s of the Model.
 * @param[in,out] model The Model whose Location%s are parameterized.
 * @param[in] region Region whose contact matrices are used.
 */
/**
 * @brief Relative transmission weight per location type, multiplied onto the heuristic contact rate
 * factors in set_local_parameters.
 *
 * CALIBRATED DEFAULTS. The defaults below were fitted with calibration/fit_contact_weights.py to the
 * sources of infection of adults in the French ComCor study (Home 42.5%, SocialEvent 36.7%, Work
 * 20.7%), keeping the growth rate of the uncalibrated model (all weights 1, k = 2.4). They belong
 * together with Config::DEFAULT_INFECTION_K = 2.11. School and shop have no target data and keep the
 * weight of work. All weights 1 reproduce the original heuristic factors.
 */
struct ContactWeights {
    double home         = 1.28;
    double school       = 1.0;
    double work         = 1.0;
    double social_event = 2.0;
    double shop         = 1.0;
};

/**
 * @brief Parse "home=1.2,work=0.8,..." into @p weights. Keys not given keep their value.
 * @return False if a key is unknown or a value is not a non-negative number.
 */
bool parse_contact_weights(const std::string& text, ContactWeights& weights);

/// @brief The weights as "home=...,school=...,work=...,social_event=...,shop=...".
std::string contact_weights_to_string(const ContactWeights& weights);

void set_local_parameters(mio::abm::Model& model, CityParameters::Region region,
                          const ContactWeights& weights = ContactWeights{});

#endif // LOCATION_SPLIT_PARAMETER_SETTER_H
