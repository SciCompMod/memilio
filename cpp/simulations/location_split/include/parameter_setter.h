#ifndef LOCATION_SPLIT_PARAMETER_SETTER_H
#define LOCATION_SPLIT_PARAMETER_SETTER_H

#include "city_parameters.h"

#include "abm/model.h"
#include "abm/parameters.h"

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
void set_local_parameters(mio::abm::Model& model, CityParameters::Region region);

#endif // LOCATION_SPLIT_PARAMETER_SETTER_H
