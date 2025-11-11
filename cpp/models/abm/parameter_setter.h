#pragma once
#include "model.h"

void set_parameters(mio::abm::Parameters& params);
void set_local_parameters(mio::abm::Model& world);
std::pair<double, double> get_my_and_sigma(std::pair<double, double> mean_and_std);
void set_local_parameters_event(mio::abm::Model& world, double contact_rate_multiplier);