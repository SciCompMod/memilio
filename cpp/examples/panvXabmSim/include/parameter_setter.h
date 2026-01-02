#pragma once
#include "abm/abm.h"

void set_parameters(mio::abm::Parameters& params);
void set_local_parameters_ger(mio::abm::World& world);
void set_local_parameters_usa(mio::abm::World& world);
void set_local_parameters_fra(mio::abm::World& world);
std::pair<double, double> get_my_and_sigma(std::pair<double, double> mean_and_std);
void set_local_parameters_event(mio::abm::World& world, double contact_rate_multiplier);