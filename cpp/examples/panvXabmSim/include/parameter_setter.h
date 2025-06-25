#pragma once
#include "abm/abm.h"

void set_parameters(mio::abm::Parameters& params);
void set_local_parameters(mio::abm::World& world);
std::pair<double, double> get_my_and_sigma(std::pair<double, double> mean_and_std);