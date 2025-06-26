#pragma once
#include "abm/abm.h"
#include <string>
#include <map>

std::map<uint32_t, bool> read_infection_data(const std::string& filename);
mio::abm::World create_world_from_file(const std::string& infection_data_file, int number_of_persons = 0);

template <typename DataFormat>
mio::abm::World read_world_from_file(const std::string& filename);