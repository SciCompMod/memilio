#ifndef LOCATION_SPLIT_FILE_UTILS_H
#define LOCATION_SPLIT_FILE_UTILS_H

#include "memilio/io/io.h"

#include <string>

/// @brief Current UTC date and time formatted as YYYY-MM-DDHHMMSS.
std::string current_date_time();

/// @brief Create the result directory and its subdirectories.
mio::IOResult<void> create_result_folders(const std::string& result_dir);

#endif // LOCATION_SPLIT_FILE_UTILS_H
