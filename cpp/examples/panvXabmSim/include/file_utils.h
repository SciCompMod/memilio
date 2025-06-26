#pragma once
#include "memilio/io/result_io.h"
#include <string>

// Get current date and time as a string in format YYYY-MM-DDHHMMSS
std::string currentDateTime();

// Create result directory structure for simulation outputs
mio::IOResult<void> create_result_folders(const std::string& result_dir, int n_params = 1);

// Copy a directory and its contents recursively
mio::IOResult<void> copy_result_folder(const std::string& from_dir, const std::string& to_dir);