#pragma once
#include "memilio/io/result_io.h"
#include <string>

// Get current date and time as a string in format YYYY-MM-DDHHMMSS
std::string currentDateTime();

// Create result directory structure for simulation outputs
mio::IOResult<void> create_result_folders(const std::string& result_dir);

// Copy a directory and its contents recursively
mio::IOResult<void> copy_result_folder(const fs::path& from_dir, const fs::path& to_dir);