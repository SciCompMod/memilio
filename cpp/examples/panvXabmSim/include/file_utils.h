#pragma once
#include "memilio/io/result_io.h"
#include <string>

const std::string currentDateTime();
mio::IOResult<bool> create_result_folders(std::string const& result_dir, int n_params = 1);
mio::IOResult<bool> copy_result_folder(std::string const& from_dir, std::string const& to_dir);