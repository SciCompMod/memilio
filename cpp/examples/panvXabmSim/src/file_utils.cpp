#include "../include/file_utils.h"
#include "boost/filesystem.hpp"
#include "memilio/io/result_io.h"
#include <ctime>
#include <iostream>
#include <string>

namespace fs = boost::filesystem;

std::string currentDateTime()
{
    std::time_t time = std::time({});
    char timeString[std::size("yyyy-mm-ddHHMMSS")];
    std::strftime(std::data(timeString), std::size(timeString), "%F%H%M%S", std::gmtime(&time));
    return std::string(timeString);
}

mio::IOResult<void> create_result_folders(const std::string& result_dir)
{
    // Define subdirectory paths
    const std::string inf_per_location_type_per_age_group = result_dir + "/infection_per_location_type_per_age_group";
    const std::string inf_state_per_age_group             = result_dir + "/infection_state_per_age_group";

    // Create main directories
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir));
    BOOST_OUTCOME_TRY(mio::create_directory(inf_per_location_type_per_age_group));
    BOOST_OUTCOME_TRY(mio::create_directory(inf_state_per_age_group));

    return mio::success();
}

mio::IOResult<void> copy_result_folder(const std::string& from_dir, const std::string& to_dir)
{
    // Create destination directory
    BOOST_OUTCOME_TRY(mio::create_directory(to_dir));

    try {
        // Copy directory recursively with overwrite
        fs::copy(from_dir, to_dir, fs::copy_options::overwrite_existing | fs::copy_options::recursive);
    }
    catch (const fs::filesystem_error& e) {
        return mio::failure(mio::StatusCode::InvalidFileFormat, "Failed to copy directory: " + std::string(e.what()));
    }

    return mio::success();
}