#include "../include/file_utils.h"
#include "boost/filesystem.hpp"
#include "memilio/io/result_io.h"
#include <ctime>
#include <iostream>

namespace fs = boost::filesystem;

const std::string currentDateTime()
{
    std::time_t time = std::time({});
    char timeString[std::size("yyyy-mm-dddddddd")];
    std::strftime(std::data(timeString), std::size(timeString), "%F%H%M%S", std::gmtime(&time));
    return timeString;
}

mio::IOResult<bool> create_result_folders(std::string const& result_dir, int n_params)
{
    std::string inf_p_loc_t_p_ag = result_dir + "/infection_per_location_type_per_age_group/";
    std::string inf_state_p_ag = result_dir + "/infection_state_per_age_group/";

    BOOST_OUTCOME_TRY(mio::create_directory(result_dir));
    BOOST_OUTCOME_TRY(mio::create_directory(inf_p_loc_t_p_ag));
    BOOST_OUTCOME_TRY(mio::create_directory(inf_state_p_ag));
    
    if (n_params > 0) {
        for (int i = 0; i < n_params; i++) {
            BOOST_OUTCOME_TRY(mio::create_directory(inf_p_loc_t_p_ag + "/" + std::to_string(i)));
            BOOST_OUTCOME_TRY(mio::create_directory(inf_state_p_ag + "/" + std::to_string(i)));
        }
    }
    return mio::success();
}

mio::IOResult<bool> copy_result_folder(std::string const& from_dir, std::string const& to_dir)
{
    BOOST_OUTCOME_TRY(mio::create_directory(to_dir));
    fs::copy(from_dir, to_dir, fs::copy_options::overwrite_existing | fs::copy_options::recursive);
    return mio::success();
}