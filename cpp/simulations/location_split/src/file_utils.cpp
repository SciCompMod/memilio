#include "file_utils.h"

#include "memilio/io/io.h"

#include <ctime>
#include <iterator>
#include <string>

std::string current_date_time()
{
    const std::time_t time = std::time({});
    char time_string[std::size("yyyy-mm-ddHHMMSS")];
    std::strftime(std::data(time_string), std::size(time_string), "%F%H%M%S", std::gmtime(&time));
    return std::string(time_string);
}

mio::IOResult<void> create_result_folders(const std::string& result_dir)
{
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir, /*create_parents=*/true));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir + "/infection_state_per_age_group"));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir + "/infection_per_location_type_per_age_group"));
    BOOST_OUTCOME_TRY(mio::create_directory(result_dir + "/total_infections"));
    return mio::success();
}
