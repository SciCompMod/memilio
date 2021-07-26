#include <epidemiology_io/io.h>
#include "epidemiology/utils/compiler_diagnostics.h"
#include "epidemiology/utils/logging.h"


MSVC_WARNING_DISABLE_PUSH(4268)
#include <boost/filesystem.hpp>
#include <hdf5.h>
MSVC_WARNING_POP

namespace epi
{

std::string get_current_dir_name()
{
    auto path = boost::filesystem::current_path();
    return path.string();
}

IOResult<bool> create_directory(std::string const& rel_path, std::string& abs_path){
    boost::filesystem::path dir(rel_path);
    abs_path = dir.string();
    boost::system::error_code ec;
    bool created =  boost::filesystem::create_directory(dir, ec);
    if (ec) {
        return failure(ec, abs_path);
    }

    if (created) {
        log_info("Directory '{:s}' was created.", dir.string());
    }
    else {
        log_info(
            "Directory '{:s}' already exists.",
            dir.string(), epi::get_current_dir_name());
    }
    
    return success(created);
}

bool file_exists(std::string const& rel_path, std::string& abs_path)
{
    boost::filesystem::path dir(rel_path);
    abs_path = dir.string();
    return boost::filesystem::exists(dir);
}

} // namespace epi
