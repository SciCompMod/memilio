#include <epidemiology_io/io.h>
#include "epidemiology/utils/compiler_diagnostics.h"
#include "epidemiology/utils/logging.h"


MSVC_WARNING_DISABLE_PUSH(4268)
#include <boost/filesystem.hpp>
#include <tixi.h>
#include <hdf5.h>
MSVC_WARNING_POP

namespace epi
{

void io_cleanup()
{
    tixiCleanup();
    H5close();
}

std::string get_current_dir_name()
{
    auto path = boost::filesystem::current_path();
    return path.string();
}

bool create_directory(std::string const& rel_path, std::string& abs_path){
    boost::filesystem::path dir(rel_path);
    abs_path = dir.string();
    bool created =  boost::filesystem::create_directory(dir);
    if (created) {
        log_info("Directory '{:s}' was created.", dir.string());
    }
    else {
        log_info(
            "Directory '{:s}' already exists.",
            dir.string(), epi::get_current_dir_name());
    }
    return created;
}

bool directory_exists(std::string const& rel_path, std::string& abs_path)
{
    boost::filesystem::path dir(rel_path);
    abs_path = dir.string();
    return boost::filesystem::exists(dir);
}

} // namespace epi
