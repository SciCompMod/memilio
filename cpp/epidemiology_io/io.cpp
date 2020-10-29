#include <epidemiology_io/io.h>
#include "epidemiology/utils/compiler_diagnostics.h"


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

} // namespace epi
