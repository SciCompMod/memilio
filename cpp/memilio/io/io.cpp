/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Wadim Koslow
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include "memilio/io/io.h"
#include "memilio/utils/compiler_diagnostics.h"
#include "memilio/utils/logging.h"

MSVC_WARNING_DISABLE_PUSH(4268)
#include <boost/filesystem.hpp>
MSVC_WARNING_POP()

namespace mio
{

std::string get_current_dir_name()
{
    auto path = boost::filesystem::current_path();
    return path.string();
}

IOResult<bool> create_directory(std::string const& rel_path, std::string& abs_path)
{
    boost::filesystem::path dir(rel_path);
    boost::system::error_code ec;
    bool created = boost::filesystem::create_directory(dir, ec);
    if (ec) {
        return failure(ec, rel_path);
    }
    abs_path = boost::filesystem::canonical(dir, ec).string();
    if (ec) {
        return failure(ec, rel_path);
    }

    if (created) {
        log_info("Directory '{:s}' was created.", dir.string());
    }
    else {
        log_info("Directory '{:s}' already exists.", dir.string(), mio::get_current_dir_name());
    }

    return success(created);
}

IOResult<bool> create_directory(std::string const& rel_path)
{
    std::string abs_path;
    return create_directory(rel_path, abs_path);
}

bool file_exists(std::string const& rel_path, std::string& abs_path)
{
    boost::filesystem::path dir(rel_path);
    abs_path = dir.string();
    return boost::filesystem::exists(dir);
}

} // namespace mio
