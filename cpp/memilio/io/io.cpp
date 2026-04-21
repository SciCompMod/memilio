/* 
* Copyright (C) 2020-2026 MEmilio
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

#include <filesystem>

namespace mio
{

std::string get_current_dir_name()
{
    auto path = std::filesystem::current_path();
    return path.string();
}

IOResult<bool> create_directory(const std::filesystem::path& rel_path, std::string& abs_path, bool create_parents)
{
    auto result = create_directory(rel_path, create_parents);
    if (result) {
        std::error_code ec;
        abs_path = std::filesystem::canonical(rel_path, ec).string();
        if (ec) {
            return failure(ec, "Failed to get absolute path of " + rel_path.string());
        }
    }
    return result;
}

IOResult<bool> create_directory(const std::filesystem::path& rel_path, bool create_parents)
{
    std::error_code ec;
    bool created;

    if (create_parents) {
        created = std::filesystem::create_directories(rel_path, ec);
    }
    else {
        created = std::filesystem::create_directory(rel_path, ec);
    }

    if (ec) {
        const std::string with_parents = create_parents ? " (with parents)" : "";
        return failure(ec, "Failed to create directory " + rel_path.string() + with_parents);
    }

    if (created) {
        log_info("Directory '{:s}' was created.", rel_path);
    }
    else {
        log_info("Directory '{:s}' already exists.", rel_path);
    }

    return success(created);
}

std::filesystem::path create_directories_or_exit(const std::filesystem::path& path, bool create_parents)
{
    std::string abs_path;
    auto result = create_directory(path, abs_path, create_parents);
    if (!result) {
        log_critical("Could not create directory \"{}\": {}", path, result.error().message());
        exit(result.error().code().value());
    }
    return abs_path;
}

bool file_exists(std::string const& rel_path, std::string& abs_path)
{
    std::filesystem::path dir(rel_path);
    abs_path = dir.string();
    return std::filesystem::exists(dir);
}

} // namespace mio
