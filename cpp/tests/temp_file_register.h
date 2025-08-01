/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele
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
#ifndef EPI_TEST_TMP_FILE_REGISTER_H
#define EPI_TEST_TMP_FILE_REGISTER_H

#include "memilio/io/io.h"
#include "memilio/utils/logging.h"
#include "boost/filesystem.hpp"

/**
 * Stores paths of files or directories that will be removed when the register is destroyed.
 */
class TempFileRegister
{
public:
    /**
     * Destructor.
     * Removes all registered files.
     */
    ~TempFileRegister()
    {
        for (auto&& file : m_files) {
            boost::system::error_code ec;
            boost::filesystem::remove_all(file, ec);
            if (ec) {
                //just log a warning, failed cleanup should not be considered a test failure.
                mio::log_warning("Failed to remove temporary file {}:{}", file.string(), ec.message());
            }
        }
    }

    /**
     * create a unique path and register it for cleanup.
     * The path will be in the system temp directory if available.
     * Otherwise it will be in the current working directory.
     * The name of the file or directory follows the specified model.
     * The model may contain placeholders one or more `%` that are replaced with random characters to make the name unique.
     * @param model name of the file or directory with placeholders for random characters.
     * @return a unique file path that follows the specified model.
     */
    std::string get_unique_path(const std::string& model = "%%%%-%%%%-%%%%-%%%%")
    {
        auto tmp_path  = get_tmp_path();
        auto file_name = boost::filesystem::unique_path(model);
        auto file_path = tmp_path / file_name;
        m_files.push_back(file_path);
        return file_path.string();
    }

private:
    //get the system temp directory if available, otherwise the current working directory
    boost::filesystem::path get_tmp_path()
    {
        boost::system::error_code ec;
        auto path = boost::filesystem::temp_directory_path(ec);
        if (ec) {
            path = boost::filesystem::current_path();
        }
        return path;
    }

    std::vector<boost::filesystem::path> m_files;
};

#endif //EPI_TEST_TMP_FILE_REGISTER_H
