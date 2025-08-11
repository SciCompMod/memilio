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
#include "test_data_dir.h"
#include "memilio/utils/stl_util.h"
#include <string>
#include <cstring>
#include <vector>
#include <fstream>

template <class String>
std::string get_test_data_file_path(String&& filename)
{
    return mio::path_join(TEST_DATA_DIR, filename);
}

/**
 * @brief read table of space separated numerical values from a file
 */
template <class Real>
std::vector<std::vector<Real>> load_test_data_csv(const std::string& filename)
{
    // File pointer
    std::fstream fin;

    // Open an existing file
    fin.open(get_test_data_file_path(filename), std::ios::in);

    // Read the Data from the file
    std::vector<std::vector<Real>> data;
    std::vector<Real> row;
    std::string line, word;
    while (getline(fin, line)) {
        row.clear();

        // ignore comments
        if (line[0] == '#') {
            continue;
        }

        //read columns in this row
        std::stringstream s(line);
        double v;
        while (s >> v) {
            row.push_back(v);
        }

        //append non empty rows
        if (row.size() > 0) {
            data.push_back(row);
        }
    }

    return data;
}
