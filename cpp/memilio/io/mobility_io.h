/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele, Wadim Koslow, Martin J. Kuehn
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
#ifndef READ_TWITTER_H
#define READ_TWITTER_H

#include "memilio/math/eigen.h"
#include "memilio/io/io.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace mio
{

/**
 * @brief Splits string into a Vector of strings according to delimiter
 * @param s string which is splitted
 * @param delimiter sign at which to split string
 */
std::vector<std::string> split(const std::string& s, char delimiter);

/**
 * @brief Counts lines of txt file
 * @param filename name of file which is counted
 */
IOResult<int> count_lines(const std::string& filename);

/**
 * @brief Reads formatted migration or contact data which is given in columns
 *          from_str	to_str	from_rs	    to_rs	count_abs
 *        and separated by tabs. Writes it into a NxN Eigen Matrix, 
 *        where N is the number of regions
 * @param filename name of file to be read
 */
IOResult<Eigen::MatrixXd> read_mobility_formatted(const std::string& filename);

/**
 * @brief Reads txt migration data or contact which is given by values only
 *        and separated by spaces. Writes it into a NxN Eigen 
 *        Matrix, where N is the number of regions
 * @param filename name of file to be read
 */
IOResult<Eigen::MatrixXd> read_mobility_plain(const std::string& filename);

} // namespace mio

#endif // READ_TWITTER_H
