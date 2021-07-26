#ifndef READ_TWITTER_H
#define READ_TWITTER_H

#include "epidemiology/utils/eigen.h"
#include "epidemiology/utils/io.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace epi
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

} // namespace epi

#endif // READ_TWITTER_H
