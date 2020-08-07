#ifndef READ_TWITTER_H
#define READ_TWITTER_H

#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Core>
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
int count_lines(const std::string& filename);

/**
 * @brief Reads raw Twitter Migration data ans writes it into a NxN Eigen Matrix, where N is the number of regions
 * @param filename name of file to be read
 */
Eigen::MatrixXi read_twitter(const std::string& filename);

} // namespace epi

#endif // READ_TWITTER_H
