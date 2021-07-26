#include <epidemiology_io/mobility_io.h>
#include "epidemiology/utils/eigen.h"
#include <epidemiology/utils/logging.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace epi
{

std::vector<std::string> split(const std::string& s, char delimitor)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimitor)) {
        tokens.push_back(token);
    }
    return tokens;
}

IOResult<int> count_lines(const std::string& filename)
{
    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, filename);
    }

    int count = 0;
    std::string line;
    while (getline(file, line)) {
        count++;
    }
    return success(count);
}

IOResult<Eigen::MatrixXd> read_mobility_formatted(const std::string& filename)
{
    BOOST_OUTCOME_TRY(num_lines, count_lines(filename));

    if (num_lines == 0) {
        return success(Eigen::MatrixXd(0, 0));
    }

    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, filename);
    }

    Eigen::MatrixXd txt_matrix(num_lines - 1, 3);
    std::vector<int> ids;

    try {
        std::string tp;
        int linenumber = 0;
        while (getline(file, tp)) {
            if (linenumber > 0) {
                auto line = split(tp, '\t');
                if (line.size() < 5) {
                    return failure(StatusCode::InvalidFileFormat,
                                   filename + ":" + std::to_string(linenumber) + ": Not enough entries in line.");
                }
                ids.push_back(std::stoi(line[2]));
                txt_matrix(linenumber - 1, 0) = std::stoi(line[2]);
                txt_matrix(linenumber - 1, 1) = std::stoi(line[3]);
                txt_matrix(linenumber - 1, 2) = std::stod(line[4]);
            }
            linenumber++;
        }
    }
    catch (std::runtime_error& ex) {
        return failure(StatusCode::InvalidFileFormat, filename + ": " + ex.what());
    }

    sort(ids.begin(), ids.end());
    std::vector<int>::iterator iter = std::unique(ids.begin(), ids.end());
    ids.resize(std::distance(ids.begin(), iter));

    Eigen::MatrixXd migration = Eigen::MatrixXd::Zero(ids.size(), ids.size());

    for (int k = 0; k < num_lines - 1; k++) {
        int row_ind = 0;
        int col_ind = 0;
        while (txt_matrix(k, 0) != ids[row_ind]) {
            row_ind++;
        }
        while (txt_matrix(k, 1) != ids[col_ind]) {
            col_ind++;
        }
        migration(row_ind, col_ind) = txt_matrix(k, 2);
    }

    return success(migration);
}

IOResult<Eigen::MatrixXd> read_mobility_plain(const std::string& filename)
{
    BOOST_OUTCOME_TRY(num_lines, count_lines(filename));

    if (num_lines == 0) {
        return success(Eigen::MatrixXd(0, 0));
    }

    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, filename);
    }

    Eigen::MatrixXd migration(num_lines, num_lines);

    try {
        std::string tp;
        int linenumber = 0;
        while (getline(file, tp)) {
            auto line = split(tp, ' ');
            if (line.size() != size_t(num_lines)) {
                return failure(StatusCode::InvalidFileFormat, filename + ": Not a square matrix.");
            }
            Eigen::Index i = static_cast<Eigen::Index>(linenumber);
            for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(line.size()); j++) {
                migration(i, j) = std::stod(line[j]);
            }
            linenumber++;
        }
    }
    catch (std::runtime_error& ex) {
        return failure(StatusCode::InvalidFileFormat, filename + ": " + ex.what());
    }

    return success(migration);
}

} // namespace epi
