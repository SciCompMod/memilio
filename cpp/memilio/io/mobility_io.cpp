#include "memilio/io/mobility_io.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace mio
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
    BOOST_OUTCOME_TRY(auto&& num_lines, count_lines(filename));

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

    Eigen::MatrixXd mobility = Eigen::MatrixXd::Zero(ids.size(), ids.size());

    for (int k = 0; k < num_lines - 1; k++) {
        int row_ind = 0;
        int col_ind = 0;
        while (txt_matrix(k, 0) != ids[row_ind]) {
            row_ind++;
        }
        while (txt_matrix(k, 1) != ids[col_ind]) {
            col_ind++;
        }
        mobility(row_ind, col_ind) = txt_matrix(k, 2);
    }

    return success(mobility);
}

IOResult<Eigen::MatrixXd> read_mobility_plain(const std::string& filename)
{
    BOOST_OUTCOME_TRY(auto&& num_lines, count_lines(filename));

    if (num_lines == 0) {
        return success(Eigen::MatrixXd(0, 0));
    }

    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, filename);
    }

    Eigen::MatrixXd mobility(num_lines, num_lines);

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
                mobility(i, j) = std::stod(line[j]);
            }
            linenumber++;
        }
    }
    catch (std::runtime_error& ex) {
        return failure(StatusCode::InvalidFileFormat, filename + ": " + ex.what());
    }

    return success(mobility);
}

IOResult<Eigen::MatrixXd> read_duration_stay(const std::string& filename)
{
    BOOST_OUTCOME_TRY(auto&& num_lines, count_lines(filename));

    if (num_lines == 0) {
        return success(Eigen::MatrixXd(0, 0));
    }

    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, filename);
    }

    Eigen::VectorXd duration(num_lines);

    try {
        std::string tp;
        int linenumber = 0;
        while (getline(file, tp)) {
            auto line      = split(tp, ' ');
            Eigen::Index i = static_cast<Eigen::Index>(linenumber);
            duration(i)    = std::stod(line[0]);
            linenumber++;
        }
    }
    catch (std::runtime_error& ex) {
        return failure(StatusCode::InvalidFileFormat, filename + ": " + ex.what());
    }

    return success(duration);
}

IOResult<std::vector<std::vector<std::vector<int>>>> read_path_mobility(const std::string& filename)
{
    BOOST_OUTCOME_TRY(auto&& num_lines, count_lines(filename));

    if (num_lines == 0) {
        std::vector<std::vector<std::vector<int>>> arr(0, std::vector<std::vector<int>>(0));
        return success(arr);
    }

    std::fstream file;
    file.open(filename, std::ios::in);
    if (!file.is_open()) {
        return failure(StatusCode::FileNotFound, filename);
    }

    const int num_nodes = static_cast<int>(std::sqrt(num_lines));
    std::vector<std::vector<std::vector<int>>> arr(num_nodes, std::vector<std::vector<int>>(num_nodes));

    try {
        std::string tp;
        while (getline(file, tp)) {
            auto line   = split(tp, ' ');
            auto indx_x = std::stoi(line[0]);
            auto indx_y = std::stoi(line[1]);
            if (indx_x != indx_y) {
                auto path = std::accumulate(line.begin() + 2, line.end(), std::string(""));

                // string -> vector of integers
                std::vector<int> path_vec;

                // Remove the square brackets and \r
                path = path.substr(1, path.size() - 3);
                std::stringstream ss(path);
                std::string token;

                // get numbers and save them in path_vec
                while (std::getline(ss, token, ',')) {
                    path_vec.push_back(std::stoi(token));
                }

                for (int number : path_vec) {
                    arr[indx_x][indx_y].push_back(number);
                }
            }
            else {
                arr[indx_x][indx_y].push_back(static_cast<int>(indx_x));
            }
        }
    }
    catch (std::runtime_error& ex) {
        return failure(StatusCode::InvalidFileFormat, filename + ": " + ex.what());
    }

    return success(arr);
}

} // namespace mio
