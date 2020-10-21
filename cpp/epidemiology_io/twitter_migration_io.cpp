#include <epidemiology_io/twitter_migration_io.h>
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

int count_lines(const std::string& filename)
{
    std::fstream file;
    file.open(filename, std::ios::in);

    int count = 0;
    if (file.is_open()) { //checking whether the file could be opened
        std::string line;
        while (getline(file, line)) {
            count++;
        }
        file.close();
    }
    else {
        epi::log_error("File could not be opened.");
    }
    return count;
}

Eigen::MatrixXi read_migration(const std::string& filename)
{
    std::fstream file;

    int num_lines = count_lines(filename);

    if (num_lines > 0) {
        file.open(filename, std::ios::in);

        Eigen::MatrixXi txt_matrix(num_lines - 1, 3);
        std::vector<int> ids;

        std::string tp;
        if (file.is_open()) { //checking whether the file could be opened
            int linenumber = 0;
            while (getline(file, tp)) {
                auto line = split(tp, '\t');
                if (linenumber > 0) {
                    ids.push_back(std::stoi(line[2]));
                    txt_matrix(linenumber - 1, 0) = std::stoi(line[2]);
                    txt_matrix(linenumber - 1, 1) = std::stoi(line[3]);
                    txt_matrix(linenumber - 1, 2) = std::stoi(line[4]);
                }
                linenumber++;
            }

            sort(ids.begin(), ids.end());

            std::vector<int>::iterator iter = std::unique(ids.begin(), ids.end());
            ids.resize(std::distance(ids.begin(), iter));

            Eigen::MatrixXi migration(ids.size(), ids.size());
            for (Eigen::Index i = 0; i < static_cast<Eigen::Index>(ids.size()); i++) {
                for (Eigen::Index j = 0; j < static_cast<Eigen::Index>(ids.size()); j++)
                    migration(i, j) = 0;
            }

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
            file.close();

            return migration;
        }
        else {
            epi::log_error("File could not be opened.");

            return Eigen::MatrixXi(0, 0);
        }
    }
    else {
        return Eigen::MatrixXi(0, 0);
    }
}

} // namespace epi
