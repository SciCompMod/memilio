#include <epidemiology_io/read_twitter.h>

#include <iostream>
#include <fstream>
#include <string>

#include <Eigen/Core>
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
    if (file.is_open()) { //checking whether the file is open
        std::string line;
        while (getline(file, line)) {
            count++;
        }
        file.close();
    }
    return count;
}

Eigen::MatrixXi read_twitter(const std::string& filename)
{
    std::fstream file;

    int nb_lines = count_lines(filename);

    file.open(filename, std::ios::in);

    Eigen::MatrixXi txtMatrix(nb_lines - 1, 3);
    std::vector<int> ids;

    std::string tp;
    if (file.is_open()) {
        int i = 0;
        while (getline(file, tp)) {
            auto line = split(tp, '\t');
            if (i > 0) {
                ids.push_back(stoi(line[2]));
                txtMatrix(i - 1, 0) = std::stoi(line[2]);
                txtMatrix(i - 1, 1) = std::stoi(line[3]);
                txtMatrix(i - 1, 2) = std::stoi(line[4]);
            }
            i++;
        }
    }
    sort(ids.begin(), ids.end());

    std::vector<int>::iterator iter = std::unique(ids.begin(), ids.end());
    ids.resize(std::distance(ids.begin(), iter));

    Eigen::MatrixXi twitter(ids.size(), ids.size());
    for (int i = 0; i < ids.size(); i++) {
        for (int j = 0; j < ids.size(); j++)
            twitter(i, j) = 0;
    }

    for (int k = 0; k < nb_lines - 1; k++) {
        int row_ind = 0;
        int col_ind = 0;
        while (txtMatrix(k, 0) != ids[row_ind]) {
            row_ind++;
        }
        while (txtMatrix(k, 1) != ids[col_ind]) {
            col_ind++;
        }
        twitter(row_ind, col_ind) = txtMatrix(k, 2);
    }
    file.close();

    return twitter;
}

} // namespace epi
