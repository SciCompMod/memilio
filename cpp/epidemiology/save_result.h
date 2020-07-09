#ifndef SAVE_RESULT_H
#define SAVE_RESULT_H

#include <vector>
//#include <epidemiology/secir.h>
#include <epidemiology/eigen_util.h>

namespace epi
{

void save_result(const std::vector<double>& times, const std::vector<Eigen::VectorXd>& secir,
                 const std::string& filename);

struct result {
    std::vector<double> time;
    std::vector<std::vector<std::vector<double>>> groups;
    std::vector<std::vector<double>> total;
};

result read_result(const std::string& filename, int nb_groups);

} // namespace epi

#endif