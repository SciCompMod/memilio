#ifndef SAVE_RESULT_H
#define SAVE_RESULT_H

#include <vector>
//#include <epidemiology/secir.h>
#include <epidemiology/eigen_util.h>

namespace epi
{
/**
 * @brief save secir simulation result to h5 file
 * @param times Vector of timesteps used during simulation
 * @param secir Results of secir simulation
 * @param filename name of file
 */
void save_result(const std::vector<double>& times, const std::vector<Eigen::VectorXd>& secir,
                 const std::string& filename);

/**
 * @brief strcut which holds secir simulation results
 * @param time Vector of timesteps used during simulation
 * @param groups Simulation Results of individual groups
 * @param total Simulation Results of the sum over all groups
 */
struct result {
    std::vector<double> time;
    std::vector<Eigen::VectorXd> groups;
    std::vector<std::vector<double>> total;
};

/**
 * @brief read secir simulation result from h5 file
 * @param filename name of file
 * @param nb_groups number of groups used during simulation
 */
result read_result(const std::string& filename, int nb_groups);

} // namespace epi

#endif
