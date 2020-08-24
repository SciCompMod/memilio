#ifndef SAVE_RESULT_H
#define SAVE_RESULT_H

//#include <epidemiology/secir.h>
#include <epidemiology/eigen_util.h>
#include <epidemiology/time_series.h>

namespace epi
{
/**
 * @brief save secir simulation result to h5 file
 * @param times Vector of timesteps used during simulation
 * @param secir Results of secir simulation
 * @param filename name of file
 */
void save_result(const TimeSeries<double>& result, const std::string& filename);

class SecirSimulationResult
{
public:
    /**
     * @brief Standard constructor of SecirSimulationResult
     */
    SecirSimulationResult(int num_groups, int num_compartments)
        : m_groups(num_groups * num_compartments), m_totals(num_compartments)
    {}

    /**
     * @brief Constructor of SecirSimulationResult storing time, groups, and total sums of all groups
     * @param groups Simulation Results of individual groups
     * @param total Simulation Results of the sum over all groups
     */
    SecirSimulationResult(const TimeSeries<double>& groups, const TimeSeries<double>& totals)
        : m_groups(groups), m_totals(totals)
    {

    }

    /**
     * @brief Simulation Results of individual groups.
     */
    const TimeSeries<double>& get_groups() const
    {
        return m_groups;
    }

    /**
     * @brief Simulation Results of the sum over all groups.
     */
    const TimeSeries<double>& get_totals() const
    {
        return m_totals;
    }

private:
    TimeSeries<double> m_groups;
    TimeSeries<double> m_totals;
};

/**
 * @brief read secir simulation result from h5 file
 * @param filename name of file
 * @param nb_groups number of groups used during simulation
 */
SecirSimulationResult read_result(const std::string& filename, int nb_groups);

} // namespace epi

#endif
