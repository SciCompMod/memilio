#ifndef PARAMETER_STUDIES_H
#define PARAMETER_STUDIES_H

#include <epidemiology/parameter_studies/parameter_space.h>
#include <epidemiology_io/secir_result_io.h>
#include <iostream>
#include <sys/stat.h>
#include <unordered_map>

namespace epi
{
using HandleSimulationResultFunction =
    std::function<void(const SecirParams&, const TimeSeries<double>& result)>;
auto dummy_func = [](const auto& params, const auto& secir_result) {};

// The function type for the kind of simulation that we want to run
using secir_simulation_function_t = std::function<TimeSeries<double>(double t0, double tmax, double dt, SecirParams const& params)>;

// TODO: document class

class ParameterStudy
{
public:
    /* 
     * @brief Constructor from file name
     * @param[in] parameter_filename filename of a file storing ranges of input parameters.
     */
    ParameterStudy(secir_simulation_function_t const& simu_func, ParameterSpace&& parameter_space, size_t num_runs,
                   double t0, double tmax);

    /* 
     * @brief Constructor from contact frequency matrix and parameter vector
     * @param[in] parameter_filename filename of a file storing ranges of input parameters.
     */
    ParameterStudy(secir_simulation_function_t const& simu_func, SecirParams const& params, double t0, double tmax,
                   size_t num_runs);

    /* 
     * @brief Constructor from contact frequency matrix and parameter vector
     * @param[in] parameter_filename filename of a file storing ranges of input parameters.
     */
    ParameterStudy(secir_simulation_function_t const& simu_func, SecirParams const& params, double t0, double tmax,
                   double dev_rel, size_t num_runs);

    /*
     * @brief Carry out all simulations in the parameter study.
     */
    std::vector<TimeSeries<double>> run(HandleSimulationResultFunction simulation_result_function = dummy_func);

    /*
     * @brief sets the number of Monte Carlo runs
     * @param[in] num_runs number of runs
     */

    void set_num_runs(size_t num_runs)
    {
        m_num_runs = num_runs;
    }

    /*
     * @brief returns the number of Monte Carlo runs
     */
    int get_num_runs() const
    {
        return static_cast<int>(m_num_runs);
    }

    /*
     * @brief sets end point in simulation
     * @param[in] tmax end point in simulation
     */
    void set_tmax(double tmax)
    {
        m_tmax = tmax;
    }

    /*
     * @brief returns end point in simulation
     */
    double get_tmax() const
    {
        return m_tmax;
    }

    void set_t0(double t0)
    {
        m_t0 = t0;
    }

    /*
     * @brief returns start point in simulation
     */
    double get_t0() const
    {
        return m_t0;
    }

    const ParameterSpace& get_parameter_space() const
    {
        return parameter_space;
    }

    ParameterSpace& get_parameter_space()
    {
        return parameter_space;
    }

private:
    // The path of the file storing the parameter ranges
    std::string parameter_file;

    // Stores the names and ranges of all parameters
    ParameterSpace parameter_space;

    // The function that carries out our simulation
    secir_simulation_function_t simulation_function;

    size_t m_num_runs = 100;

    // Start time (should be the same for all simulations)
    double m_t0 = 0.0;
    // End time (should be the same for all simulations)
    double m_tmax = 400.0;
    // adaptive time step (will be corrected if too large/small)
    double m_dt = 0.1;
};

inline ParameterStudy::ParameterStudy(const secir_simulation_function_t& simu_func, ParameterSpace&& parameter_space,
                                      size_t num_runs, double t0, double tmax)
    : simulation_function(simu_func)
    , parameter_space(std::move(parameter_space))
    , m_num_runs(num_runs)
    , m_t0{t0}
    , m_tmax{tmax}
{
}

inline ParameterStudy::ParameterStudy(secir_simulation_function_t const& simu_func, SecirParams const& params,
                                      double t0, double tmax, size_t num_runs)
    : simulation_function{simu_func}
    , parameter_space{params}
    , m_num_runs{num_runs}
    , m_t0{t0}
    , m_tmax{tmax}
{
}

inline ParameterStudy::ParameterStudy(secir_simulation_function_t const& simu_func, SecirParams const& params,
                                      double t0, double tmax, double dev_rel, size_t num_runs)
    : simulation_function{simu_func}
    , parameter_space{params, t0, tmax, dev_rel}
    , m_num_runs{num_runs}
    , m_t0{t0}
    , m_tmax{tmax}
{
}

inline std::vector<TimeSeries<double>> ParameterStudy::run(HandleSimulationResultFunction simulation_result_function)
{
    std::vector<TimeSeries<double>> ensemble_result;

    // Iterate over all parameters in the parameter space
    for (size_t i = 0; i < (*this).get_num_runs(); i++) {
        SecirParams params_sample = parameter_space.draw_sample();

        // Call the simulation function
        auto result = simulation_function((*this).m_t0, (*this).m_tmax, (*this).m_dt, params_sample);
        simulation_result_function(params_sample, result);

        ensemble_result.push_back(result);
    }

    return ensemble_result;
}

} // namespace epi

#endif // PARAMETER_STUDIES_H
