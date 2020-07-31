#ifndef PARAMETER_STUDIES_H
#define PARAMETER_STUDIES_H

#include <epidemiology/parameter_studies/parameter_space.h>
#include <epidemiology_io/secir_result_io.h>
#include <iostream>
#include <sys/stat.h>
#include <unordered_map>

namespace epi
{
using HandleSimulationResultFunction = std::function<void(const epi::ContactFrequencyMatrix&, const epi::SecirParams&,
                                                          std::vector<double>, std::vector<Eigen::VectorXd>)>;

// The function type for the kind of simulation that we want to run
using secir_simulation_function_t =
    std::function<std::vector<double>(double t0, double tmax, double dt, ContactFrequencyMatrix const& cont_freq_matrix,
                                      SecirParams const& params, std::vector<Eigen::VectorXd>& secir_result)>;

// TODO: document class
// TODO: document input file convention

class ParameterStudy
{
public:
    /* 
     * @brief Constructor from file name
     * @param[in] parameter_filename filename of a file storing ranges of input parameters.
     */
    ParameterStudy(secir_simulation_function_t const& simu_func, ParameterSpace&& parameter_space, size_t n_runs,
                   double t0, double tmax);

    /* 
     * @brief Constructor from contact frequency matrix and parameter vector
     * @param[in] parameter_filename filename of a file storing ranges of input parameters.
     */
    ParameterStudy(secir_simulation_function_t const& simu_func, ContactFrequencyMatrix const& cont_freq_matrix,
                   SecirParams const& params, double t0, double tmax, double dev_rel = 0.2, size_t nb_runs = 1);

    /*
     * @brief Carry out all simulations in the parameter study.
     */
    std::vector<std::vector<Eigen::VectorXd>> run(HandleSimulationResultFunction simulation_result_function);

    /*
     * @brief sets the number of Monte Carlo runs
     * @param[in] nb_runs number of runs
     */
    void set_nb_runs(size_t nb_runs)
    {
        m_nb_runs = nb_runs;
    }

    /*
     * @brief returns the number of Monte Carlo runs
     */
    int get_nb_runs() const
    {
        return static_cast<int>(m_nb_runs);
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

private:
    // The path of the file storing the parameter ranges
    std::string parameter_file;

    // Stores the names and ranges of all parameters
    ParameterSpace parameter_space;

    // The function that carries out our simulation
    secir_simulation_function_t simulation_function;

    size_t m_nb_runs = 100;

    // Start time (should be the same for all simulations)
    double m_t0 = 0.0;
    // End time (should be the same for all simulations)
    double m_tmax = 400.0;
    // adaptive time step (will be corrected if too large/small)
    double m_dt = 0.1;
};

inline ParameterStudy::ParameterStudy(const secir_simulation_function_t& simu_func, ParameterSpace&& parameter_space,
                                      size_t n_runs, double t0, double tmax)
    : simulation_function(simu_func)
    , parameter_space(std::move(parameter_space))
    , m_nb_runs(n_runs)
    , m_t0{t0}
    , m_tmax{tmax}
{
}

inline ParameterStudy::ParameterStudy(secir_simulation_function_t const& simu_func,
                                      ContactFrequencyMatrix const& cont_freq_matrix, SecirParams const& params,
                                      double t0, double tmax, double dev_rel, size_t nb_runs)
    : simulation_function{simu_func}
    , parameter_space{cont_freq_matrix, params, t0, tmax, dev_rel}
    , m_nb_runs{nb_runs}
    , m_t0{t0}
    , m_tmax{tmax}
{
}

inline std::vector<std::vector<Eigen::VectorXd>>
ParameterStudy::run(HandleSimulationResultFunction simulation_result_function)
{
    std::vector<std::vector<Eigen::VectorXd>> ensemble_result;

    // Iterate over all parameters in the parameter space
    for (size_t i = 0; i < (*this).get_nb_runs(); i++) {

        std::vector<Eigen::VectorXd> secir_result;

        SecirParams params_sample             = parameter_space.get_secir_params_sample();
        ContactFrequencyMatrix contact_sample = parameter_space.get_cont_freq_matrix_sample();

        // Call the simulation function
        auto time = simulation_function((*this).m_t0, (*this).m_tmax, (*this).m_dt, contact_sample, params_sample,
                                        secir_result);
        simulation_result_function(contact_sample, params_sample, time, secir_result);

        ensemble_result.push_back(std::move(secir_result));
    }

    return ensemble_result;
}

} // namespace epi

#endif // PARAMETER_STUDIES_H
