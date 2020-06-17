#ifndef PARAMETER_STUDIES_H
#define PARAMETER_STUDIES_H

#include <epidemiology/parameter_studies/parameter_space.h>
#include <iostream>
#include <unordered_map>

namespace epi
{

// The function type for the kind of simulation that we want to run
using secir_simulation_function_t = void (*)(const double t0, const double tmax, const double dt,
                                             SecirParams const& params, std::vector<std::vector<double>>& secir);

// TODO: document class
// TODO: document input file convention

class parameter_study_t
{
public:
    /* 
     * @brief Constructor from file name
     * @param[in] parameter_filename filename of a file storing ranges of input parameters.
     */
    parameter_study_t(std::string& parameter_filename);

    /* 
     * @brief Constructor from contact frequency matrix and parameter vector
     * @param[in] parameter_filename filename of a file storing ranges of input parameters.
     */
    parameter_study_t(ContactFrequencyMatrix const& cont_freq_matrix, std::vector<SecirParams> const& params, double t0,
                      double tmax, double dev_rel = 0.2, size_t nb_runs = 1);

    /*
     * @brief Carry out all simulations in the parameter study.
     */
    void run();

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
        return m_nb_runs;
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

private:
    // The path of the file storing the parameter ranges
    std::string parameter_file;

    // Stores the names and ranges of all parameters
    parameter_space_t parameter_space;

    // The function that carries out our simulation
    secir_simulation_function_t simulation_function;

    size_t m_nb_runs = 100;

    // Start time (should be the same for all simulations)
    double m_t0 = 0;
    // End time (should be the same for all simulations)
    double m_tmax = 400;
    // adaptive time step (will be corrected if too large/small)
    double dt = 0.1;
};

parameter_study_t::parameter_study_t(std::string& parameter_filename)
    : parameter_space(parameter_filename)
{
}

parameter_study_t::parameter_study_t(ContactFrequencyMatrix const& cont_freq_matrix,
                                     std::vector<SecirParams> const& params, double t0, double tmax, double dev_rel,
                                     size_t nb_runs)
    : parameter_space(cont_freq_matrix, params, t0, tmax, dev_rel)
    , m_nb_runs{nb_runs}
{
}

void parameter_study_t::run()
{
    // parameter_space.get_parameters().at("incubation time")->

    // Iterate over all parameters in the parameter space
    for (size_t i = 0; i < (*this).get_nb_runs(); i++) {
        // Get the current parameters
        // const struct seirParam<double>& params = *param_it;

        // Print the parameters if we are in debug mode
        // TODO: Should we get an own debugging mode use it instead of NDEBUG
        // TODO: Replace cout with logging function once we have one.
#ifndef NDEBUG
        std::cout << "Starting simulation with params:\n" << std::endl;
        // epi::print_secir_params(params, cont_freq);
#endif
        // The vector in which we store the result.
        /* TODO: In the current version we do not use the result.
         *       This is of course only temporarily until we have a 
         *       mechanism to collect the results.
         */
        // std::vector<std::vector<T>> result_vector;
        // Call the simulation function
        // simulation_function(paramter_space.t0, paramter_space.tmax, paramter_space.dt, params, result_vector);
    }
}

} // namespace epi

#endif // PARAMETER_STUDIES_H