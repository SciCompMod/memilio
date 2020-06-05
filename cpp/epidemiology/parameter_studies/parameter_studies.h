#ifndef PARAMETER_STUDIES_H
#define PARAMETER_STUDIES_H

#include <epidemiology/parameter_studies/parameter_space.h>
#include <iostream>
#include <unordered_map>

namespace epi
{

// The function type for the kind of simulation that we want to run
using secir_simulation_function_t = void (*)(const double t0, const double tmax, const double dt,
                                             SecirParams const& params, std::vector<std::vector<double>>& seir);

// TODO: document class
// TODO: document input file convention

class parameter_study_t
{
public:
    /* Constructor
   * \param [in] paramter_filename filename of a file storing ranges of input
   * parameters.
   */
    parameter_study_t(std::string& parameter_filename);

    // Carry out all simulations in the parameter study.
    void run();

    /*
     * @brief sets end point in simulation and 
     * limits the maximum number of dampings by one per 10 days (in the mean)
     * @param[in] tmax end point in simulation
     */
    void set_tmax(double tmax)
    {
        m_tmax            = tmax;
        m_max_nb_dampings = (int)(tmax / 10);
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

    // Start time (should be the same for all simulations)
    double m_t0 = 0;
    // End time (should be the same for all simulations)
    double m_tmax = 400;
    // adaptive time step (will be corrected if too large/small)
    double dt = 0.1;
    // maximum number of dampings; to avoid overfitting only allow one damping for every 10 days simulated
    int m_max_nb_dampings = 40;
};

parameter_study_t::parameter_study_t(std::string& parameter_filename)
    : parameter_space(parameter_filename)
{
}

void parameter_study_t::run()
{
// Iterate over all parameters in the parameter space
#if 0
    for(parameter_space_t::ConstIterator param_it = parameter_space.BeginConst(iCell); param_it.IsValid(); ++param_it){
        // Get the current parameters
        const struct seirParam<double> &params = *param_it;

        // Print the parameters if we are in debug mode
        // TODO: Should we get an own debugging mode use it instead of NDEBUG
        // TODO: Replace cout with logging function once we have one.
#ifndef NDEBUG
        std::cout << "Starting simulation with params:\n" << std::endl;
        printSeirParams (params);
#endif
        // The vector in which we store the result.
        /* TODO: In the current version we do not use the result.
         *       This is of course only temporarily until we have a 
         *       mechanism to collect the results.
         */
        std::vector<std::vector<T>> result_vector;
        // Call the simulation function
        simulation_function (paramter_space.t0, paramter_space.tmax, paramter_space.dt, params, result_vector);
    }
#endif
}

} // namespace epi

#endif // PARAMETER_STUDIES_H