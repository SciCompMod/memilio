#ifndef PARAMETER_STUDIES_H
#define PARAMETER_STUDIES_H

#include <iostream>

// The function type for the kind of simulation that we want to run
template<typename T>
using seir_simulation_function_t = void (*)(const T t0, const T tmax, const T dt, struct seirParam<T> const& params, std::vector<std::vector<T> >& seir)

// TODO: document class
// TODO: document input file convention

// template parameter is the datatype of the associated seirParams
template <typename T>
class parameter_study_t
{
public:
    /* Constructor
     * \param [in] paramter_filename filename of a file storing ranges of input parameters.
     */
    parameter_study_t (std::string &parameter_filename);

    // Carry out all simulations in the parameter study.
    void run();
private:
    // The path of the file storing the parameter ranges
    std::string parameter_file;

    // Stores the names and ranges of all parameters
    parameter_space_t parameter_space;

    // The function that carries out our simulation
    seir_simulation_function_t simulation_function;

    // Start time (should be the same for all simulations)
    T t0;
    // End time (should be the same for all simulations)
    T tmax;
    // time step (should be the same for all simulations)
    T dt;
}

template <typename T>
void parameter_study_t::run ()
{
    // Iterate over all parameters in the parameter space
    for(parameter_space_t::ConstIterator param_it = parameter_space.BeginConst(iCell); param_it.IsValid(); ++param_it){
        // Get the current parameters
        const struct seirParam<T> &params = *param_it;

        // Print the parameters if we are in debug mode
        // TODO: Should we get an own debugging mode use it instead of NDEBUG
        // TODO: Replace cout with logging function once we have one.
#ifndef NDEBUG
        std::cout << "Starting simulation with params:\n" << std::end;
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
}

#endif // PARAMETER_STUDIES_H