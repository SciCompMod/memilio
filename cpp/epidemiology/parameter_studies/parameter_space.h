#ifndef PARAMETER_SPACE_H
#define PARAMETER_SPACE_H

#include <vector>
#include <string>

/* The class parameter_space_t stores ranges of parameters
 * together with information on step sizes,
 * a start and end time as well as an initial time step.
 * The class provides an iterator that iterates over all 
 * generated parameter combinations.
 * 
 * Currently all parameters are of type double.
 */
class parameter_space_t
{
public:
    using ConstIterator = int; // TODO: placeholder. Replace with correct iterator once it is implemented

    /* Constructor
     * \param [in] paramter_filename filename of a file storing ranges of input parameters.
     * Reads parameter names and values from an input file.
     */
    parameter_space_t (std::string &parameter_filename);

    // Return an iterator that starts at the first parameter set
    ConstIterator Begin();
private:
    // The names of all stored parameters
    std::vector<std::string> parameter_names;
    
    // The start values of the parameters
    std::vector<double> parameter_start_values;
    // The end values of the parameters
    std::vector<double> parameter_end_values;
    // The step size values of the parameters
    std::vector<double> parameter_step_values;

    // Start time (should be the same for all simulations)
    T t0;
    // End time (should be the same for all simulations)
    T tmax;
    // time step (should be the same for all simulations)
    T dt;
}


#endif // PARAMETER_SPACE_H