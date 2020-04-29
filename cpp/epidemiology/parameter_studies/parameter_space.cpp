#include "parameter_space.h"
#ifndef NDEBUG
// TODO: If we have an own debug mode and assertion switch to this
#include <assert.h>
#endif

parameter_space_t::parameter_space_t(std::string &parameter_filename)
{
#ifndef NDEBUG
    // For debugging we check that we get the correct number
    // of parameters.
    const size_t num_params = 26;
#endif

    // TODO: Open file and read parameter names and values

#ifndef NDEBUG
    // For debugging we check that we get the correct number
    // of parameters.
    assert (parameter_names.size() == num_params);
#endif
}