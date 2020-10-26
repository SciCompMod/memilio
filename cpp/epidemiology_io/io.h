#ifndef IO_H
#define IO_H

#include <string>

namespace epi
{

/**
 * @brief Function to clean up all memory allocated by tixi and hdf5
 */
void io_cleanup();

/**
 * @brief Returns the current working directory name
 */ 
std::string get_current_dir_name();

} // namespace epi

#endif // IO_H
