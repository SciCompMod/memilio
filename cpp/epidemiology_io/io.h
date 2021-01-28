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

/**
 * @brief Creates a directory on the file system
 * returns true if the creation was successful.
 */
bool create_directory(std::string const& rel_path, std::string& abs_path);

bool directory_exists(std::string const& rel_path, std::string& abs_path);

} // namespace epi

#endif // IO_H
