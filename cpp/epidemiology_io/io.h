#ifndef IO_H
#define IO_H

#include "epidemiology/utils/io.h"
#include <string>

namespace epi
{

/**
 * @brief Returns the current working directory name
 */ 
std::string get_current_dir_name();

/**
 * @brief Creates a directory on the file system
 * returns true if the creation was successful.
 */
IOResult<bool> create_directory(std::string const& rel_path, std::string& abs_path);

bool file_exists(std::string const& rel_path, std::string& abs_path);

} // namespace epi

#endif // IO_H
