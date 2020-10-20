#ifndef EPI_UTILS_EIGEN_H
#define EPI_UTILS_EIGEN_H

#include "epidemiology/utils/compiler_diagnostics.h"

/* this file wraps includes from eigen3 library to disable warnings. */

GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wint-in-bool-context")
GCC_CLANG_DIAGNOSTIC(ignored "-Wshadow")

#include <Eigen/Core>

GCC_CLANG_DIAGNOSTIC(pop)

#endif //EPI_UTILS_EIGEN_H