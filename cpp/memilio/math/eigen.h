/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Daniel Abele
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef EPI_UTILS_EIGEN_H
#define EPI_UTILS_EIGEN_H

#include "memilio/utils/compiler_diagnostics.h"

/* this file wraps includes from eigen3 library to disable warnings. */

//C4996: some std functions that have been deprecated in c++17; maybe fixed in new eigen versions? 
MSVC_WARNING_DISABLE_PUSH(4996)

GCC_CLANG_DIAGNOSTIC(push)
GCC_CLANG_DIAGNOSTIC(ignored "-Wint-in-bool-context")
GCC_CLANG_DIAGNOSTIC(ignored "-Wshadow")

#include <Eigen/Core>

GCC_CLANG_DIAGNOSTIC(pop)

MSVC_WARNING_POP()

#endif //EPI_UTILS_EIGEN_H
