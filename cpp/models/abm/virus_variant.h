/* 
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*
* Authors: David Kerkmann
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
#ifndef EPI_VIRUS_VARIANT_H
#define EPI_VIRUS_VARIANT_H

#include <cstdint>

namespace mio
{
namespace abm
{
// Virus Variant handling to be discussed for better solutions
// Ultimately, one would like to read in all Virus Variants from
// a config file.

/**
 * Virus variants in ABM.
 * can be used as 0-based index
*/
enum class VirusVariant : std::uint32_t
{
    Wildtype = 0,

    Count // last!!
};

} // namespace abm
} // namespace mio

#endif
