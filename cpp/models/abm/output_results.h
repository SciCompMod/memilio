/*
* Copyright (C) 2020-2021 German Aerospace Center (DLR-SC)
*        & Helmholtz Centre for Infection Research (HZI)
*
* Authors: Sascha Korf
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
#ifndef EPI_ABM_OUTPUT_RESULTS_H
#define EPI_ABM_OUTPUT_RESULTS_H

namespace mio
{

namespace abm
{

class OutputResults
{
public:
    OutputResults(bool print_results, bool print_location_results)
        : m_print_results(print_results)
        , m_print_location_results(print_location_results){};

private:
    bool m_print_results;
    bool m_print_location_results;
};
} // namespace abm
} // namespace mio

#endif // EPI_ABM_OUTPUT_RESULTS_H