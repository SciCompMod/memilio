/* 
* Copyright (C) 2020-2023 German Aerospace Center (DLR-SC)
*
* Authors: Lena Ploetzke
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
#ifndef LCTSECIR_PARAMETERS_IO_H
#define LCTSECIR_PARAMETERS_IO_H

#include "memilio/config.h"

#ifdef MEMILIO_HAS_JSONCPP
#include "lct_secir/parameters.h"
#include "lct_secir/infection_state.h"

#include "memilio/io/epi_data.h"
#include "memilio/io/io.h"
#include "memilio/utils/date.h"
#include "memilio/math/eigen.h"
#include "memilio/utils/logging.h"
#include <string>

namespace mio
{
namespace lsecir
{

IOResult<Eigen::VectorXd> get_initial_data_from_file(std::string const& path, Date date, InfectionState infectionState,
                                                     Parameters&& parameters, ScalarType scale_confirmed_cases = 1.)
{
    BOOST_OUTCOME_TRY(rki_data, mio::read_confirmed_cases_germany(path));
    auto max_date_entry = std::max_element(rki_data.begin(), rki_data.end(), [](auto&& a, auto&& b) {
        return a.date < b.date;
    });
    if (max_date_entry == rki_data.end()) {
        log_error("RKI data file is empty.");
        return failure(StatusCode::InvalidFileFormat, path + ", file is empty.");
    }
    auto max_date = max_date_entry->date;
    if (max_date < date) {
        log_error("Specified date does not exist in RKI data");
        return failure(StatusCode::OutOfRange, path + ", specified date does not exist in RKI data.");
    }
    // was soll die 6 da?
    auto days_surplus = std::min(get_offset_in_days(max_date, date) - 6, 0);

    Eigen::VectorXd init(infectionState.get_count());
    return init;
}

} // namespace lsecir
} // namespace mio
#endif // MEMILIO_HAS_JSONCPP

#endif // LCTSECIR_PARAMETERS_IO_H