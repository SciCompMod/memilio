/* 
* Copyright (C) 2020-2024 MEmilio
*
* Authors: Daniel Abele, Elisabeth Kluth, Carlotta Gerstein, Martin J. Kuehn, Khoa Nguyen, David Kerkmann
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
#include "abm/mask_type.h"
#include "abm/location.h"
#include "abm/parameters.h"
#include "abm/random_events.h"
#include "abm/infection.h"

namespace mio
{
namespace abm
{

Location::Location(LocationId loc_id, size_t num_agegroups, uint32_t num_cells)
    : m_id(loc_id)
    , m_parameters(num_agegroups)
    , m_cells(num_cells)
    , m_required_mask(MaskType::Community)
    , m_npi_active(false)
{
    assert(num_cells > 0 && "Number of cells has to be larger than 0.");
}

Location Location::copy_location_without_persons(size_t) const
{
    return *this;
}

/*
For every cell in a location we have a transmission factor that is nomalized to m_capacity.volume / m_capacity.persons of 
the location "Home", which is 66. We multiply this rate with the individual size of each cell to obtain a "space per person" factor.
*/
ScalarType Cell::compute_space_per_person_relative()
{
    if (m_capacity.volume != 0) {
        return 66.0 / m_capacity.volume;
    }
    else {
        return 1.0;
    }
}

} // namespace abm
} // namespace mio
