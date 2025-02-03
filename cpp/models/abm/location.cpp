/* 
* Copyright (C) 2020-2025 MEmilio
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
#include "abm/location_type.h"
#include "abm/intervention_type.h"
#include "abm/location.h"
#include "abm/random_events.h"

namespace mio
{
namespace abm
{

Location::Location(LocationType loc_type, LocationId loc_id, size_t num_agegroups, int model_id, uint32_t num_cells)
    : m_type(loc_type)
    , m_id(loc_id)
    , m_parameters(num_agegroups)
    , m_cells(num_cells)
    , m_required_mask(MaskType::None)
    , m_model_id(model_id)
{
    assert(num_cells > 0 && "Number of cells has to be larger than 0.");
}

/*
For every cell in a location we have a transmission factor that is nomalized to m_capacity.volume / m_capacity.persons of 
the location "Home", which is 66. We multiply this rate with the individual size of each cell to obtain a "space per person" factor.
*/
ScalarType Cell::compute_space_per_person_relative() const
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
