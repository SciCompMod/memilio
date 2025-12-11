/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Julia Bicker
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
#ifndef SEASON_H
#define SEASON_H

#include "memilio/utils/index.h"

namespace mio
{

/**
 * @brief The Season struct is used as a dynamically
 * sized tag for all season dependent categories
 */
struct Season : public Index<Season> {
    Season(size_t val)
        : Index<Season>(val)
    {
    }

    /**
     * Override deserialize of base class
     * @see mio::Index::deserialize
     */
    template <class IOContext>
    static IOResult<Season> deserialize(IOContext& io)
    {
        BOOST_OUTCOME_TRY(auto&& i, mio::deserialize(io, Tag<size_t>{}));
        return success(Season(i));
    }
};

} // namespace mio

#endif
