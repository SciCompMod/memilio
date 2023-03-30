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
#ifndef MIO_IO_SERIALIZER_BASE_H
#define MIO_IO_SERIALIZER_BASE_H

#include "memilio/io/io.h"
#include <memory>

namespace mio
{

/**
 * Base class for implementations of serialization framework concepts.
 * Stores status and flags.
 * See io.h for documentation of the serialization framework.
 */
class SerializerBase
{
public:
    /**
     * Constructor that sets status and flags.
     */
    SerializerBase(std::shared_ptr<IOStatus> status, int flags)
        : m_status(status)
        , m_flags(flags)
    {
        assert(status && "Status must not be null.");
    }

    /**
     * Flags that determine the behavior of serialization.
     * @see mio::IOFlags
     */
    int flags() const
    {
        return m_flags;
    }

    /**
     * Set flags that determine the behavior of serialization.
     * @see mio::IOFlags
     */
    void set_flags(int f)
    {
        m_flags = f;
    }

    /**
     * The current status of serialization.
     * Contains errors that occurred.
     */
    const IOStatus& status() const
    {
        return *m_status;
    }

    /**
     * Set the current status of serialization.
     */
    void set_error(const IOStatus& status)
    {
        if (*m_status) {
            *m_status = status;
        }
    }

protected:
    std::shared_ptr<IOStatus> m_status;
    int m_flags;
};

} // namespace mio

#endif
