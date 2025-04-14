#ifndef MEMILIO_TIMER_BASIC_TIMER_H
#define MEMILIO_TIMER_BASIC_TIMER_H

#include "memilio/timer/definitions.h"

namespace mio
{
namespace timing
{

class BasicTimer
{
public:
    void start()
    {
        m_start_time = get_time_now();
    }

    void stop()
    {
        const TimeType stop_time = get_time_now();
        m_elapsed_time += stop_time - m_start_time;
    }

    DurationType get_elapsed_time() const
    {
        return m_elapsed_time;
    }

private:
    TimeType m_start_time;
    DurationType m_elapsed_time;
};

} // namespace timing
} // namespace mio

#endif // MEMILIO_TIMER_BASIC_TIMER_H
