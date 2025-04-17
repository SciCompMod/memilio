#ifndef MEMILIO_TIMER_REGISTRATION_H
#define MEMILIO_TIMER_REGISTRATION_H

#include "memilio/timer/basic_timer.h"

#include <list>
#include <string>

namespace mio
{
namespace timing
{

struct TimerRegistration {
    std::string name;
    BasicTimer& timer;
    int thread_id;
};

struct Printer {
    virtual void print(const std::list<TimerRegistration>&, std::ostream&) = 0;
};

} // namespace timing
} // namespace mio

#endif // MEMILIO_TIMER_REGISTRATION_H
