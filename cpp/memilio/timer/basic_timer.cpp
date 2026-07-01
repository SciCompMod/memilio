#include "memilio/timer/basic_timer.h"
#include "memilio/utils/logging.h"

namespace mio
{
namespace timing
{

void BasicTimer::set_running(bool new_state)
{
#ifndef NDEBUG
    m_is_running = new_state;
#else
    mio::unused(new_state);
#endif
}

void BasicTimer::should_be_running(bool expected, const std::string_view function) const
{
#ifndef NDEBUG
    if (m_is_running != expected) {
        mio::log_error("A BasicTimer was {}running while expected to be {}. "
                       "The offending call was {}. "
                       "Consider using an AutoTimer with name (and scope) to avoid this.",
                       m_is_running ? "" : "not ", expected ? "started" : "stopped", function);
    }
    // else: everything ok.

#else
    mio::unused(expected, function);
#endif
}

} // namespace timing
} // namespace mio
