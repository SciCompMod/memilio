#ifndef MEMILIO_TIMER_TIMERS_H
#define MEMILIO_TIMER_TIMERS_H

#include "memilio/utils/mioomp.h"
#include "memilio/timer/basic_timer.h"
#include "memilio/timer/definitions.h"
#include "memilio/timer/timer_registrar.h"
#include <boost/core/demangle.hpp>
#include <type_traits>

namespace mio
{
namespace timing
{

// second attempt, templates dont look nice at all, but are indep. of the char[] used, only using its content
template <class Tag, char... Chars>
class Timer : public BasicTimer
{
public:
    static constexpr char Name[] = {Chars..., '\0'};

    static Timer& get_instance()
    {
        // create and return a timer
        // the constructor is only called the first time this function is used (for each thread)
        thread_local static Timer t;
        return t;
    }

private:
    Timer()
    {
        const auto name = std::is_same_v<Tag, void> ? Name : boost::core::demangle(typeid(Tag).name()) + "::" + Name;
        TimerRegistrar::get_instance().add_timer({name, *this, mio::get_omp_thread_id()});
    }
};

namespace details
{

template <class Tag, typename Identifier, std::size_t... I>
constexpr auto& get_timer_by_name(Identifier timer_name, std::index_sequence<I...>)
{
    (void)timer_name; // if the name is "", i.e. empty, the Identifier shows up as unused
    return Timer<Tag, timer_name()[I]...>::get_instance();
}

} // namespace details

template <class Tag = void, class Identifier>
constexpr auto& get_timer(Identifier timer_name)
{
    static_assert(std::is_same_v<decltype(std::declval<Identifier>()()), const char*>, "message");
    return details::get_timer_by_name<Tag>(timer_name,
                                           std::make_index_sequence<::mio::details::constexpr_strlen(timer_name())>{});
}

} // namespace timing

template <class Tag = void>
class ScopedTimer
{
public:
    template <class Identifier>
    ScopedTimer(Identifier timer_name)
        : m_timer(timing::get_timer<Tag>(timer_name))
    {
        m_timer.start();
    }

    ScopedTimer(timing::BasicTimer& timer)
        : m_timer(timer)
    {
        m_timer.start();
    }

    ~ScopedTimer()
    {
        m_timer.stop();
    }

private:
    timing::BasicTimer& m_timer;
};

} // namespace mio

#endif // MEMILIO_TIMER_TIMERS_H
