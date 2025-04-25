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
class NamedTimer : public BasicTimer
{
public:
    static constexpr char Name[] = {Chars..., '\0'};

    /**
     * @brief Get a thread local instance of this timer.
     * 
     * This is the only way to obtain a NamedTimer, as its constructors are disabled. The benefit is, that NamedTimer is
     * inherently threadsafe, the drawback is that 
     *
     */
    static NamedTimer& get_instance()
    {
        // create and return a timer
        // the constructor is only called the first time this function is used (for each thread)
        thread_local static NamedTimer t;
        return t;
    }

private:
    NamedTimer()
    {
        const auto name = std::is_same_v<Tag, void> ? Name : boost::core::demangle(typeid(Tag).name()) + "::" + Name;
        TimerRegistrar::get_instance().add_timer({name, *this, mio::get_omp_thread_id()});
    }

    NamedTimer(NamedTimer&)  = delete;
    NamedTimer(NamedTimer&&) = delete;
};

namespace details
{

/**
 * @brief Get the NamedTimer identified by timer_name.
 * @param timer_name[in] A callable returning a constexpr char, for example `CONST_LITERAL("example")`.
 * @tparam Tag Any type.
 * @tparam Identifier The type of timer_name.
 * @tparam I Index sequence from 0 to strlen(timer_name()).
 * @return A reference to the specified NamedTimer.
 */
template <class Tag, class Identifier, std::size_t... I>
constexpr auto& get_named_timer_impl(Identifier timer_name, std::index_sequence<I...>)
{
    (void)timer_name; // if the name is "", i.e. empty, the Identifier shows up as unused
    return NamedTimer<Tag, timer_name()[I]...>::get_instance();
}

} // namespace details

/**
 * @brief Get the NamedTimer identified by timer_name.
 * @param timer_name[in] A callable returning a constexpr char, for example `CONST_LITERAL("example")`.
 * @tparam Tag Any type, defaults to void.
 * @tparam Identifier The type of timer_name.
 * @return A reference to the specified NamedTimer.
 */
template <class Tag = void, class Identifier>
constexpr auto& get_named_timer(Identifier timer_name)
{
    static_assert(
        std::is_same_v<decltype(std::declval<Identifier>()()), const char*>,
        "Identifier must evaluate to a const char*. Use for example 'get_named_timer(CONST_LITERAL(\"example\"))'");
    // use an index sequence to "explode" the char* givne by timer_name() into a list of single char template values.
    return details::get_named_timer_impl<Tag>(
        timer_name, std::make_index_sequence<::mio::details::constexpr_strlen(timer_name())>{});
}

} // namespace timing

template <class Tag = void>
class AutoTimer
{
public:
    template <class Identifier>
    AutoTimer(Identifier timer_name)
        : m_timer(timing::get_named_timer<Tag>(timer_name))
    {
        m_timer.start();
    }

    AutoTimer(timing::BasicTimer& timer)
        : m_timer(timer)
    {
        m_timer.start();
    }

    ~AutoTimer()
    {
        m_timer.stop();
    }

private:
    timing::BasicTimer& m_timer;
};

} // namespace mio

#endif // MEMILIO_TIMER_TIMERS_H
