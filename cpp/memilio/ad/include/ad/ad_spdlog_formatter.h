#ifndef AD_SPDLOG_FORMATTER_H
#define AD_SPDLOG_FORMATTER_H

#include "ad/ad.hpp"
#include "spdlog/fmt/bundled/format.h"
#include "memilio/utils/uncertain_value.h"

template <>
struct fmt::formatter<ad::gt1s<double>::type> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return ctx.end();
    }
    template <typename FormatContext>
    auto format(const ad::gt1s<double>::type& input, FormatContext& ctx) const -> decltype(ctx.out())
    {
        return format_to(ctx.out(), "{}", ad::value(input));
    }
};

template <>
struct fmt::formatter<ad::ga1s<double>::type> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return ctx.end();
    }
    template <typename FormatContext>
    auto format(const ad::ga1s<double>::type& input, FormatContext& ctx) const -> decltype(ctx.out())
    {
        return format_to(ctx.out(), "{}", ad::value(input));
    }
};

template <>
struct fmt::formatter<mio::UncertainValue<ad::gt1s<double>::type>> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return ctx.end();
    }
    template <typename FormatContext>
    auto format(const mio::UncertainValue<ad::gt1s<double>::type>& input, FormatContext& ctx) const
        -> decltype(ctx.out())
    {
        return format_to(ctx.out(), "{}", ad::value(input.value()));
    }
};

template <>
struct fmt::formatter<mio::UncertainValue<ad::ga1s<double>::type>> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return ctx.end();
    }
    template <typename FormatContext>
    auto format(const mio::UncertainValue<ad::ga1s<double>::type>& input, FormatContext& ctx) const
        -> decltype(ctx.out())
    {
        return format_to(ctx.out(), "{}", ad::value(input.value()));
    }
};

#endif // AD_SPDLOG_FORMATTER_H