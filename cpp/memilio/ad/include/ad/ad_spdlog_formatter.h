#ifndef AD_SPDLOG_FORMATTER_H
#define AD_SPDLOG_FORMATTER_H

#include "ad/ad.hpp"
#include "memilio/utils/logging.h"

template <>
struct fmt::formatter<ad::gt1s<double>::type> {
    constexpr auto parse(format_parse_context& ctx) -> decltype(ctx.begin())
    {
        return ctx.end();
    }
    template <typename FormatContext>
    auto format(const ad::gt1s<double>::type& input, FormatContext& ctx) -> decltype(ctx.out())
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
    auto format(const ad::ga1s<double>::type& input, FormatContext& ctx) -> decltype(ctx.out())
    {
        return format_to(ctx.out(), "{}", ad::value(input));
    }
};

#endif // AD_SPDLOG_FORMATTER_H
