#ifndef EPI_UTILS_UNUSED_H
#define EPI_UTILS_UNUSED_H

namespace epi
{
    /**
     * does nothing, can be used to mark variables as not used.
     * e.g. for avoiding compiler warnings/error about unused variables.
     */
    template<class... T>
    void unused(T&&...) {}

#define QUOTE(x) #x

//defines warning macros for MSVC
#ifdef _MSC_VER
#define MSVC_WARNING_DISABLE_PUSH(list) __pragma(warning(push)) __pragma(warning(disable : list))
#define MSVC_WARNING_POP __pragma(warning(pop))
#else
#define MSVC_WARNING_DISABLE_PUSH(...)
#define MSVC_WARNING_POP
#endif

//defines warning macros for gcc and clang
#if defined(__GNUC__) || defined(__clang__)

#ifdef __GNUC__
#define COMPILER GCC
#else 
#define COMPILER clang
#endif

//there's probably a way to disable multiple warnings in one macro with __VA_ARGS__
#define GCC_CLANG_DIAGNOSTIC_(x) _Pragma(QUOTE(x))
#define GCC_CLANG_DIAGNOSTIC(x) GCC_CLANG_DIAGNOSTIC_(COMPILER diagnostic x)

#else //gcc, clang

#define GCC_CLANG_DIAGNOSTIC(...)

#endif //gcc, clang

//macro for if with compile time constant condition
#define IF_CONSTEXPR(expr)                                                                                             \
    MSVC_WARNING_DISABLE_PUSH(4127)                                                                                    \
    if (expr)                                                                                                          \
    MSVC_WARNING_POP
}

#endif //EPI_UTILS_UNUSED_H
