// Taken from: https://github.com/google/googletest/pull/3086/files
// Sanitizer Integration
// The
// [Undefined Behavior Sanitizer](https://clang.llvm.org/docs/UndefinedBehaviorSanitizer.html),
// [Address Sanitizer](https://github.com/google/sanitizers/wiki/AddressSanitizer),
// and
// [Thread Sanitizer](https://github.com/google/sanitizers/wiki/ThreadSanitizerCppManual)
// all provide weak functions that you can override to trigger explicit failures
// when they detect sanitizer errors, such as creating a reference from `nullptr`.
// To override these functions, place definitions for them in a source file that
// you compile as part of your main binary
#include <gtest/gtest.h>

extern "C" {
void __ubsan_on_report()
{
    FAIL() << "Encountered an undefined behavior sanitizer error";
}
void __asan_on_error()
{
    FAIL() << "Encountered an address sanitizer error";
}
void __tsan_on_report()
{
    FAIL() << "Encountered a thread sanitizer error";
}
} // extern "C"
