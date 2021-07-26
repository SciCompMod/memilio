#include "epidemiology/utils/io.h"
#include "matchers.h"
#include "gtest/gtest.h"

namespace iotest
{
struct Context {
    void set_error(const epi::IOStatus& err)
    {
        status = err;
    }
    epi::IOStatus status;
};
} // namespace iotest

TEST(IO, apply)
{
    epi::IOResult<int> i{epi::Tag<int>{}, 3};
    epi::IOResult<double> d{epi::Tag<double>{}, 4.0};
    iotest::Context io;
    auto r = epi::apply(
        io,
        [](int i_, double d_) {
            return i_ * d_;
        },
        i, d);
    EXPECT_EQ(r, epi::IOResult<double>(epi::success(12.0)));
}

TEST(IO, apply_with_validation_success)
{
    epi::IOResult<int> i{epi::Tag<int>{}, 3};
    epi::IOResult<double> d{epi::Tag<double>{}, 4.0};
    iotest::Context io;
    auto r = epi::apply(
        io,
        [](int i_, double d_) -> epi::IOResult<double> {
            return epi::success(i_ * d_);
        },
        i, d);
    EXPECT_EQ(r, epi::IOResult<double>(epi::success(12.0)));
}

TEST(IO, apply_with_arg_error)
{
    epi::IOResult<int> i{epi::success(3)};
    epi::IOResult<double> d{epi::failure(epi::StatusCode::OutOfRange, "")};
    iotest::Context io;
    auto r = epi::apply(
        io,
        [](int i_, double d_) {
            return i_ * d_;
        },
        i, d);
    EXPECT_EQ(r, epi::IOResult<double>(epi::failure(epi::StatusCode::OutOfRange, "")));
}

TEST(IO, apply_with_multiple_arg_errors)
{
    epi::IOResult<int> i{epi::failure(epi::StatusCode::InvalidValue, "")};
    epi::IOResult<double> d{epi::failure(epi::StatusCode::OutOfRange, "")};
    iotest::Context io;
    auto r = epi::apply(
        io,
        [](int i_, double d_) {
            return i_ * d_;
        },
        i, d);
    EXPECT_EQ(r, epi::IOResult<double>(epi::failure(epi::StatusCode::InvalidValue, "")));
}

TEST(IO, apply_with_validation_error)
{
    epi::IOResult<int> i{epi::Tag<int>{}, 3};
    epi::IOResult<double> d{epi::Tag<double>{}, 4.0};
    iotest::Context io;
    auto r = epi::apply(
        io,
        [](int, double) -> epi::IOResult<double> {
            return epi::failure(epi::StatusCode::InvalidValue, "");
        },
        i, d);
    EXPECT_EQ(r, epi::IOResult<double>(epi::failure(epi::StatusCode::InvalidValue, "")));
    EXPECT_EQ(io.status, epi::IOStatus(epi::StatusCode::InvalidValue, ""));
}