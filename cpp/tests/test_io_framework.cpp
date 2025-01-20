/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele, Wadim Koslow
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#include "memilio/io/io.h"
#include "matchers.h"
#include "gtest/gtest.h"

namespace iotest
{
struct Context {
    void set_error(const mio::IOStatus& err)
    {
        status = err;
    }
    mio::IOStatus status;
};
} // namespace iotest

TEST(IO, apply)
{
    mio::IOResult<int> i{mio::Tag<int>{}, 3};
    mio::IOResult<double> d{mio::Tag<double>{}, 4.0};
    iotest::Context io;
    auto r = mio::apply(
        io,
        [](int i_, double d_) {
            return i_ * d_;
        },
        i, d);
    EXPECT_EQ(r, mio::IOResult<double>(mio::success(12.0)));
}

TEST(IO, apply_with_validation_success)
{
    mio::IOResult<int> i{mio::Tag<int>{}, 3};
    mio::IOResult<double> d{mio::Tag<double>{}, 4.0};
    iotest::Context io;
    auto r = mio::apply(
        io,
        [](int i_, double d_) -> mio::IOResult<double> {
            return mio::success(i_ * d_);
        },
        i, d);
    EXPECT_EQ(r, mio::IOResult<double>(mio::success(12.0)));
}

TEST(IO, apply_with_arg_error)
{
    mio::IOResult<int> i{mio::success(3)};
    mio::IOResult<double> d{mio::failure(mio::StatusCode::OutOfRange, "")};
    iotest::Context io;
    auto r = mio::apply(
        io,
        [](int i_, double d_) {
            return i_ * d_;
        },
        i, d);
    EXPECT_EQ(r, mio::IOResult<double>(mio::failure(mio::StatusCode::OutOfRange, "")));
}

TEST(IO, apply_with_multiple_arg_errors)
{
    mio::IOResult<int> i{mio::failure(mio::StatusCode::InvalidValue, "")};
    mio::IOResult<double> d{mio::failure(mio::StatusCode::OutOfRange, "")};
    iotest::Context io;
    auto r = mio::apply(
        io,
        [](int i_, double d_) {
            return i_ * d_;
        },
        i, d);
    EXPECT_EQ(r, mio::IOResult<double>(mio::failure(mio::StatusCode::InvalidValue, "")));
}

TEST(IO, apply_with_validation_error)
{
    mio::IOResult<int> i{mio::Tag<int>{}, 3};
    mio::IOResult<double> d{mio::Tag<double>{}, 4.0};
    iotest::Context io;
    auto r = mio::apply(
        io,
        [](int, double) -> mio::IOResult<double> {
            return mio::failure(mio::StatusCode::InvalidValue, "");
        },
        i, d);
    EXPECT_EQ(r, mio::IOResult<double>(mio::failure(mio::StatusCode::InvalidValue, "")));
    EXPECT_EQ(io.status, mio::IOStatus(mio::StatusCode::InvalidValue, ""));
}
