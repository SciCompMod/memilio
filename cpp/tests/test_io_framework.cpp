/*
* Copyright (C) 2020-2026 MEmilio
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
#include "memilio/utils/logging.h"
#include "temp_file_register.h"
#include "utils.h"
#include "gmock/gmock.h"
#include "gtest/gtest.h"
#include <filesystem>

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

TEST(IO, create_directory)
{
    mio::RedirectLogger logger(mio::LogLevel::info);
    mio::IOResult<bool> result = mio::success(true);
    TempFileRegister tmp;
    std::string model          = "%%%%-%%%%-%%%%-%%%%";
    const auto dir             = std::filesystem::path(tmp.get_unique_path(model));
    const auto dangling_subdir = std::filesystem::path(tmp.get_unique_path(model + "/" + model));

    logger.capture();
    // case group: do not create parents
    // case: create new directory; expect success, value true
    result = mio::create_directory(dir);
    ASSERT_THAT(print_wrap(result), IsSuccess());
    EXPECT_TRUE(result.value());
    EXPECT_THAT(logger.read(), ::testing::HasSubstr("was created"));
    // case: recreate existing directory; expect success, value false
    result = mio::create_directory(dir);
    ASSERT_THAT(print_wrap(result), IsSuccess());
    EXPECT_FALSE(result.value());
    EXPECT_THAT(logger.read(), ::testing::HasSubstr("already exists"));
    // case: create a directory in a non-existent subdirectory; expect failure
    result = mio::create_directory(dangling_subdir);
    ASSERT_TRUE(result.has_failure()); // error code may be OS dependant, so we don't use IsFailure
    EXPECT_TRUE(logger.read().empty());
    EXPECT_THAT(result.error().message(), ::testing::HasSubstr("Failed to create"));

    // case group: create parents
    // case: create new directory; expect success, value true
    result = mio::create_directory(dangling_subdir, true);
    ASSERT_THAT(print_wrap(result), IsSuccess());
    EXPECT_TRUE(result.value());
    EXPECT_THAT(logger.read(), ::testing::HasSubstr("was created"));
    // case: recreate existing directory; expect success, value false
    result = mio::create_directory(dangling_subdir, true);
    ASSERT_THAT(print_wrap(result), IsSuccess());
    EXPECT_FALSE(result.value());
    EXPECT_THAT(logger.read(), ::testing::HasSubstr("already exists"));
    // omit case for failure, as that would require file system problems (e.g. insufficient space or permissions)

    logger.release();
}

TEST(IO, create_directories_or_exit)
{
    mio::set_death_test_mode();
    TempFileRegister tmp;
    std::string model          = "%%%%-%%%%-%%%%-%%%%";
    const auto dangling_subdir = std::filesystem::path(tmp.get_unique_path(model + "/" + model));

    // test cases here do not cover all create_directory cases
    // case: create new dangling directory without creating parents; expect exit
    EXPECT_DEATH(
        {
            mio::LogLevelOverride llo(mio::LogLevel::off);
            mio::create_directories_or_exit(dangling_subdir, false);
        },
        ::testing::IsEmpty());
    // case: create new dangling directory, now with creating parents: expect directory to be created
    std::filesystem::path result_path = mio::create_directories_or_exit(dangling_subdir, true);
    EXPECT_TRUE(result_path.is_absolute());
    EXPECT_TRUE(std::filesystem::exists(result_path));
}
