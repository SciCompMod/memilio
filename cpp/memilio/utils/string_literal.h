/* 
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Rene Schmieding
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
#ifndef MIO_UTILS_STRING_LITERAL_H
#define MIO_UTILS_STRING_LITERAL_H

#include <algorithm>
#include <cstddef>
#include <string_view>

namespace mio
{

/// @brief Wrapper for string literals, that allows passing them as template arguments. Should be used with constexpr.
template <size_t Size>
struct StringLiteral {
    using value_type       = char;
    using pointer          = value_type*;
    using const_pointer    = const value_type*;
    using size_type        = size_t;
    using string_view_type = std::basic_string_view<value_type>;

    /// @brief The length of the StringLiteral (not counting the trailing '\0').
    static constexpr size_type size()
    {
        return Size;
    }

    /**
     * @brief Construct a StringLiteral. Use as `StringLiteral{"example string literal"}`.
     * The type of string is a const reference to a fixed-size char array.
     * @param[in] string_literal Any string literal. 
     */
    constexpr StringLiteral(const value_type (&string_literal)[size() + 1])
    {
        std::copy_n(string_literal, size() + 1, value);
    }

    /// @brief Create a string filled with '\0'. Mainly used by StringLiteral internally.
    constexpr StringLiteral() = default;

    /// @brief Check whether the StringLiteral is empty.
    constexpr bool empty() const
    {
        return size() == 0;
    }

    /// @brief Access the underlying array. Modification is only possible during compile time.
    constexpr pointer data()
    {
        return value;
    }

    /// @brief Access the underlying array.
    constexpr const_pointer data() const
    {
        return value;
    }

    /// @brief Implicit conversion from StringLiteral to string_view
    constexpr operator string_view_type() const
    {
        return string_view_type(data(), size());
    }

    /// @brief Get a string_view of this literal. Be mindful of the lifetime of the view object.
    constexpr string_view_type string_view() const
    {
        return {*this};
    }

    /**
     * @brief Join two literals.
     * @param[in] left, right The strings to join. Operands can be StringLiteral%s or regular string literals.
     * @return A new StringLiteral containing the joined string.
     * @{
     */
    template <size_type N>
    friend constexpr auto operator+(const StringLiteral& left, const StringLiteral<N>& right)
    {
        constexpr size_type new_size = left.size() + right.size();
        StringLiteral<new_size> new_string;
        pointer val = new_string.data();
        val         = std::copy_n(left.data(), left.size(), val); // do not copy the terminating \0
        std::copy_n(right.data(), right.size() + 1, val);
        return new_string;
    }
    template <size_type N>
    friend constexpr auto operator+(const StringLiteral& left, const value_type (&right)[N])
    {
        auto r = StringLiteral<N - 1>(right);
        return left + r;
    }
    template <size_type N>
    friend constexpr auto operator+(const value_type (&left)[N], const StringLiteral& right)
    {
        auto l = StringLiteral<N - 1>(left);
        return l + right;
    }
    /** @} */

    /**
     * @brief Compare two literals.
     * @param[in] left, right The strings to compare. Operands can be StringLiteral%s or regular string literals.
     * @return Whether both strings are the same.
     * @{
     */
    template <size_type N>
    friend constexpr auto operator==(const StringLiteral& left, const StringLiteral<N>& right)
    {
        return left.string_view() == right.string_view();
    }
    template <size_type N>
    friend constexpr auto operator==(const StringLiteral& left, const value_type (&right)[N])
    {
        auto r = StringLiteral<N - 1>(right);
        return left == r;
    }
    /** @} */

    value_type value[size() + 1] = {}; ///< Contains the actual characters. Access this through data() or string_view().
};

// Deduction guide to help the compiler with setting the correct size, i.e. not counting the terminating '\0' char.
template <size_t Size>
StringLiteral(const char (&string)[Size]) -> StringLiteral<Size - 1>;

} // namespace mio

#endif // MIO_UTILS_STRING_LITERAL_H
