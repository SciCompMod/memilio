/*
* Copyright (C) 2020-2025 MEmilio
*
* Authors: Daniel Abele
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

#include "memilio/utils/metaprogramming.h"
#include <type_traits>

class Foo
{
public:
    void fox(int);
};
class Bar
{
public:
    void baz(int);
};
template <class T>
using fox_member_function_t = decltype(std::declval<T&>().fox(0));
static_assert(mio::is_expression_valid<fox_member_function_t, Foo>::value,
              "Expression Valid: Foo has member function fox");
static_assert(!mio::is_expression_valid<fox_member_function_t, Bar>::value,
              "Expression Not Valid: Bar does not have member function fox");

static_assert(!mio::negation_v<std::is_same<int, int>>, "Negation: Not True = False");
static_assert(mio::negation_v<std::is_same<int, double>>, "Negation: Not False = True");

template <class T>
struct InvalidPredicate {
};
static_assert(mio::conjunction_v<std::is_same<int, int>, std::is_same<double, double>>,
              "Conjunction: True && True = True");
static_assert(!mio::conjunction_v<std::is_same<int, double>, std::is_same<int, int>>,
              "Conjunction: False && True = False");
static_assert(!mio::conjunction_v<std::is_same<int, int>, std::is_same<double, int>>,
              "Conjunction: True && False = False");
static_assert(!mio::conjunction_v<std::is_same<int, double>, std::is_same<double, int>>,
              "Conjunction: False && False = False");
static_assert(!mio::conjunction_v<std::is_same<int, double>, InvalidPredicate<int>>,
              "Conjunction: False && Unevaluated = False (Short Circuit)");

struct NC : mio::not_copyable_if_t<true> {
};
struct C : mio::not_copyable_if_t<false> {
};
static_assert(std::is_copy_constructible<C>::value, "Copyable");
static_assert(!std::is_copy_constructible<NC>::value, "Not copyable");
