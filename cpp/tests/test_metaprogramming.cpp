/*
* Copyright (C) 2020-2026 MEmilio
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

struct NC : mio::not_copyable_if_t<true> {
};
struct C : mio::not_copyable_if_t<false> {
};
static_assert(std::is_copy_constructible<C>::value, "Copyable");
static_assert(!std::is_copy_constructible<NC>::value, "Not copyable");
