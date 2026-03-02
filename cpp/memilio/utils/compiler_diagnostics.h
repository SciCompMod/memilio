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
#ifndef EPI_UTILS_UNUSED_H
#define EPI_UTILS_UNUSED_H

namespace mio
{
/**
* does nothing, can be used to mark variables as not used.
* e.g. for avoiding compiler warnings/error about unused variables.
*/
template <class... T>
void unused(T&&...)
{
}

} // namespace mio

#define QUOTE(x) #x

//defines warning macros for MSVC
#ifdef _MSC_VER
#define MSVC_WARNING_DISABLE_PUSH(list) __pragma(warning(push)) __pragma(warning(disable : list))
#define MSVC_WARNING_POP() __pragma(warning(pop))
#else
#define MSVC_WARNING_DISABLE_PUSH(...)
#define MSVC_WARNING_POP()
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

//MSVC has a long standing bug that breaks empty base optimization (EBO)
//if the class has more than one empty base class.
//types that rely on empty base optimization must add this declaration
//e.g. struct MEMILIO_ENABLE_EBO Foo : EmptyBase1, EmptyBase2, ...
//see https://en.cppreference.com/w/cpp/language/ebo
//see https://stackoverflow.com/questions/12701469/why-is-the-empty-base-class-optimization-ebo-is-not-working-in-msvc
#ifdef _MSC_VER
#define MEMILIO_ENABLE_EBO __declspec(empty_bases)
#else
#define MEMILIO_ENABLE_EBO
#endif

#endif //EPI_UTILS_UNUSED_H
