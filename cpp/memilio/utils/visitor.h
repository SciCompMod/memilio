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
#ifndef VISITOR_H
#define VISITOR_H

namespace mio
{

/**
 * @brief A generic visitor inspired by Fedor Pikus
 *
 * Fedor Pikus: "C++ Design Patterns: from C++03 to C++17", CppCon 2019
 *
 * Derived class must implement the visit function for all Types
 * that are connected with the visitor.
 *
 * Example usage:
 *
 *   using PetVisitor = Visitor<class Dog, class Cat>;
 *
 *   struct PetMakeNoiseVisitor : public PetVisitor
 *   {
 *       void visit(Cat& c) {
 *           std::cout << "Miau" << std::endl;
 *       }
 *
         void visit(Dog& d) {
 *           std::cout << "Wooof" << std::endl;
 *       }
 *
 *   };
 */

template <typename... Types>
struct Visitor;

template <typename T>
struct Visitor<T> {
    virtual void visit(T& t) = 0;

    virtual ~Visitor() = default;
};

template <typename T, typename... Types>
struct Visitor<T, Types...> : public Visitor<Types...> // Recursive
{
    using Visitor<Types...>::visit;
    virtual void visit(T& t) = 0;

    virtual ~Visitor() = default;
};

template <typename... Types>
struct ConstVisitor;

template <typename T>
struct ConstVisitor<T> {
    virtual void visit(const T& t) = 0;

    virtual ~ConstVisitor() = default;
};

template <typename T, typename... Types>
struct ConstVisitor<T, Types...> : public ConstVisitor<Types...> // Recursive
{
    using ConstVisitor<Types...>::visit;
    virtual void visit(const T& t) = 0;

    virtual ~ConstVisitor() = default;
};

template <class Derived, class Base, class Visitor, class ConstVisitor>
struct Visitable : public Base {
    using Base::Base;

    void accept(Visitor& visitor)
    {
        visitor.visit(*static_cast<Derived*>(this));
    }
    void accept(ConstVisitor& visitor) const
    {
        visitor.visit(*static_cast<Derived const*>(this));
    }
};

} // namespace mio

#endif // VISITOR_H
