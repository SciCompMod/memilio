#ifndef VISITOR_H
#define VISITOR_H

namespace epi
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

} // namespace epi

#endif // VISITOR_H
