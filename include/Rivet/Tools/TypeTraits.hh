// -*- C++ -*-
#ifndef RIVET_TypeTraits_HH
#define RIVET_TypeTraits_HH

#include <type_traits>

namespace Rivet {

  /// Mechanisms to allow references and pointers to templated types
  /// to be distinguished from one another (since C++ doesn't allow
  /// partial template specialisation for functions.
  /// Traits methods use specialisation of class/struct templates, and
  /// some trickery with typedefs and static const integral types (or
  /// enums) to implement partial function specialisation as a work-around.

  /// @cond INTERNAL

  namespace SFINAE {
    /// C++11 equivalent of C++17 std::void_t
    template <typename ...>
    using void_t = void;
  }


  struct RefType { };

  struct PtrType { };

  template <typename T>
  struct TypeTraits;

  template <typename U>
  struct TypeTraits<const U&> {
    typedef RefType ArgType;
  };

  template <typename U>
  struct TypeTraits<const U*> {
    typedef PtrType ArgType;
  };



  /// SFINAE definition of dereferenceability trait, cf. Boost has_dereference
  template <typename T, typename=void>
  struct Derefable : std::false_type {};
  //
  template <typename T>
  struct Derefable<T, SFINAE::void_t< decltype(*std::declval<T>())> > : std::true_type {};


  /// SFINAE check for non-const iterability trait
  // template <typename T, typename=void>
  // struct Iterable : std::false_type {};
  // //
  // template <typename T>
  // struct Iterable<T, SFINAE::void_t< decltype(*std::declval<T>())> > : std::true_type {};
  template <typename T>
  using ConstIterable = pretty_print::is_container<T>


  /// @endcond

}

#endif
