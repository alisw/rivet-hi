#ifndef RIVET_RivetSTL_HH
#define RIVET_RivetSTL_HH

#include <string>
#include <array>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <utility>
#include <tuple>
#include <algorithm>
#include <type_traits>
#include <stdexcept>
#include <cassert>
#include <memory>
#include <typeinfo>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <functional>


#ifndef foreach
/// @decl A foreach macro for backward compatibility with BOOST_FOREACH
#define foreach(value, container) for (value : container)
#endif


namespace Rivet {


  /// We implicitly use STL entities in the Rivet namespace
  using namespace std;


  /// @name Streaming containers as string reps
  /// @todo Make these named toStr rather than operator<<
  /// @todo Make these generic to any iterable
  //@{

  /// Convenient function for streaming out vectors of any streamable object.
  template<typename T>
  inline std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    os << "[ ";
    for (size_t i=0; i<vec.size(); ++i) {
      os << vec[i] << " ";
    }
    os << "]";
    return os;
  }

  /// Convenient function for streaming out lists of any streamable object.
  template<typename T>
  inline std::ostream& operator<<(std::ostream& os, const std::list<T>& vec) {
    os << "[ ";
    for (size_t i=0; i<vec.size(); ++i) {
      os << vec[i] << " ";
    }
    os << "]";
    return os;
  }

  //@}


  /// @name Boolean-return container searching
  //@{

  /// @todo Use SFINAE, Boost.Range, or other template trickery for more generic container matching?

  /// Does @a s contain @a sub as a substring?
  inline bool contains(const std::string& s, const std::string& sub) {
    return s.find(sub) != string::npos;
  }

  /// Does the vector @a v contain @a x?
  template <typename T>
  inline bool contains(const std::vector<T>& v, const T& x) {
    return find(v.begin(), v.end(), x) != v.end();
  }

  /// Does the list @a l contain @a x?
  template <typename T>
  inline bool contains(const std::list<T>& l, const T& x) {
    return find(l.begin(), l.end(), x) != l.end();
  }

  /// Does the set @a s contain @a x?
  template <typename T>
  inline bool contains(const std::set<T>& s, const T& x) {
    return find(s.begin(), s.end(), x) != s.end();
  }

  /// Does the map @a m contain the key @a key?
  template <typename K, typename T>
  inline bool has_key(const std::map<K, T>& m, const K& key) {
    return m.find(key) != m.end();
  }

  /// Does the map @a m contain the value @a val?
  template <typename K, typename T>
  inline bool has_value(const std::map<K, T>& m, const T& val) {
    for (typename std::map<K,T>::const_iterator it = m.begin(); it != m.end(); ++it) {
      if (it->second == val) return true;
    }
    return false;
  }

  //@}


}

namespace std {


  /// @name Container filling and merging
  //@{

  /// Append a single item to vector @a v
  template <typename T>
  inline void operator+=(std::vector<T>& v, const T& x) { v.push_back(x); }

  /// Append all the items from vector @a v2 to vector @a v1
  template <typename T>
  inline void operator+=(std::vector<T>& v1, const std::vector<T>& v2) {
    for (const auto& x : v2) v1.push_back(x);
  }

  /// Create a new vector from the concatenated items in vectors @a v1 and @a v2
  template <typename T>
  inline std::vector<T> operator+(const std::vector<T>& v1, const std::vector<T>& v2) {
    std::vector<T> rtn(v1);
    rtn += v2;
    return rtn;
  }


  /// Merge the contents of set @a s2 into @a s1
  template <typename T>
  inline void operator+=(std::set<T>& s1, const std::set<T>& s2) {
    for (const auto& x : s2) s1.insert(x);
  }

  /// Merge the contents of sets @a s1 and @a s2
  template <typename T>
  inline std::set<T> operator+(const std::set<T>& s1, const std::set<T>& s2) {
    std::set<T> rtn(s1);
    rtn += s2;
    return rtn;
  }

  //@}


  /// @name Function helpers
  //@{

  /// Get a function pointer / hash integer from an std::function
  template<typename T, typename... U>
  inline size_t get_address(std::function<T(U...)> f) {
    typedef T(fnType)(U...);
    fnType ** fnPointer = f.template target<fnType*>();
    return (fnPointer != nullptr) ? reinterpret_cast<size_t>(*fnPointer) : 0;
  }

  //@}


}

#endif
