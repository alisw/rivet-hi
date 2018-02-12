// -*- C++ -*-
#ifndef RIVET_Utils_HH
#define RIVET_Utils_HH

#include "Rivet/Tools/RivetSTL.hh"
#include "Rivet/Tools/PrettyPrint.hh"
#include <ostream>
#include <iostream>
#include <cctype>
#include <cerrno>
#include <stdexcept>
#include <numeric>
#include <limits>
#include <climits>
#include <cfloat>
#include <cmath>


// // Macro to help with overzealous compiler warnings
// /// @note It's easier and better to just not give an arg name to args which won't be used, when possible.
// #ifdef UNUSED
// #elif defined(__GNUC__)
// # define UNUSED(x) UNUSED_ ## x __attribute__((unused))
// #elif defined(__LCLINT__)
// # define UNUSED(x) /*@unused@*/ x
// #else
// # define UNUSED(x) x
// #endif


/// Macro to help mark code as deprecated to produce compiler warnings
#ifndef DEPRECATED
#if __GNUC__ && __cplusplus && RIVET_NO_DEPRECATION_WARNINGS == 0
#define GCC_VERSION (__GNUC__ * 10000 + __GNUC_MINOR__ * 100 + __GNUC_PATCHLEVEL__)
#if GCC_VERSION >= 40500
  #if __cplusplus > 201103L
  #define DEPRECATED(x) [[deprecated(x)]]
  #else
  #define DEPRECATED(x) __attribute__((deprecated(x)))
  #endif
#else
  #define DEPRECATED(x) __attribute__((deprecated))
#endif
#else
  #define DEPRECATED(x)
#endif
#endif


namespace Rivet {


  /// Convenient const for getting the double NaN value
  static constexpr double DBL_NAN = std::numeric_limits<double>::quiet_NaN();


  /// @name String utils
  //@{

  struct bad_lexical_cast : public std::runtime_error {
    bad_lexical_cast(const std::string& what) : std::runtime_error(what) {}
  };

  /// @brief Convert between any types via stringstream
  template<typename T, typename U>
  T lexical_cast(const U& in) {
    try {
      std::stringstream ss;
      ss << in;
      T out;
      ss >> out;
      return out;
    } catch (const std::exception& e) {
      throw bad_lexical_cast(e.what());
    }
  }

  /// @brief Convert any object to a string
  ///
  /// Just a convenience wrapper for the more general Boost lexical_cast
  template <typename T>
  inline string to_str(const T& x) {
    return lexical_cast<string>(x);
  }

  /// @brief Convert any object to a string
  ///
  /// An alias for to_str() with a more "Rivety" mixedCase name.
  template <typename T>
  inline string toString(const T& x) {
    return to_str(x);
  }

  /// Replace the first instance of patt with repl
  inline string& replace_first(string& str, const string& patt, const string& repl) {
    if (!contains(str, patt)) return str; //< contains from RivetSTL
    str.replace(str.find(patt), patt.size(), repl);
    return str;
  }

  /// @brief Replace all instances of patt with repl
  ///
  /// @note Finding is interleaved with replacement, so the second search happens after
  /// first replacement, etc. This could lead to infinite loops and other counterintuitive
  /// behaviours if not careful.
  inline string& replace_all(string& str, const string& patt, const string& repl) {
    if (!contains(str, patt)) return str; //< contains from RivetSTL
    while (true) {
      string::size_type it = str.find(patt);
      if (it == string::npos) break;
      str.replace(it, patt.size(), repl);
    }
    return str;
  }


  /// Case-insensitive string comparison function
  inline int nocase_cmp(const string& s1, const string& s2) {
    string::const_iterator it1 = s1.begin();
    string::const_iterator it2 = s2.begin();
    while ( (it1 != s1.end()) && (it2 != s2.end()) ) {
      if(::toupper(*it1) != ::toupper(*it2)) { // < Letters differ?
        // Return -1 to indicate smaller than, 1 otherwise
        return (::toupper(*it1) < ::toupper(*it2)) ? -1 : 1;
      }
      // Proceed to the next character in each string
      ++it1;
      ++it2;
    }
    size_t size1 = s1.size(), size2 = s2.size(); // Cache lengths
    // Return -1,0 or 1 according to strings' lengths
    if (size1 == size2) return 0;
    return (size1 < size2) ? -1 : 1;
  }


  /// Case-insensitive string equality function
  inline bool nocase_equals(const string& s1, const string& s2) {
    return nocase_cmp(s1, s2) == 0;
  }


  /// Convert a string to lower-case
  inline string toLower(const string& s) {
    string out = s;
    std::transform(out.begin(), out.end(), out.begin(), (int(*)(int)) std::tolower);
    return out;
  }


  /// Convert a string to upper-case
  inline string toUpper(const string& s) {
    string out = s;
    std::transform(out.begin(), out.end(), out.begin(), (int(*)(int)) std::toupper);
    return out;
  }


  /// Check whether a string @a start is found at the start of @a s
  inline bool startsWith(const string& s, const string& start) {
    if (s.length() < start.length()) return false;
    return s.substr(0, start.length()) == start;
  }


  /// Check whether a string @a end is found at the end of @a s
  inline bool endsWith(const string& s, const string& end) {
    if (s.length() < end.length()) return false;
    return s.substr(s.length() - end.length()) == end;
  }


  /// Make a string containing the string representations of each item in v, separated by sep
  template <typename T>
  inline string join(const vector<T>& v, const string& sep=" ") {
    string rtn;
    for (size_t i = 0; i < v.size(); ++i) {
      if (i != 0) rtn += sep;
      rtn += to_str(v[i]);
    }
    return rtn;
  }

  /// Make a string containing the string representations of each item in s, separated by sep
  template <typename T>
  inline string join(const set<T>& s, const string& sep=" ") {
    string rtn;
    for (const T& x : s) {
      if (rtn.size() > 0) rtn += sep;
      rtn += to_str(x);
    }
    return rtn;
  }

  /// @brief Split a string on a specified separator string
  inline vector<string> split(const string& s, const string& sep) {
    vector<string> dirs;
    string tmp = s;
    while (true) {
      const size_t delim_pos = tmp.find(sep);
      if (delim_pos == string::npos) break;
      const string dir = tmp.substr(0, delim_pos);
      if (dir.length()) dirs.push_back(dir); // Don't insert "empties"
      tmp.replace(0, delim_pos+1, "");
    }
    if (tmp.length()) dirs.push_back(tmp); // Don't forget the trailing component!
    return dirs;
  }

  //@}


  /// @name Path utils
  //@{

  /// @brief Split a path string with colon delimiters
  ///
  /// Ignores zero-length substrings. Designed for getting elements of filesystem paths, naturally.
  inline vector<string> pathsplit(const string& path) {
    return split(path, ":");
  }


  /// @brief Join several filesystem paths together with the standard ':' delimiter
  ///
  /// Note that this does NOT join path elements together with a platform-portable
  /// directory delimiter, cf. the Python @c {os.path.join}!
  inline string pathjoin(const vector<string>& paths) {
    return join(paths, ":");
  }

  /// Operator for joining strings @a a and @a b with filesystem separators
  inline string operator / (const string& a, const string& b) {
    // Ensure that a doesn't end with a slash, and b doesn't start with one, to avoid "//"
    const string anorm = (a.find("/") != string::npos) ? a.substr(0, a.find_last_not_of("/")+1) : a;
    const string bnorm = (b.find("/") != string::npos) ? b.substr(b.find_first_not_of("/")) : b;
    return anorm + "/" + bnorm;
  }

  /// Get the basename (i.e. terminal file name) from a path @a p
  inline string basename(const string& p) {
    if (!contains(p, "/")) return p;
    return p.substr(p.rfind("/")+1);
  }

  /// Get the dirname (i.e. path to the penultimate directory) from a path @a p
  inline string dirname(const string& p) {
    if (!contains(p, "/")) return "";
    return p.substr(0, p.rfind("/"));
  }

  /// Get the stem (i.e. part without a file extension) from a filename @a f
  inline string file_stem(const string& f) {
    if (!contains(f, ".")) return f;
    return f.substr(0, f.rfind("."));
  }

  /// Get the file extension from a filename @a f
  inline string file_extn(const string& f) {
    if (!contains(f, ".")) return "";
    return f.substr(f.rfind(".")+1);
  }

  //@}


  /// @name Container utils
  //@{

  /// Return number of elements in the container @a c for which @c f(x) is true.
  template <typename CONTAINER>
  inline unsigned int count(const CONTAINER& c) {
    // return std::count_if(std::begin(c), std::end(c), [](const typename CONTAINER::value_type& x){return bool(x);});
    unsigned int rtn = 0;
    for (const auto& x : c) if (bool(x)) rtn += 1;
    return rtn;
  }

  /// Return number of elements in the container @a c for which @c f(x) is true.
  template <typename CONTAINER, typename FN>
  inline unsigned int count(const CONTAINER& c, const FN& f) {
    return std::count_if(std::begin(c), std::end(c), f);
  }

  /// Return true if x is true for any x in container c, otherwise false.
  template <typename CONTAINER>
  inline bool any(const CONTAINER& c) {
    // return std::any_of(std::begin(c), std::end(c), [](const auto& x){return bool(x);});
    for (const auto& x : c) if (bool(x)) return true;
    return false;
  }

  /// Return true if f(x) is true for any x in container c, otherwise false.
  template <typename CONTAINER, typename FN>
  inline bool any(const CONTAINER& c, const FN& f) {
    return std::any_of(std::begin(c), std::end(c), f);
  }

  /// Return true if @a x is true for all @c x in container @a c, otherwise false.
  template <typename CONTAINER>
  inline bool all(const CONTAINER& c) {
    // return std::all_of(std::begin(c), std::end(c), [](const auto& x){return bool(x);});
    for (const auto& x : c) if (!bool(x)) return false;
    return true;
  }

  /// Return true if @a f(x) is true for all @c x in container @a c, otherwise false.
  template <typename CONTAINER, typename FN>
  inline bool all(const CONTAINER& c, const FN& f) {
    return std::all_of(std::begin(c), std::end(c), f);
  }

  /// Return true if @a x is false for all @c x in container @a c, otherwise false.
  template <typename CONTAINER>
  inline bool none(const CONTAINER& c) {
    // return std::none_of(std::begin(c), std::end(c), [](){});
    for (const auto& x : c) if (bool(x)) return false;
    return true;
  }

  /// Return true if @a f(x) is false for all @c x in container @a c, otherwise false.
  template <typename CONTAINER, typename FN>
  inline bool none(const CONTAINER& c, const FN& f) {
    return std::none_of(std::begin(c), std::end(c), f);
  }


  /// A single-container-arg version of std::transform, aka @c map
  template <typename C1, typename C2, typename FN>
  inline const C2& transform(const C1& in, C2& out, const FN& f) {
    out.clear(); out.resize(in.size());
    std::transform(in.begin(), in.end(), out.begin(), f);
    return out;
  }

  /// A single-container-arg version of std::accumulate, aka @c reduce
  template <typename C1, typename T, typename FN>
  inline T accumulate(const C1& in, const T& init, const FN& f) {
    const T rtn = std::accumulate(in.begin(), in.end(), init, f);
    return rtn;
  }

  /// Generic sum function, adding @c x for all @c x in container @a c, starting with @a start
  template <typename CONTAINER, typename T>
  inline T sum(const CONTAINER& c, const T& start=T()) {
    T rtn = start;
    for (const auto& x : c) rtn += x;
    return rtn;
  }

  /// Generic sum function, adding @a fn(@c x) for all @c x in container @a c, starting with @a start
  template <typename CONTAINER, typename FN, typename T>
  inline T sum(const CONTAINER& c, const FN& f, const T& start=T()) {
    T rtn = start;
    for (const auto& x : c) rtn += f(x);
    return rtn;
  }


  /// Filter a collection in-place, removing the subset that passes the supplied function
  template <typename CONTAINER, typename FN>
  inline CONTAINER& ifilter_discard(CONTAINER& c, const FN& f) {
    const auto newend = std::remove_if(std::begin(c), std::end(c), f);
    c.erase(newend, c.end());
    return c;
  }

  /// Filter a collection by copy, removing the subset that passes the supplied function
  template <typename CONTAINER, typename FN>
  inline CONTAINER filter_discard(const CONTAINER& c, const FN& f) {
    CONTAINER rtn = c;
    return ifilter_discard(rtn, f); ///< @todo More efficient would be copy_if with back_inserter...
  }

  /// Filter a collection by copy into a supplied container, removing the subset that passes the supplied function
  /// @note New container will be replaced, not appended to
  template <typename CONTAINER, typename FN>
  inline CONTAINER& filter_discard(const CONTAINER& c, const FN& f, CONTAINER& out) {
    out = filter_discard(c, f);
    return out;
  }


  /// Filter a collection in-place, keeping the subset that passes the supplied function
  template <typename CONTAINER, typename FN>
  inline CONTAINER& ifilter_select(CONTAINER& c, const FN& f) {
    //using value_type = typename std::remove_reference<decltype(*std::begin(std::declval<typename std::add_lvalue_reference<CONTAINER>::type>()))>::type;
    auto invf = [&](const typename CONTAINER::value_type& x){ return !f(x); };
    return ifilter_discard(c, invf);
  }

  /// Filter a collection by copy, keeping the subset that passes the supplied function
  template <typename CONTAINER, typename FN>
  inline CONTAINER filter_select(const CONTAINER& c, const FN& f) {
    CONTAINER rtn = c;
    return ifilter_select(rtn, f); ///< @todo More efficient would be copy_if with back_inserter ... but is that equally container agnostic?
  }

  /// Filter a collection by copy into a supplied container, keeping the subset that passes the supplied function
  /// @note New container will be replaced, not appended to
  template <typename CONTAINER, typename FN>
  inline CONTAINER& filter_select(const CONTAINER& c, const FN& f, CONTAINER& out) {
    out = filter_select(c, f);
    return out;
  }


  /// @brief Slice of the container elements cf. Python's [i:j] syntax
  ///
  /// The element at the @j index is not included in the returned container.
  /// @a i and @a j can be negative, treated as backward offsets from the end of the container.
  template <typename CONTAINER>
  inline CONTAINER slice(const CONTAINER& c, int i, int j) {
    CONTAINER rtn;
    const size_t off1 = (i >= 0) ? i : c.size() + i;
    const size_t off2 = (j >= 0) ? j : c.size() + j;
    if (off1 > c.size() || off2 > c.size()) throw RangeError("Attempting to slice beyond requested offsets");
    if (off2 < off1) throw RangeError("Requested offsets in invalid order");
    rtn.resize(off2 - off1);
    std::copy(c.begin()+off1, c.begin()+off2, rtn.begin());
    return rtn;
  }

  /// @brief Tail slice of the container elements cf. Python's [i:] syntax
  ///
  /// Single-index specialisation of @c slice(c, i, j)
  template <typename CONTAINER>
  inline CONTAINER slice(const CONTAINER& c, int i) {
    return slice(c, i, c.size());
  }

  /// @brief Head slice of the @a n first container elements
  ///
  /// Negative @a n means to take the head excluding the @a{n}-element tail
  template <typename CONTAINER>
  inline CONTAINER head(const CONTAINER& c, int n) {
    // if (n > c.size()) throw RangeError("Requested head longer than container");
    if (n < 0) n = std::max(0, (int)c.size()+n);
    n = std::min(n, (int)c.size());
    return slice(c, 0, n);
  }

  /// @brief Tail slice of the @a n last container elements
  ///
  /// Negative @a n means to take the tail from after the @a{n}th element
  template <typename CONTAINER>
  inline CONTAINER tail(const CONTAINER& c, int n) {
    // if (n > c.size()) throw RangeError("Requested tail longer than container");
    if (n < 0) n = std::max(0, (int)c.size()+n);
    n = std::min(n, (int)c.size());
    return slice(c, c.size()-n);
  }


  using std::min;
  using std::max;

  /// Find the minimum value in the vector
  inline double min(const vector<double>& in, double errval=DBL_NAN) {
    const auto e = std::min_element(in.begin(), in.end());
    return e != in.end() ? *e : errval;
  }

  /// Find the maximum value in the vector
  inline double max(const vector<double>& in, double errval=DBL_NAN) {
    const auto e = std::max_element(in.begin(), in.end());
    return e != in.end() ? *e : errval;
  }

  /// Find the minimum and maximum values in the vector
  inline pair<double,double> minmax(const vector<double>& in, double errval=DBL_NAN) {
    const auto e = std::minmax_element(in.begin(), in.end());
    const double rtnmin = e.first != in.end() ? *e.first : errval;
    const double rtnmax = e.second != in.end() ? *e.first : errval;
    return std::make_pair(rtnmin, rtnmax);
  }


  /// Find the minimum value in the vector
  inline int min(const vector<int>& in, int errval=-1) {
    const auto e = std::min_element(in.begin(), in.end());
    return e != in.end() ? *e : errval;
  }

  /// Find the maximum value in the vector
  inline int max(const vector<int>& in, int errval=-1) {
    const auto e = std::max_element(in.begin(), in.end());
    return e != in.end() ? *e : errval;
  }

  /// Find the minimum and maximum values in the vector
  inline pair<int,int> minmax(const vector<int>& in, int errval=-1) {
    const auto e = std::minmax_element(in.begin(), in.end());
    const double rtnmin = e.first != in.end() ? *e.first : errval;
    const double rtnmax = e.second != in.end() ? *e.first : errval;
    return std::make_pair(rtnmin, rtnmax);
  }

  //@}


  /// @brief Get a parameter from a named environment variable, with automatic type conversion
  /// @note Return @a fallback if the variable is not defined, otherwise convert its string to the template type
  /// @todo Should the param name have to be specific to an analysis? Can specialise as an Analysis member fn.
  template <typename T>
  T getEnvParam(const std::string name, const T& fallback) {
    char* env = getenv(name.c_str());
    return env ? lexical_cast<T>(env) : fallback;
  }


}

#endif
