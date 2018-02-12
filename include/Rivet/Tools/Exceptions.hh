#ifndef RIVET_EXCEPTIONS_HH
#define RIVET_EXCEPTIONS_HH

#include <string>
#include <exception>
#include <stdexcept>

namespace Rivet {


  /// @brief Generic runtime Rivet error.
  struct Error : public std::runtime_error {
    Error(const std::string& what) : std::runtime_error(what) {}
  };


  /// @brief Rivet::Exception is a synonym for Rivet::Error.
  typedef Error Exception;


  /// @brief Error for e.g. use of invalid bin ranges.
  struct RangeError : public Error {
    RangeError(const std::string& what) : Error(what) {}
  };


  /// @brief Error specialisation for places where alg logic has failed.
  struct LogicError : public Error {
    LogicError(const std::string& what) : Error(what) {}
  };


  /// @brief Error specialisation for failures relating to particle ID codes.
  struct PidError : public Error {
    PidError(const std::string& what) : Error(what) {}
  };


  /// @brief Error specialisation for failures relating to analysis info.
  struct InfoError : public Error {
    InfoError(const std::string& what) : Error(what) {}
  };


  /// @brief Errors relating to event/bin weights
  ///
  /// Arises in computing statistical quantities because e.g. the bin
  /// weight is zero or negative.
  struct WeightError : public Error {
    WeightError(const std::string& what) : Error(what) {}
  };


  /// @brief Error specialisation for where the problem is between the chair and the computer.
  struct UserError : public Error {
    UserError(const std::string& what) : Error(what) {}
  };


}

#endif
