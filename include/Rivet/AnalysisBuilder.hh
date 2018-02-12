// -*- C++ -*-
#ifndef RIVET_AnalysisBuilder_HH
#define RIVET_AnalysisBuilder_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/Logging.hh"

namespace Rivet {


  // Forward declaration
  class Analysis;


  /// @cond ANALYSIS_PLUGIN_DETAILS

  /// @brief Interface for analysis builders
  class AnalysisBuilderBase {
  public:

    /// Default constructor
    AnalysisBuilderBase() = default;

    /// Constructor with alias name
    AnalysisBuilderBase(const string& alias) : _alias(alias) {}

    /// Destructor
    virtual ~AnalysisBuilderBase() = default;

    /// Factory method, to be implemented by the analysis-specific derived class
    virtual unique_ptr<Analysis> mkAnalysis() const = 0;

    /// Get the analysis' name, by asking it directly
    /// @todo Could avoid this slow lookup by passing it via the constructor... at the cost of potential inconsistency
    string name() const {
      auto a = mkAnalysis();
      return a->name();
    }

    /// @brief Get any optional alias name attached to this builder
    ///
    /// @note An empty string is returned if there is no alias
    const string& alias() const {
      return _alias;
    }

  protected:

    /// The trick: the builder is able to register itself with Rivet's loader
    void _register() {
      AnalysisLoader::_registerBuilder(this);
    }

  private:

    /// Optional alias name
    string _alias;

  };


  /// @brief Self-registering analysis plugin builder
  template <typename T>
  class AnalysisBuilder : public AnalysisBuilderBase {
  public:

    AnalysisBuilder() {
      _register();
    }

    AnalysisBuilder(const string& alias)
      : AnalysisBuilderBase(alias)
    {
      _register();
    }

    unique_ptr<Analysis> mkAnalysis() const {
      return unique_ptr<T>(new T);
    }

  };

  /// @endcond

}

#endif
