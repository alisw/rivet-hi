// -*- C++ -*-
#ifndef RIVET_RivetHandler_HH
#define RIVET_RivetHandler_HH

#include "Rivet/Config/RivetCommon.hh"
#include "Rivet/Particle.hh"
#include "Rivet/AnalysisLoader.hh"
#include "Rivet/Tools/RivetYODA.hh"

namespace Rivet {


  // Forward declaration and smart pointer for Analysis
  class Analysis;
  typedef std::shared_ptr<Analysis> AnaHandle;


  // Needed to make smart pointers compare equivalent in the STL set
  struct CmpAnaHandle {
    bool operator() (const AnaHandle& a, const AnaHandle& b) {
      return a.get() < b.get();
    }
  };


  /// A class which handles a number of analysis objects to be applied to
  /// generated events. An {@link Analysis}' AnalysisHandler is also responsible
  /// for handling the final writing-out of histograms.
  class AnalysisHandler {
  public:

    /// @name Constructors and destructors. */
    //@{

    /// Preferred constructor, with optional run name.
    AnalysisHandler(const string& runname="");

    /// @brief Destructor
    /// The destructor is not virtual, as this class should not be inherited from.
    ~AnalysisHandler();

    //@}


  private:

    /// Get a logger object.
    Log& getLog() const;


  public:

    /// @name Run properties
    //@{

    /// Get the name of this run.
    string runName() const { return _runname; }


    /// Get the number of events seen. Should only really be used by external
    /// steering code or analyses in the finalize phase.
    size_t numEvents() const { return _eventcounter.numEntries(); }

    /// @brief Access the sum of the event weights seen
    ///
    /// This is the weighted equivalent of the number of events. It should only
    /// be used by external steering code or analyses in the finalize phase.
    double sumW() const { return _eventcounter.sumW(); }
    /// Access to the sum of squared-weights
    double sumW2() const { return _eventcounter.sumW2(); }

    /// @brief Compatibility alias for sumOfWeights
    ///
    /// @deprecated Prefer sumW
    double sumOfWeights() const { return sumW(); }

    /// @brief Set the sum of weights
    ///
    /// This is useful if Rivet is steered externally and
    /// the analyses are run for a sub-contribution of the events
    /// (but of course have to be normalised to the total sum of weights)
    ///
    /// @todo What about the sumW2 term? That needs to be set coherently. Need a
    /// new version, with all three N,sumW,sumW2 numbers (or a counter)
    /// supplied.
    ///
    /// @deprecated Weight sums are no longer tracked this way...
    void setSumOfWeights(const double& sum) {
      //_sumOfWeights = sum;
      throw Error("Can't set sum of weights independently, since it's now tracked by a Counter. "
                  "Please contact the Rivet authors if you need this.");
    }


    /// Is cross-section information required by at least one child analysis?
    /// @deprecated Should no-longer be an issue: does any generator not write the cross-section?
    bool needCrossSection() const;
    /// Whether the handler knows about a cross-section.
    /// @deprecated Should no-longer be an issue: does any generator not write the cross-section?
    bool hasCrossSection() const;

    /// Get the cross-section known to the handler.
    double crossSection() const { return _xs; }

    /// Set the cross-section for the process being generated.
    /// @todo What about the xsec uncertainty? Add a second, optional argument?
    AnalysisHandler& setCrossSection(double xs);


    /// Set the beam particles for this run
    AnalysisHandler& setRunBeams(const ParticlePair& beams) {
      _beams = beams;
      MSG_DEBUG("Setting run beams = " << beams << " @ " << sqrtS()/GeV << " GeV");
      return *this;
    }

    /// Get the beam particles for this run, usually determined from the first event.
    const ParticlePair& beams() const { return _beams; }

    /// Get beam IDs for this run, usually determined from the first event.
    /// @deprecated Use standalone beamIds(ah.beams()), to clean AH interface
    PdgIdPair beamIds() const;

    /// Get energy for this run, usually determined from the first event.
    /// @deprecated Use standalone sqrtS(ah.beams()), to clean AH interface
    double sqrtS() const;

    /// Setter for _ignoreBeams
    void setIgnoreBeams(bool ignore=true);

    //@}


    /// @name Handle analyses
    //@{

    /// Get a list of the currently registered analyses' names.
    std::vector<std::string> analysisNames() const;

    /// Get the collection of currently registered analyses.
    const std::set<AnaHandle, CmpAnaHandle>& analyses() const {
      return _analyses;
    }

    /// Get a registered analysis by name.
    const AnaHandle analysis(const std::string& analysisname) const;


    /// Add an analysis to the run list by object
    AnalysisHandler& addAnalysis(Analysis* analysis);

    /// @brief Add an analysis to the run list using its name.
    ///
    /// The actual Analysis to be used will be obtained via
    /// AnalysisLoader::getAnalysis(string).  If no matching analysis is found,
    /// no analysis is added (i.e. the null pointer is checked and discarded.
    AnalysisHandler& addAnalysis(const std::string& analysisname);

    /// @brief Add analyses to the run list using their names.
    ///
    /// The actual {@link Analysis}' to be used will be obtained via
    /// AnalysisHandler::addAnalysis(string), which in turn uses
    /// AnalysisLoader::getAnalysis(string). If no matching analysis is found
    /// for a given name, no analysis is added, but also no error is thrown.
    AnalysisHandler& addAnalyses(const std::vector<std::string>& analysisnames);


    /// Remove an analysis from the run list using its name.
    AnalysisHandler& removeAnalysis(const std::string& analysisname);

    /// Remove analyses from the run list using their names.
    AnalysisHandler& removeAnalyses(const std::vector<std::string>& analysisnames);

    //@}


    /// @name Main init/execute/finalise
    //@{

    /// Initialize a run, with the run beams taken from the example event.
    void init(const GenEvent& event);

    /// @brief Analyze the given \a event by reference.
    ///
    /// This function will call the AnalysisBase::analyze() function of all
    /// included analysis objects.
    void analyze(const GenEvent& event);

    /// @brief Analyze the given \a event by pointer.
    ///
    /// This function will call the AnalysisBase::analyze() function of all
    /// included analysis objects, after checking the event pointer validity.
    void analyze(const GenEvent* event);

    /// Finalize a run. This function calls the AnalysisBase::finalize()
    /// functions of all included analysis objects.
    void finalize();

    //@}


    /// @name Histogram / data object access
    //@{

    /// Add a vector of analysis objects to the current state.
    void addData(const std::vector<AnalysisObjectPtr>& aos);

    /// Read analysis plots into the histo collection (via addData) from the named file.
    void readData(const std::string& filename);

    /// Get all analyses' plots as a vector of analysis objects.
    std::vector<AnalysisObjectPtr> getData() const;

    /// Write all analyses' plots (via getData) to the named file.
    void writeData(const std::string& filename) const;

    //@}


  private:

    /// The collection of Analysis objects to be used.
    set<AnaHandle, CmpAnaHandle> _analyses;


    /// @name Run properties
    //@{

    /// Run name
    std::string _runname;

    /// Event counter
    Counter _eventcounter;

    /// Cross-section known to AH
    double _xs, _xserr;

    /// Beams used by this run.
    ParticlePair _beams;

    /// Flag to check if init has been called
    bool _initialised;

    /// Flag whether input event beams should be ignored in compatibility check
    bool _ignoreBeams;

    //@}


  private:

    /// The assignment operator is private and must never be called.
    /// In fact, it should not even be implemented.
    AnalysisHandler& operator=(const AnalysisHandler&);

    /// The copy constructor is private and must never be called.  In
    /// fact, it should not even be implemented.
    AnalysisHandler(const AnalysisHandler&);

  };


}

#endif
