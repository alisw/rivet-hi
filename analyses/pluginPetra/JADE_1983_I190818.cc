// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class JADE_1983_I190818 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(JADE_1983_I190818);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const ChargedFinalState cfs;
      declare(cfs, "CFS");
      if( !(fuzzyEquals(sqrtS()/GeV,12.0) ||
	    fuzzyEquals(sqrtS()/GeV,30.0) ||
	    fuzzyEquals(sqrtS()/GeV,35.0) )) {
        MSG_WARNING("CoM energy of events sqrt(s) = " << sqrtS()/GeV
                    << " doesn't match any available analysis energy .");
      }
      _hist = bookProfile1D(1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const FinalState& cfs = apply<FinalState>(event, "CFS");
      MSG_DEBUG("Total charged multiplicity = " << cfs.size());
      _hist->fill(sqrtS(),cfs.size(),event.weight());
    }


    /// Normalise histograms etc., after the run
    void finalize() {}
    //@}

  private:

    // Histogram
    Profile1DPtr _hist;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(JADE_1983_I190818);

}
