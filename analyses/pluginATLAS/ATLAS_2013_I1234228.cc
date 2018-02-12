// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class ATLAS_2013_I1234228 : public Analysis {
  public:


    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2013_I1234228);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      const FinalState fs;
      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT > 25*GeV;
      ZFinder zfinder(fs, cuts, PID::ELECTRON, 116*GeV, 1500*GeV, 0.1);      
      declare(zfinder, "ZFinder");
      
      _hist_mll = bookHisto1D(1, 1, 2);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

    	const double weight = event.weight();
      const ZFinder& zfinder = apply<ZFinder>(event, "ZFinder");

    	if (zfinder.bosons().size() != 1)  vetoEvent;
	
      double mass = zfinder.bosons()[0].mass();
	    _hist_mll->fill(mass, weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double sf = crossSection()/sumOfWeights();
      scale(_hist_mll, sf);
    }

    //@}

  private:


    /// @name Histograms
    //@{
    Histo1DPtr _hist_mll;
    //@}

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1234228);

}
