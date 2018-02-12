// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"

namespace Rivet {


  /// Generic analysis looking att pp pseudorepidity distributions at
  /// two collision energies. It usess the possibility to read in a
  /// pre-exixsting yoda file, and if histograms for both energies are
  /// filled when finalize() is called, a ratio plot is produced.
  class MC_REENTRANT : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(MC_REENTRANT);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Projections
      const FinalState fs(Cuts::abseta < 5 && Cuts::pT > 500*MeV);
      declare(fs, "FS");
      declare(ChargedFinalState(fs), "CFS");

      // Histograms. Booked for both 900 GeV and 7 TeV and their ratio.
      _histEta70    = bookHisto1D("Eta70", 50, -5, 5);
      _histEta09    = bookHisto1D("Eta09", 50, -5, 5);
      _histEtaR     = bookScatter2D("EtaR", 50, -5, 5);
      fill70 = fill09 = false;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      if (fuzzyEquals(sqrtS()/GeV, 900))
        fill09 = true;
      else if (fuzzyEquals(sqrtS()/GeV, 7000))
        fill70 = true;

      const FinalState& cfs = apply<FinalState>(event, "CFS");
      for (const Particle& p : cfs.particles()) {
        if (fuzzyEquals(sqrtS()/GeV, 900))
          _histEta09->fill(p.eta(), weight);
        else if (fuzzyEquals(sqrtS()/GeV, 7000))
          _histEta70->fill(p.eta(), weight);
      }
    }


    /// Finalize
    void finalize() {
      if ( fill70 ) scale(_histEta70, 1.0/sumOfWeights());
      if ( fill09 ) scale(_histEta09, 1.0/sumOfWeights());
      if ( _histEta70->numEntries() > 0 && _histEta09->numEntries() > 0 )
        divide(_histEta70, _histEta09, _histEtaR);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histEta09, _histEta70;
    Scatter2DPtr _histEtaR;
    //@}

    bool fill09, fill70;
    
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_REENTRANT);

}
