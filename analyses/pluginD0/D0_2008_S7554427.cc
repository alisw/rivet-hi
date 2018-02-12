// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief D0 Run II Z \f$ p_\perp \f$ differential cross-section shape
  /// @author Andy Buckley
  /// @author Gavin Hesketh
  /// @author Frank Siegert
  class D0_2008_S7554427 : public Analysis {
  public:

    /// Default constructor.
    D0_2008_S7554427()
      : Analysis("D0_2008_S7554427")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      FinalState fs;
      ZFinder zfinder(fs, Cuts::open(), PID::ELECTRON,
                      40*GeV, 200*GeV, 0.2, ZFinder::CLUSTERNODECAY, ZFinder::TRACK);
      declare(zfinder, "ZFinder");

      _h_ZpT         = bookHisto1D(1, 1, 1);
      _h_forward_ZpT = bookHisto1D(3, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& e) {
      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size() != 1) {
        MSG_DEBUG("No unique lepton pair found.");
        vetoEvent;
      }
      const double yZ = fabs(zfinder.bosons()[0].rapidity());
      const double pTZ = zfinder.bosons()[0].pT();
      _h_ZpT->fill(pTZ, e.weight());
      if (yZ > 2) _h_forward_ZpT->fill(pTZ, e.weight());
    }


    // Finalize
    void finalize() {
      normalize(_h_ZpT);
      normalize(_h_forward_ZpT);
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_ZpT, _h_forward_ZpT;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(D0_2008_S7554427);

}
