// -*- C++ -*-
#include "Rivet/Analyses/MC_JetSplittings.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {




  /// @brief MC validation analysis for higgs [-> tau tau] + jets events
  class MC_HKTSPLITTINGS : public MC_JetSplittings {
  public:

    /// Default constructor
    MC_HKTSPLITTINGS()
      : MC_JetSplittings("MC_HKTSPLITTINGS", 4, "Jets")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      Cut cut = Cuts::abseta < 3.5 && Cuts::pT > 25*GeV;
      /// @todo Urk, abuse! Need explicit HiggsFinder and TauFinder
      ZFinder hfinder(FinalState(), cut, PID::TAU, 115*GeV, 135*GeV, 0.0, ZFinder::NOCLUSTER, ZFinder::NOTRACK, 125*GeV);
      declare(hfinder, "Hfinder");
      FastJets jetpro(hfinder.remainingFinalState(), FastJets::KT, 0.6);
      declare(jetpro, "Jets");

      MC_JetSplittings::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {
      const ZFinder& hfinder = apply<ZFinder>(e, "Hfinder");
      if (hfinder.bosons().size() != 1) vetoEvent;
      MC_JetSplittings::analyze(e);
    }


    /// Finalize
    void finalize() {
      MC_JetSplittings::finalize();
    }

    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_HKTSPLITTINGS);

}
