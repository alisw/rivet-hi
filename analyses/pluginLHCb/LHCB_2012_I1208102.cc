// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// Differential cross-sections of $\mathrm{Z}/\gamma^* \to e^{+}e^{-}$ vs rapidity and $\phi^*$
  class LHCB_2012_I1208102 : public Analysis {
  public:


    /// Constructor
    LHCB_2012_I1208102()
      : Analysis("LHCB_2012_I1208102")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      ZFinder zeefinder(FinalState(), Cuts::etaIn(2.0, 4.5) && Cuts::pT > 20*GeV, PID::ELECTRON, 60*GeV, 120*GeV);
      declare(zeefinder, "ZeeFinder");

      _h_sigma_vs_y = bookHisto1D(2, 1, 1);
      _h_sigma_vs_phi = bookHisto1D(3, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& e) {
      const ZFinder& zeefinder = apply<ZFinder>(e, "ZeeFinder");
      if (zeefinder.empty()) vetoEvent;
      if (zeefinder.bosons().size() > 1)
        MSG_WARNING("Found multiple (" << zeefinder.bosons().size() << ") Z -> e+ e- decays!");

      // Z momenta
      const FourMomentum& zee = zeefinder.bosons()[0].momentum();
      const Particle& pozitron = zeefinder.constituents()[0];
      const Particle& electron = zeefinder.constituents()[1];

      // Calculation of the angular variable
      const double diffphi = deltaPhi(pozitron, electron);
      const double diffpsd = deltaEta(pozitron, electron);
      const double accphi = M_PI - diffphi;
      const double angular = tan(accphi/2) / cosh(diffpsd/2);

      // Fill histograms
      _h_sigma_vs_y->fill(zee.rapidity(), e.weight());
      _h_sigma_vs_phi->fill(angular, e.weight());
    }


    /// Finalize
    void finalize() {
      const double xs = crossSection()/picobarn;
      scale(_h_sigma_vs_y, xs/sumOfWeights());
      scale(_h_sigma_vs_phi, xs/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_sigma_vs_y, _h_sigma_vs_phi;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(LHCB_2012_I1208102);

}
