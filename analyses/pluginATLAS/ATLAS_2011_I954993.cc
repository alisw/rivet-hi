// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief WZ fiducial cross-section measurement
  class ATLAS_2011_I954993 : public Analysis {
  public:

    /// Default constructor
    ATLAS_2011_I954993()
      : Analysis("ATLAS_2011_I954993")
    {
      setNeedsCrossSection(true);
    }


    /// @name Analysis methods
    //@{

    /// Projection and histogram setup
    void init() {
      FinalState fs;
      Cut cuts = Cuts::abseta < 2.5 && Cuts::pT > 15*GeV;

      ZFinder zfinder_e(fs, cuts, PID::ELECTRON, 81.1876*GeV, 101.1876*GeV, 0.1, ZFinder::CLUSTERNODECAY);
      declare(zfinder_e, "ZFinder_e");
      ZFinder zfinder_mu(fs, cuts, PID::MUON, 81.1876*GeV, 101.1876*GeV, 0.1, ZFinder::CLUSTERNODECAY);
      declare(zfinder_mu, "ZFinder_mu");

      VetoedFinalState weinput;
      weinput.addVetoOnThisFinalState(zfinder_e);
      WFinder wfinder_e(weinput, cuts, PID::ELECTRON, 0*GeV, 1000*GeV, 25*GeV, 0.1, WFinder::CLUSTERNODECAY);
      declare(wfinder_e, "WFinder_e");

      VetoedFinalState wminput;
      wminput.addVetoOnThisFinalState(zfinder_mu);
      WFinder wfinder_mu(wminput,cuts, PID::MUON, 0*GeV, 1000*GeV, 25*GeV, 0.1, WFinder::CLUSTERNODECAY);
      declare(wfinder_mu, "WFinder_mu");

      // Histograms
      _h_fiducial = bookHisto1D(1,1,1);
    }


    /// Do the analysis
    void analyze(const Event& e) {
      const ZFinder& zfinder_e = apply<ZFinder>(e, "ZFinder_e");
      const ZFinder& zfinder_mu = apply<ZFinder>(e, "ZFinder_mu");
      const WFinder& wfinder_e = apply<WFinder>(e, "WFinder_e");
      const WFinder& wfinder_mu = apply<WFinder>(e, "WFinder_mu");

      // Looking for a Z, exit if not found
      if (zfinder_e.bosons().size() != 1 && zfinder_mu.bosons().size() != 1) {
        MSG_DEBUG("No Z boson found, vetoing event");
        vetoEvent;
      }

      // Looking for a W, exit if not found
      if (wfinder_e.bosons().size()!= 1 && wfinder_mu.bosons().size() != 1) {
        MSG_DEBUG("No W boson found, vetoing event");
        vetoEvent;
      }

      // If we find a W, make fiducial acceptance cuts and exit if not found
      if (wfinder_e.bosons().size() == 1) {
        const FourMomentum We = wfinder_e.constituentLeptons()[0];
        const FourMomentum Wenu = wfinder_e.constituentNeutrinos()[0];
        const double mT = wfinder_e.mT();
        if (Wenu.pT() < 25*GeV || We.pT() < 20*GeV || mT < 20*GeV) {
          MSG_DEBUG("Wnu pT = " << Wenu.pT()/GeV << " GeV, Wl pT = " << We.pT()/GeV << " GeV, mT = " << mT/GeV << " GeV");
          vetoEvent;
        }
      } else if (wfinder_mu.bosons().size() == 1) {
        const FourMomentum Wmu = wfinder_mu.constituentLeptons()[0];
        const FourMomentum Wmunu = wfinder_mu.constituentNeutrinos()[0];
        const double mT = wfinder_mu.mT();
        if (Wmunu.pT() < 25*GeV || Wmu.pT() < 20*GeV || mT < 20*GeV) {
          MSG_DEBUG("Wnu pT = " << Wmunu.pT()/GeV << ", Wl pT = " << Wmu.pT()/GeV << " GeV, mT = " << mT/GeV << " GeV");
          vetoEvent;
        }
      } else {
        MSG_DEBUG("No W boson found: vetoing event");
        vetoEvent;
      }

      // Update the fiducial cross-section histogram
      _h_fiducial->fill(7000, e.weight());
    }


    /// Finalize
    void finalize() {
      scale(_h_fiducial, crossSection()/femtobarn/sumOfWeights());
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_fiducial;
    //@}

  };


  //// The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_I954993);

}
