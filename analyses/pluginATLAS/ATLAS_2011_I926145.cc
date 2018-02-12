// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  /// @brief Measurement of electron and muon differential cross section from heavy flavour production
  ///
  /// Lepton cross sections differential in pT
  ///
  /// @author Paul Bell, Holger Schulz
  class ATLAS_2011_I926145 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2011_I926145);


    /// Book histograms and initialise projections before the run
    void init() {

      // Electrons and muons
      Cut cuts = (Cuts::abseta < 1.37 || Cuts::absetaIn(1.52,  2.00)) && Cuts::pT > 7*GeV;
      IdentifiedFinalState elecs(cuts, {PID::ELECTRON, PID::POSITRON});
      declare(elecs, "elecs");
      IdentifiedFinalState muons(Cuts::abseta < 2 && Cuts::pT > 7*GeV, {PID::MUON, PID::ANTIMUON});
      declare(muons, "muons");
      IdentifiedFinalState muons_full(Cuts::abseta < 2.5 && Cuts::pT > 4*GeV, {PID::MUON, PID::ANTIMUON});
      declare(muons_full, "muons_full");

	  Cut cut20 = Cuts::abseta < 2.0;
	  Cut cut25 = Cuts::abseta < 2.5;
      const FinalState fs20(cut20);
      const FinalState fs25(cut25);

      /// @todo Bare Zs ...
      ZFinder zfinder_e(fs20, cut20, PID::ELECTRON, 66.0*GeV, 116.0*GeV, 0.1, ZFinder::NOCLUSTER);
      declare(zfinder_e, "ZFinder_e");
      ZFinder zfinder_mu(fs20, cut20, PID::MUON, 66.0*GeV, 116.0*GeV, 0.1, ZFinder::NOCLUSTER);
      declare(zfinder_mu, "ZFinder_mu");
      ZFinder zfinder_mufull(fs25, cut25, PID::MUON, 66.0*GeV, 116.0*GeV, 0.1, ZFinder::NOCLUSTER);
      declare(zfinder_mufull, "ZFinder_mufull");

      /// @todo ... but dressed Ws?
      WFinder wfinder_e(fs20, cut20, PID::ELECTRON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      declare(wfinder_e, "WFinder_e");
      WFinder wfinder_mu(fs20, cut20, PID::MUON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      declare(wfinder_mu, "WFinder_mu");
      WFinder wfinder_mufull(fs25, cut25, PID::MUON, 60.0*GeV, 100.0*GeV, 25.0*GeV, 0.2);
      declare(wfinder_mufull, "WFinder_mufull");

      // Book histograms
      _histPt_elecs      = bookHisto1D(1 ,1 ,1);
      _histPt_muons      = bookHisto1D(2 ,1 ,1);
      _histPt_muons_full = bookHisto1D(3 ,1 ,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Veto event if no lepton is present
      const FinalState& elecs      = apply<FinalState>(event, "elecs");
      const FinalState& muons      = apply<FinalState>(event, "muons");
      const FinalState& muons_full = apply<FinalState>(event, "muons_full");
      if (elecs.empty() && muons.empty() && muons_full.empty()) vetoEvent;

      // Z veto
      const ZFinder& zfinder_e      = apply<ZFinder>(event, "ZFinder_e");
      const ZFinder& zfinder_mu     = apply<ZFinder>(event, "ZFinder_mu");
      const ZFinder& zfinder_mufull = apply<ZFinder>(event, "ZFinder_mufull");
      if (zfinder_e.bosons().size() > 0 || zfinder_mu.bosons().size() > 0 || zfinder_mufull.bosons().size() > 0) {
        MSG_DEBUG("Num elec Z-bosons found: " << zfinder_e.bosons().size());
        MSG_DEBUG("Num muon Z-bosons found: " << zfinder_mu.bosons().size());
        MSG_DEBUG("Num muon Z-bosons found (|eta|<2.5): " << zfinder_mufull.bosons().size());
        vetoEvent;
      }

      // W veto
      const WFinder& wfinder_e      = apply<WFinder>(event, "WFinder_e");
      const WFinder& wfinder_mu     = apply<WFinder>(event, "WFinder_mu");
      const WFinder& wfinder_mufull = apply<WFinder>(event, "WFinder_mufull");
      if (wfinder_e.bosons().size() > 0 || wfinder_mu.bosons().size() > 0 || wfinder_mufull.bosons().size() > 0) {
        MSG_DEBUG("Num elec W-bosons found: " << wfinder_e.bosons().size());
        MSG_DEBUG("Num muon W-bosons found: " << wfinder_mu.bosons().size());
        MSG_DEBUG("Num muon W-bosons found (|eta|<2.5): " << wfinder_mufull.bosons().size());
        vetoEvent;
      }

      // Electron histogram
      if (elecs.size() > 0) {
        for (const Particle& ele : elecs.particles()) {
          if (ele.pT() < 26.0*GeV) _histPt_elecs->fill(ele.pT()*GeV, weight);
        }
      }

      // Muon histogram
      if (muons.size() > 0) {
        for (const Particle& muo : muons.particles()) {
          if (muo.pT() < 26.0*GeV) _histPt_muons->fill(muo.pT()*GeV, weight);
        }
      }

      // Muon full histogram
      if (muons_full.size() > 0) {
        for (const Particle& muo : muons_full.particles()) {
          if (muo.pT() < 100.0*GeV) _histPt_muons_full->fill(muo.pT()*GeV, weight);
        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_histPt_elecs,      crossSection()/nanobarn/sumOfWeights());
      scale(_histPt_muons,      crossSection()/nanobarn/sumOfWeights());
      scale(_histPt_muons_full, crossSection()/nanobarn/sumOfWeights());
    }


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _histPt_elecs, _histPt_muons, _histPt_muons_full;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2011_I926145);

}
