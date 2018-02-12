// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief MC validation analysis for Z[ee]Z[mumu] + jets events
  class MC_ZZJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_ZZJETS()
      : MC_JetAnalysis("MC_ZZJETS", 4, "Jets")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      Cut cut = Cuts::abseta < 3.5 && Cuts::pT > 25*GeV;
      ZFinder zeefinder(FinalState(), cut, PID::ELECTRON, 65*GeV, 115*GeV, 0.2, ZFinder::CLUSTERNODECAY, ZFinder::TRACK);
      declare(zeefinder, "ZeeFinder");

      VetoedFinalState zmminput;
      zmminput.addVetoOnThisFinalState(zeefinder);
      ZFinder zmmfinder(zmminput, cut, PID::MUON, 65*GeV, 115*GeV, 0.2, ZFinder::CLUSTERNODECAY, ZFinder::TRACK);
      declare(zmmfinder, "ZmmFinder");

      VetoedFinalState jetinput;
      jetinput
          .addVetoOnThisFinalState(zeefinder)
          .addVetoOnThisFinalState(zmmfinder);
      FastJets jetpro(jetinput, FastJets::ANTIKT, 0.4);
      declare(jetpro, "Jets");

      // Correlations with jets
      _h_ZZ_jet1_deta = bookHisto1D("ZZ_jet1_deta", 70, -7.0, 7.0);
      _h_ZZ_jet1_dR = bookHisto1D("ZZ_jet1_dR", 25, 1.5, 7.0);
      _h_Ze_jet1_dR = bookHisto1D("Ze_jet1_dR", 25, 0.0, 7.0);

      // Global stuff
      _h_HT = bookHisto1D("HT", logspace(100, 100.0, 0.5*(sqrtS()>0.?sqrtS():14000.)/GeV));

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event& e) {
      const double weight = e.weight();

      const ZFinder& zeefinder = apply<ZFinder>(e, "ZeeFinder");
      if (zeefinder.bosons().size() != 1) vetoEvent;

      const ZFinder& zmmfinder = apply<ZFinder>(e, "ZmmFinder");
      if (zmmfinder.bosons().size() != 1) vetoEvent;

      // Z momenta
      const FourMomentum& zee = zeefinder.bosons()[0].momentum();
      const FourMomentum& zmm = zmmfinder.bosons()[0].momentum();
      const FourMomentum zz = zee + zmm;
      // Lepton momenta
      const FourMomentum& ep = zeefinder.constituents()[0].momentum();
      const FourMomentum& em = zeefinder.constituents()[1].momentum();
      const FourMomentum& mp = zmmfinder.constituents()[0].momentum();
      const FourMomentum& mm = zmmfinder.constituents()[1].momentum();

      const Jets& jets = apply<FastJets>(e, "Jets").jetsByPt(_jetptcut);
      if (jets.size() > 0) {
        const FourMomentum j0 = jets[0].momentum();
        _h_ZZ_jet1_deta->fill(zz.eta()-j0.eta(), weight);
        _h_ZZ_jet1_dR->fill(deltaR(zz, j0), weight);
        _h_Ze_jet1_dR->fill(deltaR(ep, j0), weight);
      }

      double HT = ep.pT() + em.pT() + mp.pT() + mm.pT();
      foreach (const Jet& jet, jets) HT += jet.pT();
      if (HT > 0.0) _h_HT->fill(HT/GeV, weight);

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      const double s = crossSection()/picobarn/sumOfWeights();
      scale(_h_ZZ_jet1_deta, s);
      scale(_h_ZZ_jet1_dR, s);
      scale(_h_Ze_jet1_dR, s);
      scale(_h_HT, s);
      MC_JetAnalysis::finalize();
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_ZZ_jet1_deta;
    Histo1DPtr _h_ZZ_jet1_dR;
    Histo1DPtr _h_Ze_jet1_dR;
    Histo1DPtr _h_HT;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_ZZJETS);

}
