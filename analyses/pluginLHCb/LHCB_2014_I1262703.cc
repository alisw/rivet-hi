// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// @brief Study of forward Z + jet production at 7 TeV at LHCb
  /// @author W. Barter, A. Bursche, M. Sirendi (Rivet implementation)
  class LHCB_2014_I1262703 : public Analysis {
  public:

    /// Default constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(LHCB_2014_I1262703);


    /// Initialise histograms and projections
    void init() {

      // Projections
      const Cut mycut = Cuts::eta >= 2.0 && Cuts::eta <= 4.5 && Cuts::pT > 20*GeV;
      ZFinder zfinder(FinalState(), mycut, PID::MUON, 60*GeV, 120*GeV, 0., ZFinder::NOCLUSTER);
      declare(zfinder, "ZFinder");
      FastJets jetpro(zfinder.remainingFinalState(), FastJets::ANTIKT, 0.5);
      declare(jetpro, "Jets");

      // Histograms
      _h_jet_pT    = bookHisto1D(3, 1, 1);
      _h_jet_eta20 = bookHisto1D(4, 1, 1);
      _h_jet_eta10 = bookHisto1D(4, 1, 2);
      _h_Z_y20     = bookHisto1D(5, 1, 1);
      _h_Z_y10     = bookHisto1D(5, 1, 2);
      _h_Z_pT20    = bookHisto1D(6, 1, 1);
      _h_Z_pT10    = bookHisto1D(6, 1, 2);
      _h_dphi20    = bookHisto1D(7, 1, 1);
      _h_dphi10    = bookHisto1D(7, 1, 2);
      _h_dy20      = bookHisto1D(8, 1, 1);
      _h_dy10      = bookHisto1D(8, 1, 2);
    }


    /// Do the analysis
    void analyze(const Event & e) {

      const ZFinder& zfinder = apply<ZFinder>(e, "ZFinder");
      if (zfinder.bosons().size() != 1) vetoEvent;
      const ParticleVector& leptons = zfinder.constituents();

      const Cut jetSelector = Cuts::eta >= 2.0 && Cuts::eta <= 4.5 && Cuts::pT > 10*GeV;
      const Jets jets = apply<FastJets>(e, "Jets").jetsByPt(jetSelector);

      // Clean the jets against the lepton candidates with a deltaR cut of 0.4
      const Jets cleanedJets = filter_discard(jets, [&](const Jet& j) { return any(leptons, deltaRLess(j, 0.4)); });
      // vector<const Jet*> cleanedJets;
      // for (size_t i = 0; i < jets.size(); i++) {
      //   bool isolated = true;
      //   for (size_t j = 0; j < 2; j++) {
      //     if (deltaR(leptons[j], jets[i]) < 0.4) {
      //       isolated = false;
      //       break;
      //     }
      //   }
      //   if (isolated) cleanedJets.push_back(&jets[i]);
      // }

      // Require at least 1 survivor and note if it is above a 20 GeV jet pT threshold
      if (cleanedJets.empty()) vetoEvent;
      const bool above20 = cleanedJets[0].pT() > 20*GeV;
      const double dphi = deltaPhi(zfinder.boson(), cleanedJets[0]);
      const double drap = zfinder.boson().rap() - cleanedJets[0].rap();

      // Fill histograms
      const double weight = e.weight();
      _h_jet_pT->fill(cleanedJets[0].pT()/GeV, weight);
      _h_jet_eta10->fill(cleanedJets[0].eta(), weight);
      _h_Z_y10->fill(zfinder.boson().rap(), weight);
      _h_Z_pT10->fill(zfinder.boson().pT()/GeV, weight);
      _h_dphi10->fill(dphi, weight);
      _h_dy10->fill(drap, weight);
      if (above20) {
        _h_jet_eta20->fill(cleanedJets[0].eta(), weight);
        _h_Z_y20->fill(zfinder.boson().rap(), weight);
        _h_Z_pT20->fill(zfinder.boson().pT()/GeV, weight);
        _h_dphi20->fill(dphi, weight);
        _h_dy20->fill(drap, weight);
      }

    }


    /// Finalize
    void finalize() {
      normalize({_h_jet_pT, _h_jet_eta20, _h_jet_eta10, _h_Z_y20, _h_Z_y10, _h_Z_pT20, _h_Z_pT10, _h_dphi20, _h_dphi10, _h_dy20, _h_dy10});
    }


    /// Histograms
    Histo1DPtr _h_jet_pT, _h_jet_eta20, _h_jet_eta10, _h_Z_y20, _h_Z_y10, _h_Z_pT20, _h_Z_pT10, _h_dphi20, _h_dphi10, _h_dy20, _h_dy10;

  };


  DECLARE_RIVET_PLUGIN(LHCB_2014_I1262703);

}
