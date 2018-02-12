// -*- C++ -*-
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief MC validation analysis for higgs pairs events (stable Higgses)
  class MC_HHJETS : public MC_JetAnalysis {
  public:

    /// Default constructor
    MC_HHJETS()
      : MC_JetAnalysis("MC_HHJETS", 4, "Jets")
    {    }


    /// @name Analysis methods
    //@{

    /// Book histograms
    void init() {
      IdentifiedFinalState ifs(Cuts::abseta < 10.0 && Cuts::pT > 0*GeV);
      ifs.acceptId(25);
      declare(ifs,"IFS");

      VetoedFinalState vfs;
      vfs.addVetoPairId(25);
      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "Jets");

      _h_HH_mass = bookHisto1D("HH_mass", 250, 240, 4000.0);
      _h_HH_dR = bookHisto1D("HH_dR", 25, 0.5, 10.0);
      _h_HH_dPhi = bookHisto1D("HH_dPhi", 64, 0, 3.2);
      _h_HH_deta= bookHisto1D("HH_deta", 50, -5, 5);
      _h_H_pT = bookHisto1D("H_pT", 50, 0, 2000.0);
      _h_HH_pT = bookHisto1D("HH_pT", 200, 0, 2000.0);
      _h_H_pT1 = bookHisto1D("H_pT1", 200, 0, 2000.0);
      _h_H_pT2 = bookHisto1D("H_pT2", 200, 0, 2000.0);
      _h_H_eta = bookHisto1D("H_eta", 50, -5.0, 5.0);
      _h_H_eta1 = bookHisto1D("H_eta1", 50, -5.0, 5.0);
      _h_H_eta2 = bookHisto1D("H_eta2", 50, -5.0, 5.0);
      _h_H_phi = bookHisto1D("H_phi", 25, 0.0, TWOPI);
      _h_H_jet1_deta = bookHisto1D("H_jet1_deta", 50, -5.0, 5.0);
      _h_H_jet1_dR = bookHisto1D("H_jet1_dR", 25, 0.5, 7.0);

      MC_JetAnalysis::init();
    }



    /// Do the analysis
    void analyze(const Event & e) {

      const IdentifiedFinalState& ifs = apply<IdentifiedFinalState>(e, "IFS");
      Particles allp = ifs.particlesByPt();
      if (allp.empty()) vetoEvent;

      const double weight = e.weight();

      FourMomentum hmom = allp[0].momentum();
      if (allp.size() > 1) {
        FourMomentum hmom2(allp[1].momentum());
        _h_HH_dR->fill(deltaR(hmom, hmom2), weight);
        _h_HH_dPhi->fill(deltaPhi(hmom, hmom2), weight);
        _h_HH_deta->fill(hmom.eta()-hmom2.eta(), weight);
        _h_HH_pT->fill((hmom+hmom2).pT(), weight);
        _h_HH_mass->fill((hmom+hmom2).mass(), weight);

        if (hmom.pT() > hmom2.pT()) {
          _h_H_pT1->fill(hmom.pT(), weight);
          _h_H_eta1->fill(hmom.eta(), weight);
          _h_H_pT2->fill(hmom2.pT(), weight);
          _h_H_eta2->fill(hmom2.eta(), weight);
        } else {
          _h_H_pT1->fill(hmom2.pT(), weight);
          _h_H_eta1->fill(hmom2.eta(), weight);
          _h_H_pT2->fill(hmom.pT(), weight);
          _h_H_eta2->fill(hmom.eta(), weight);
        }
      }
      _h_H_pT->fill(hmom.pT(), weight);
      _h_H_eta->fill(hmom.eta(), weight);
      _h_H_phi->fill(hmom.azimuthalAngle(), weight);


      // Get the jet candidates
      Jets jets = apply<FastJets>(e, "Jets").jetsByPt(20.0*GeV);
      if (!jets.empty()) {
        _h_H_jet1_deta->fill(deltaEta(hmom, jets[0]), weight);
        _h_H_jet1_dR->fill(deltaR(hmom, jets[0]), weight);
      }

      MC_JetAnalysis::analyze(e);
    }


    /// Finalize
    void finalize() {
      scale(_h_HH_mass, crossSection()/sumOfWeights());
      scale(_h_HH_dR, crossSection()/sumOfWeights());
      scale(_h_HH_deta, crossSection()/sumOfWeights());
      scale(_h_HH_dPhi, crossSection()/sumOfWeights());
      scale(_h_H_pT, crossSection()/sumOfWeights());
      scale(_h_H_pT1, crossSection()/sumOfWeights());
      scale(_h_H_pT2, crossSection()/sumOfWeights());
      scale(_h_HH_pT, crossSection()/sumOfWeights());
      scale(_h_H_eta, crossSection()/sumOfWeights());
      scale(_h_H_eta1, crossSection()/sumOfWeights());
      scale(_h_H_eta2, crossSection()/sumOfWeights());
      scale(_h_H_phi, crossSection()/sumOfWeights());
      scale(_h_H_jet1_deta, crossSection()/sumOfWeights());
      scale(_h_H_jet1_dR, crossSection()/sumOfWeights());

      MC_JetAnalysis::finalize();
    }

    //@}


  private:

    /// @name Histograms
    //@{
    Histo1DPtr _h_HH_mass, _h_HH_pT, _h_HH_dR, _h_HH_deta, _h_HH_dPhi;
    Histo1DPtr _h_H_pT, _h_H_pT1, _h_H_pT2, _h_H_eta, _h_H_eta1, _h_H_eta2, _h_H_phi;
    Histo1DPtr _h_H_jet1_deta, _h_H_jet1_dR;
    //@}

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_HHJETS);

}
