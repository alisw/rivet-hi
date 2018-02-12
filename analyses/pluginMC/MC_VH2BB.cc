// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/ZFinder.hh"
#include "Rivet/Projections/WFinder.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Math/LorentzTrans.hh"

namespace Rivet {


  class MC_VH2BB : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    MC_VH2BB()
      : Analysis("MC_VH2BB")
    {    }

    //@}


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs;
      Cut cut = Cuts::abseta < 3.5 && Cuts::pT > 25*GeV;
      ZFinder zeefinder(fs, cut, PID::ELECTRON, 65*GeV, 115*GeV, 0.2);
      declare(zeefinder, "ZeeFinder");
      ZFinder zmmfinder(fs, cut, PID::MUON, 65*GeV, 115*GeV, 0.2);
      declare(zmmfinder, "ZmmFinder");
      WFinder wefinder(fs, cut, PID::ELECTRON, 60*GeV, 100*GeV, 25*GeV, 0.2);
      declare(wefinder, "WeFinder");
      WFinder wmfinder(fs, cut, PID::MUON, 60*GeV, 100*GeV, 25*GeV, 0.2);
      declare(wmfinder, "WmFinder");

      declare(fs, "FinalState");
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "AntiKT04");
      declare(FastJets(fs, FastJets::ANTIKT, 0.5), "AntiKT05");
      declare(FastJets(fs, FastJets::ANTIKT, 0.6), "AntiKT06");

      /// Book histograms
      _h_jet_bb_Delta_eta = bookHisto1D("jet_bb_Delta_eta", 50, 0, 4);
      _h_jet_bb_Delta_phi = bookHisto1D("jet_bb_Delta_phi", 50, 0, M_PI);
      _h_jet_bb_Delta_pT = bookHisto1D("jet_bb_Delta_pT", 50,0, 500);
      _h_jet_bb_Delta_R = bookHisto1D("jet_bb_Delta_R", 50, 0, 5);
      _h_jet_b_jet_eta = bookHisto1D("jet_b_jet_eta", 50, -4, 4);
      _h_jet_b_jet_multiplicity = bookHisto1D("jet_b_jet_multiplicity", 11, -0.5, 10.5);
      _h_jet_b_jet_phi = bookHisto1D("jet_b_jet_phi", 50, 0, 2.*M_PI);
      _h_jet_b_jet_pT = bookHisto1D("jet_b_jet_pT", 50, 0, 500);
      _h_jet_H_eta_using_bb = bookHisto1D("jet_H_eta_using_bb", 50, -4, 4);
      _h_jet_H_mass_using_bb = bookHisto1D("jet_H_mass_using_bb", 50, 50, 200);
      _h_jet_H_phi_using_bb = bookHisto1D("jet_H_phi_using_bb", 50, 0, 2.*M_PI);
      _h_jet_H_pT_using_bb = bookHisto1D("jet_H_pT_using_bb", 50, 0, 500);
      _h_jet_eta = bookHisto1D("jet_eta", 50, -4, 4);
      _h_jet_multiplicity = bookHisto1D("jet_multiplicity", 11, -0.5, 10.5);
      _h_jet_phi = bookHisto1D("jet_phi", 50, 0, 2.*M_PI);
      _h_jet_pT = bookHisto1D("jet_pT", 50, 0, 500);
      _h_jet_VBbb_Delta_eta = bookHisto1D("jet_VBbb_Delta_eta", 50, 0, 4);
      _h_jet_VBbb_Delta_phi = bookHisto1D("jet_VBbb_Delta_phi", 50, 0, M_PI);
      _h_jet_VBbb_Delta_pT = bookHisto1D("jet_VBbb_Delta_pT", 50, 0, 500);
      _h_jet_VBbb_Delta_R = bookHisto1D("jet_VBbb_Delta_R", 50, 0, 8);

      _h_VB_eta = bookHisto1D("VB_eta", 50, -4, 4);
      _h_VB_mass = bookHisto1D("VB_mass", 50, 60, 110);
      _h_Z_multiplicity = bookHisto1D("Z_multiplicity", 11, -0.5, 10.5);
      _h_W_multiplicity = bookHisto1D("W_multiplicity", 11, -0.5, 10.5);
      _h_VB_phi = bookHisto1D("VB_phi", 50, 0, 2.*M_PI);
      _h_VB_pT = bookHisto1D("VB_pT", 50, 0, 500);

      _h_jet_bVB_angle_Hframe = bookHisto1D("jet_bVB_angle_Hframe", 50, 0, M_PI);
      _h_jet_bVB_cosangle_Hframe = bookHisto1D("jet_bVB_cosangle_Hframe", 50, -1, 1);
      _h_jet_bb_angle_Hframe = bookHisto1D("jet_bb_angle_Hframe", 50, 0, M_PI);
      _h_jet_bb_cosangle_Hframe = bookHisto1D("jet_bb_cosangle_Hframe", 50, -1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const double JETPTCUT = 30*GeV;

      const ZFinder& zeefinder = apply<ZFinder>(event, "ZeeFinder");
      const ZFinder& zmmfinder = apply<ZFinder>(event, "ZmmFinder");
      const WFinder& wefinder = apply<WFinder>(event, "WeFinder");
      const WFinder& wmfinder = apply<WFinder>(event, "WmFinder");
      const Particles vectorBosons = zeefinder.bosons() + zmmfinder.bosons() + wefinder.bosons() + wmfinder.bosons();
      _h_Z_multiplicity->fill(zeefinder.bosons().size() + zmmfinder.bosons().size(), weight);
      _h_W_multiplicity->fill(wefinder.bosons().size() + wmfinder.bosons().size(), weight);

      const Jets jets = apply<FastJets>(event, "AntiKT04").jetsByPt(JETPTCUT);
      _h_jet_multiplicity->fill(jets.size(), weight);

      // Identify the b-jets
      Jets bjets;
      foreach (const Jet& jet, jets) {
        const double jetEta = jet.eta();
        const double jetPhi = jet.phi();
        const double jetPt = jet.pT();
        _h_jet_eta->fill(jetEta, weight);
        _h_jet_phi->fill(jetPhi, weight);
        _h_jet_pT->fill(jetPt/GeV, weight);

        if (jet.bTagged() && jet.pT() > JETPTCUT) {
          bjets.push_back(jet);
          _h_jet_b_jet_eta->fill( jetEta , weight );
          _h_jet_b_jet_phi->fill( jetPhi , weight );
          _h_jet_b_jet_pT->fill( jetPt , weight );
        }
      }
      _h_jet_b_jet_multiplicity->fill(bjets.size(), weight);

      // Plot vector boson properties
      foreach (const Particle& v, vectorBosons) {
        _h_VB_phi->fill(v.phi(), weight);
        _h_VB_pT->fill(v.pT(), weight);
        _h_VB_eta->fill(v.eta(), weight);
        _h_VB_mass->fill(v.mass(), weight);
      }

      // rest of analysis requires at least 1 b jets
      if(bjets.empty()) vetoEvent;

      // Construct Higgs candidates from pairs of b-jets
      for (size_t i = 0; i < bjets.size()-1; ++i) {
        for (size_t j = i+1; j < bjets.size(); ++j) {
          const Jet& jet1 = bjets[i];
          const Jet& jet2 = bjets[j];

          const double deltaEtaJJ = fabs(jet1.eta() - jet2.eta());
          const double deltaPhiJJ = deltaPhi(jet1.momentum(), jet2.momentum());
          const double deltaRJJ = deltaR(jet1.momentum(), jet2.momentum());
          const double deltaPtJJ = fabs(jet1.pT() - jet2.pT());
          _h_jet_bb_Delta_eta->fill(deltaEtaJJ, weight);
          _h_jet_bb_Delta_phi->fill(deltaPhiJJ, weight);
          _h_jet_bb_Delta_pT->fill(deltaPtJJ, weight);
          _h_jet_bb_Delta_R->fill(deltaRJJ, weight);

          const FourMomentum phiggs = jet1.momentum() + jet2.momentum();
          _h_jet_H_eta_using_bb->fill(phiggs.eta(), weight);
          _h_jet_H_mass_using_bb->fill(phiggs.mass(), weight);
          _h_jet_H_phi_using_bb->fill(phiggs.phi(), weight);
          _h_jet_H_pT_using_bb->fill(phiggs.pT(), weight);

          foreach (const Particle& v, vectorBosons) {
            const double deltaEtaVH = fabs(phiggs.eta() - v.eta());
            const double deltaPhiVH = deltaPhi(phiggs, v.momentum());
            const double deltaRVH = deltaR(phiggs, v.momentum());
            const double deltaPtVH = fabs(phiggs.pT() - v.pT());
            _h_jet_VBbb_Delta_eta->fill(deltaEtaVH, weight);
            _h_jet_VBbb_Delta_phi->fill(deltaPhiVH, weight);
            _h_jet_VBbb_Delta_pT->fill(deltaPtVH, weight);
            _h_jet_VBbb_Delta_R->fill(deltaRVH, weight);

            // Calculate boost angles
            const vector<double> angles = boostAngles(jet1.momentum(), jet2.momentum(), v.momentum());
            _h_jet_bVB_angle_Hframe->fill(angles[0], weight);
            _h_jet_bb_angle_Hframe->fill(angles[1], weight);
            _h_jet_bVB_cosangle_Hframe->fill(cos(angles[0]), weight);
            _h_jet_bb_cosangle_Hframe->fill(cos(angles[1]), weight);
          }

        }
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_jet_bb_Delta_eta, crossSection()/sumOfWeights());
      scale(_h_jet_bb_Delta_phi, crossSection()/sumOfWeights());
      scale(_h_jet_bb_Delta_pT, crossSection()/sumOfWeights());
      scale(_h_jet_bb_Delta_R, crossSection()/sumOfWeights());
      scale(_h_jet_b_jet_eta, crossSection()/sumOfWeights());
      scale(_h_jet_b_jet_multiplicity, crossSection()/sumOfWeights());
      scale(_h_jet_b_jet_phi, crossSection()/sumOfWeights());
      scale(_h_jet_b_jet_pT, crossSection()/sumOfWeights());
      scale(_h_jet_H_eta_using_bb, crossSection()/sumOfWeights());
      scale(_h_jet_H_mass_using_bb, crossSection()/sumOfWeights());
      scale(_h_jet_H_phi_using_bb, crossSection()/sumOfWeights());
      scale(_h_jet_H_pT_using_bb, crossSection()/sumOfWeights());
      scale(_h_jet_eta, crossSection()/sumOfWeights());
      scale(_h_jet_multiplicity, crossSection()/sumOfWeights());
      scale(_h_jet_phi, crossSection()/sumOfWeights());
      scale(_h_jet_pT, crossSection()/sumOfWeights());
      scale(_h_jet_VBbb_Delta_eta, crossSection()/sumOfWeights());
      scale(_h_jet_VBbb_Delta_phi, crossSection()/sumOfWeights());
      scale(_h_jet_VBbb_Delta_pT, crossSection()/sumOfWeights());
      scale(_h_jet_VBbb_Delta_R, crossSection()/sumOfWeights());

      scale(_h_VB_eta, crossSection()/sumOfWeights());
      scale(_h_VB_mass, crossSection()/sumOfWeights());
      scale(_h_Z_multiplicity, crossSection()/sumOfWeights());
      scale(_h_W_multiplicity, crossSection()/sumOfWeights());
      scale(_h_VB_phi, crossSection()/sumOfWeights());
      scale(_h_VB_pT, crossSection()/sumOfWeights());

      scale(_h_jet_bVB_angle_Hframe, crossSection()/sumOfWeights());
      scale(_h_jet_bb_angle_Hframe, crossSection()/sumOfWeights());
      scale(_h_jet_bVB_cosangle_Hframe, crossSection()/sumOfWeights());
      scale(_h_jet_bb_cosangle_Hframe, crossSection()/sumOfWeights());
    }


    /// This should take in the four-momenta of two b's (jets/hadrons) and a vector boson, for the process VB*->VBH with H->bb
    /// It should return the smallest angle between the virtual vector boson and one of the b's, in the rest frame of the Higgs boson.
    /// It should also return (as the second element of the vector) the angle between the b's, in the rest frame of the Higgs boson.
    vector<double> boostAngles(const FourMomentum& b1, const FourMomentum& b2, const FourMomentum& vb) {
      const FourMomentum higgsMomentum = b1 + b2;
      const FourMomentum virtualVBMomentum = higgsMomentum + vb;
      const LorentzTransform lt = LorentzTransform::mkFrameTransformFromBeta(higgsMomentum.betaVec());

      const FourMomentum virtualVBMomentumBOOSTED = lt.transform(virtualVBMomentum);
      const FourMomentum b1BOOSTED = lt.transform(b1);
      const FourMomentum b2BOOSTED = lt.transform(b2);

      const double angle1 = b1BOOSTED.angle(virtualVBMomentumBOOSTED);
      const double angle2 = b2BOOSTED.angle(virtualVBMomentumBOOSTED);

      const double anglebb = b1BOOSTED.angle(b2BOOSTED);

      vector<double> rtn;
      rtn.push_back(angle1 < angle2 ? angle1 : angle2);
      rtn.push_back(anglebb);
      return rtn;
    }

    //@}


  private:

    /// @name Histograms
    //@{

    Histo1DPtr _h_Z_multiplicity, _h_W_multiplicity;
    Histo1DPtr _h_jet_bb_Delta_eta, _h_jet_bb_Delta_phi, _h_jet_bb_Delta_pT, _h_jet_bb_Delta_R;
    Histo1DPtr _h_jet_b_jet_eta, _h_jet_b_jet_multiplicity, _h_jet_b_jet_phi, _h_jet_b_jet_pT;
    Histo1DPtr _h_jet_H_eta_using_bb, _h_jet_H_mass_using_bb, _h_jet_H_phi_using_bb, _h_jet_H_pT_using_bb;
    Histo1DPtr _h_jet_eta, _h_jet_multiplicity, _h_jet_phi, _h_jet_pT;
    Histo1DPtr _h_jet_VBbb_Delta_eta, _h_jet_VBbb_Delta_phi, _h_jet_VBbb_Delta_pT, _h_jet_VBbb_Delta_R;
    Histo1DPtr _h_VB_eta, _h_VB_mass, _h_VB_phi, _h_VB_pT;
    Histo1DPtr _h_jet_bVB_angle_Hframe, _h_jet_bb_angle_Hframe, _h_jet_bVB_cosangle_Hframe, _h_jet_bb_cosangle_Hframe;
    //Histo1DPtr _h_jet_cuts_bb_deltaR_v_HpT;

    //@}

  };


  // This global object acts as a hook for the plugin system
  DECLARE_RIVET_PLUGIN(MC_VH2BB);

}
