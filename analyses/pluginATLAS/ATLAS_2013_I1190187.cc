// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// ATLAS Wee Wemu Wmumu analysis at Z TeV
  class ATLAS_2013_I1190187 : public Analysis {
  public:

    /// Default constructor
    ATLAS_2013_I1190187()
      : Analysis("ATLAS_2013_I1190187")
    {    }


    void init() {
      FinalState fs;

      Cut etaRanges_EL = (Cuts::abseta < 1.37 || Cuts::absetaIn(1.52, 2.47)) && Cuts::pT > 20*GeV;
      Cut etaRanges_MU = Cuts::abseta < 2.4 && Cuts::pT > 20*GeV;

      MissingMomentum met(fs);
      declare(met, "MET");

      IdentifiedFinalState Photon(fs);
      Photon.acceptIdPair(PID::PHOTON);

      IdentifiedFinalState bare_EL(fs);
      bare_EL.acceptIdPair(PID::ELECTRON);

      IdentifiedFinalState bare_MU(fs);
      bare_MU.acceptIdPair(PID::MUON);

      IdentifiedFinalState neutrinoFS(fs);
      neutrinoFS.acceptNeutrinos();
      declare(neutrinoFS, "Neutrinos");

      ////////////////////////////////////////////////////////
      // DRESSED LEPTONS
      //    3.arg: 0.1      = dR(lep,phot)
      //    4.arg: true     = do clustering
      //    7.arg: false    = ignore photons from hadron or tau
      //
      //////////////////////////////////////////////////////////
      DressedLeptons electronFS(Photon, bare_EL, 0.1, etaRanges_EL);
      declare(electronFS, "ELECTRON_FS");

      DressedLeptons muonFS(Photon, bare_MU, 0.1, etaRanges_MU);
      declare(muonFS, "MUON_FS");

      VetoedFinalState jetinput;
      jetinput.addVetoOnThisFinalState(bare_MU);
      jetinput.addVetoOnThisFinalState(neutrinoFS);

      FastJets jetpro(jetinput, FastJets::ANTIKT, 0.4);
      declare(jetpro, "jet");

      // Book histograms
      _h_Wl1_pT_mumu = bookHisto1D(1, 1, 2);
      _h_Wl1_pT_ee = bookHisto1D(1, 1, 1);
      _h_Wl1_pT_emu = bookHisto1D(1, 1, 3);
      _h_Wl1_pT_inclusive = bookHisto1D(4, 1, 1);
    }


    /// Do the analysis
    void analyze(const Event& e) {

      const  vector<DressedLepton>& muonFS = apply<DressedLeptons>(e, "MUON_FS").dressedLeptons();
      const  vector<DressedLepton>& electronFS = apply<DressedLeptons>(e, "ELECTRON_FS").dressedLeptons();
      const MissingMomentum& met = apply<MissingMomentum>(e, "MET");

      vector<DressedLepton> dressed_lepton, isolated_lepton, fiducial_lepton;
      dressed_lepton.insert(dressed_lepton.end(), muonFS.begin(), muonFS.end());
      dressed_lepton.insert(dressed_lepton.end(), electronFS.begin(), electronFS.end());

      ////////////////////////////////////////////////////////////////////////////
      // OVERLAP REMOVAL
      //    -electrons with dR(e,mu)<0.1 are removed
      //    -lower pT electrons with dR(e,e)<0.1 are removed
      //
      ////////////////////////////////////////////////////////////////////////////
      foreach (DressedLepton& l1, dressed_lepton) {
        bool l_isolated = true;
        foreach (DressedLepton& l2, dressed_lepton) {
          if (l1 != l2 && l2.constituentLepton().abspid() == PID::ELECTRON) {
            double overlapControl_ll= deltaR(l1.constituentLepton(),l2.constituentLepton());
            if (overlapControl_ll < 0.1) {
              l_isolated = false;
              // e/e overlap removal
              if (l1.constituentLepton().abspid() == PID::ELECTRON) {
                if (l1.constituentLepton().pT()>l2.constituentLepton().pT()) {
                  isolated_lepton.push_back(l1);//keep e with highest pT
                } else {
                  isolated_lepton.push_back(l2);//keep e with highest pT
                }
              }
              // e/mu overlap removal
              if (l1.constituentLepton().abspid() == PID::MUON) isolated_lepton.push_back(l1); //keep mu
            }
          }
        }
        if (l_isolated) isolated_lepton.push_back(l1);
      }


      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      // PRESELECTION:
      // "isolated_lepton:"
      //    * electron: pt>20 GeV, |eta|<1.37, 1.52<|eta|<2.47, dR(electron,muon)>0.1
      //    * muon:     pt>20 GeV, |eta|<2.4
      //        *   dR(l,l)>0.1
      //
      // "fiducial_lepton"= isolated_lepton with
      //                * 2 leptons (e or mu)
      //              * leading lepton pt (pT_l1) >25 GeV
      //            * opposite charged leptons
      //
      ///////////////////////////////////////////////////////////////////////////////////////////////////////////
      if (isolated_lepton.size() != 2) vetoEvent;
      sort(isolated_lepton.begin(), isolated_lepton.end(), cmpMomByPt);
      if (isolated_lepton[0].pT() > 25*GeV && threeCharge(isolated_lepton[0]) != threeCharge(isolated_lepton[1])) {
        fiducial_lepton.insert(fiducial_lepton.end(), isolated_lepton.begin(), isolated_lepton.end());
      }
      if (fiducial_lepton.size() == 0) vetoEvent;
      double pT_l1 = fiducial_lepton[0].pT();
      double M_l1l2 = (fiducial_lepton[0].momentum() + fiducial_lepton[1].momentum()).mass();
      double pT_l1l2 = (fiducial_lepton[0].momentum() + fiducial_lepton[1].momentum()).pT();


      /////////////////////////////////////////////////////////////////////////
      // JETS
      //    -"alljets": found by "jetpro" projection && pT()>25 GeV && |y|<4.5
      //    -"vetojets": "alljets"  && dR(electron,jet)>0.3
      //
      /////////////////////////////////////////////////////////////////////////
      Jets alljets, vetojets;
      foreach (const Jet& j, apply<FastJets>(e, "jet").jetsByPt(25)) {
        if (j.absrap() > 4.5 ) continue;
        alljets.push_back(j);
        bool deltaRcontrol = true;
        foreach (DressedLepton& fl,fiducial_lepton) {
          if (fl.constituentLepton().abspid() == PID::ELECTRON) { //electrons
            double deltaRjets = deltaR(fl.constituentLepton().momentum(), j.momentum(), RAPIDITY);
            if (deltaRjets <= 0.3) deltaRcontrol = false; //false if at least one electron is in the overlap region
          }
        }
        if (deltaRcontrol) vetojets.push_back(j);
      }


      /////////////////////////////////////////////////////////////////////////////////////////////////
      // MISSING ETrel
      //    -"mismom": fourvector of invisible momentum found by "met" projection
      //    -"delta_phi": delta phi between mismom and the nearest "fiducial_lepton" or "vetojet"
      //    -"MET_rel": missing transverse energy defined as:
      //            *"mismom"   for "delta_phi" >= (0.5*pi)
      //            *"mismom.pT()*sin(delta_phi)"   for "delta_phi" < (0.5*pi)
      //
      /////////////////////////////////////////////////////////////////////////////////////////////////
      FourMomentum mismom;
      double MET_rel = 0, delta_phi = 0;
      mismom = -met.visibleMomentum();
      vector<double> vL_MET_angle, vJet_MET_angle;
      vL_MET_angle.push_back(fabs(deltaPhi(fiducial_lepton[0].momentum(), mismom)));
      vL_MET_angle.push_back(fabs(deltaPhi(fiducial_lepton[1].momentum(), mismom)));
      foreach (double& lM, vL_MET_angle) if (lM > M_PI) lM = 2*M_PI - lM;

      std::sort(vL_MET_angle.begin(), vL_MET_angle.end());
      if (vetojets.size() == 0) delta_phi = vL_MET_angle[0];
      if (vetojets.size() > 0) {
        foreach (Jet& vj, vetojets) {
          double jet_MET_angle = fabs(deltaPhi(vj.momentum(), mismom));
          if (jet_MET_angle > M_PI) jet_MET_angle = 2*M_PI - jet_MET_angle;
          vJet_MET_angle.push_back(jet_MET_angle);
        }
        std::sort(vJet_MET_angle.begin(), vJet_MET_angle.end());
        if (vL_MET_angle[0] <= vJet_MET_angle[0]) delta_phi = vL_MET_angle[0];
        if (vL_MET_angle[0] > vJet_MET_angle[0]) delta_phi = vJet_MET_angle[0];
      }

      if (delta_phi >= (0.5*M_PI)) delta_phi = 0.5*M_PI;
      MET_rel = mismom.pT()*sin(delta_phi);


      ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      // CUTS
      //        -jetveto: event with at least one vetojet is vetoed
      //        -M_Z: Z mass M_Z=91.1876*GeV
      //
      //    * ee   channel: MET_rel > 45 GeV, M_l1l2 > 15 GeV, |M_l1l2-M_Z| > 15 GeV, jetveto, pT_l1l2 > 30 GeV
      //    * mumu channel: MET_rel > 45 GeV, M_l1l2 > 15 GeV, |M_l1l2-M_Z| > 15 GeV, jetveto, pT_l1l2 > 30 GeV
      //        * emu  channel: MET_rel > 25 GeV, M_l1l2 > 10 GeV, |M_l1l2-M_Z| > 0  GeV, jetveto, pT_l1l2 > 30 GeV
      //
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Get event weight for histo filling
      const double weight = e.weight();

      // ee channel
      if (fiducial_lepton[0].abspid() == PID::ELECTRON && fiducial_lepton[1].abspid() == PID::ELECTRON) {
        if (MET_rel <= 45*GeV) vetoEvent;
        if (M_l1l2 <= 15*GeV) vetoEvent;
        if (fabs(M_l1l2 - 91.1876*GeV) <= 15*GeV) vetoEvent;
        if (vetojets.size() != 0) vetoEvent;
        if (pT_l1l2 <= 30*GeV) vetoEvent;
        _h_Wl1_pT_ee->fill(sqrtS()*GeV, weight);
        _h_Wl1_pT_inclusive->fill(pT_l1, weight);
      }

      // mumu channel
      else if (fiducial_lepton[0].abspid() == PID::MUON && fiducial_lepton[1].abspid() == PID::MUON) {
        if (MET_rel <= 45*GeV) vetoEvent;
        if (M_l1l2 <= 15*GeV) vetoEvent;
        if (fabs(M_l1l2-91.1876*GeV) <= 15*GeV) vetoEvent;
        if (vetojets.size() != 0) vetoEvent;
        if (pT_l1l2 <= 30*GeV) vetoEvent;
        _h_Wl1_pT_mumu->fill(sqrtS()*GeV, weight);
        _h_Wl1_pT_inclusive->fill(pT_l1, weight);
      }

      // emu channel
      else if (fiducial_lepton[0].abspid() != fiducial_lepton[1].abspid()) {
        if (MET_rel <= 25*GeV) vetoEvent;
        if (M_l1l2 <= 10*GeV) vetoEvent;
        if (vetojets.size() != 0) vetoEvent;
        if (pT_l1l2 <= 30*GeV) vetoEvent;
        _h_Wl1_pT_emu->fill(sqrtS()*GeV, weight);
        _h_Wl1_pT_inclusive->fill(pT_l1, weight);
      }
    }


    /// Finalize
    void finalize() {
      const double norm = crossSection()/sumOfWeights()/femtobarn;
      scale(_h_Wl1_pT_ee, norm);
      scale(_h_Wl1_pT_mumu, norm);
      scale(_h_Wl1_pT_emu, norm);
      normalize(_h_Wl1_pT_inclusive, 1);
    }


  private:

    Histo1DPtr _h_Wl1_pT_ee, _h_Wl1_pT_mumu, _h_Wl1_pT_emu, _h_Wl1_pT_inclusive;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2013_I1190187);

}
