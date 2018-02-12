// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief Measurement of the WZ production cross section at 13 TeV
  class ATLAS_2016_I1469071 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1469071);

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Lepton cuts
      Cut FS_Zlept = Cuts::abseta < 2.5 && Cuts::pT > 15*GeV;

      FinalState fs;
      Cut fs_z = Cuts::abseta < 2.5 && Cuts::pT > 15*GeV;
      Cut fs_j = Cuts::abseta < 4.5 && Cuts::pT > 25*GeV;

      // Get photons to dress leptons
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);

      // Electrons and muons in Fiducial PS
      PromptFinalState leptons(FinalState(fs_z && (Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON)));
      leptons.acceptTauDecays(false);
      DressedLeptons dressedleptons(photons, leptons, 0.1, FS_Zlept, true);
      addProjection(dressedleptons, "DressedLeptons");

      // Electrons and muons in Total PS
      PromptFinalState leptons_total(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      leptons_total.acceptTauDecays(false);
      DressedLeptons dressedleptonsTotal(photons, leptons_total, 0.1, Cuts::open(), true);
      addProjection(dressedleptonsTotal, "DressedLeptonsTotal");

      // Promot neutrinos (yikes!)
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(false);
      declare(neutrinos, "Neutrinos");
      MSG_WARNING("\033[91;1mLIMITED VALIDITY - check info file for details!\033[m");

      // Jets
      VetoedFinalState veto;
      veto.addVetoOnThisFinalState(dressedleptons);
      FastJets jets(veto, FastJets::ANTIKT, 0.4);
      declare(jets, "Jets");

      // Book histograms
      _h_eee          = bookHisto1D(1, 1, 1);
      _h_mee          = bookHisto1D(1, 1, 2);
      _h_emm          = bookHisto1D(1, 1, 3);
      _h_mmm          = bookHisto1D(1, 1, 4);
      _h_fid          = bookHisto1D(1, 1, 5);
      _h_eee_Plus     = bookHisto1D(2, 1, 1);
      _h_mee_Plus     = bookHisto1D(2, 1, 2);
      _h_emm_Plus     = bookHisto1D(2, 1, 3);
      _h_mmm_Plus     = bookHisto1D(2, 1, 4);
      _h_fid_Plus     = bookHisto1D(2, 1, 5);
      _h_eee_Minus    = bookHisto1D(3, 1, 1);
      _h_mee_Minus    = bookHisto1D(3, 1, 2);
      _h_emm_Minus    = bookHisto1D(3, 1, 3);
      _h_mmm_Minus    = bookHisto1D(3, 1, 4);
      _h_fid_Minus    = bookHisto1D(3, 1, 5);
      _h_total        = bookHisto1D(6, 1, 1);
      _h_Njets        = bookHisto1D(8, 1, 1);

    }


    void analyze(const Event& event) {

      const double weight = event.weight();

      const vector<DressedLepton>& dressedleptons = apply<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
      const vector<DressedLepton>& dressedleptonsTotal = apply<DressedLeptons>(event, "DressedLeptonsTotal").dressedLeptons();
      const Particles& neutrinos = apply<PromptFinalState>(event, "Neutrinos").particlesByPt();
      Jets jets = apply<JetAlg>(event, "Jets").jetsByPt( (Cuts::abseta < 4.5) && (Cuts::pT > 25*GeV) );

      if (dressedleptonsTotal.size() < 3 || neutrinos.size() < 1) vetoEvent;

      //---Total PS: assign leptons to W and Z bosons using Resonant shape algorithm
      // NB: This resonant shape algorithm assumes the Standard Model and can therefore
      //     NOT be used for any kind of reinterpretation in terms of new-physics models..

      int i, j, k;
      double MassZ01 = 0., MassZ02 = 0., MassZ12 = 0.;
      double MassW0 = 0., MassW1 = 0., MassW2 = 0.;
      double WeightZ1, WeightZ2, WeightZ3;
      double WeightW1, WeightW2, WeightW3;
      double M1, M2, M3;
      double WeightTotal1, WeightTotal2, WeightTotal3;

      //try Z pair of leptons 01
      if ( (dressedleptonsTotal[0].pid() ==-(dressedleptonsTotal[1].pid())) && (dressedleptonsTotal[2].abspid()==neutrinos[0].abspid()-1)){
        MassZ01 = (dressedleptonsTotal[0].momentum()+dressedleptonsTotal[1].momentum()).mass();
        MassW2 = (dressedleptonsTotal[2].momentum()+neutrinos[0].momentum()).mass();
      }
      //try Z pair of leptons 02
      if ( (dressedleptonsTotal[0].pid()==-(dressedleptonsTotal[2].pid())) && (dressedleptonsTotal[1].abspid()==neutrinos[0].abspid()-1)){
        MassZ02 = (dressedleptonsTotal[0].momentum()+dressedleptonsTotal[2].momentum()).mass();
        MassW1 = (dressedleptonsTotal[1].momentum()+neutrinos[0].momentum()).mass();
      }
      //try Z pair of leptons 12
      if ( (dressedleptonsTotal[1].pid()==-(dressedleptonsTotal[2].pid())) && (dressedleptonsTotal[0].abspid()==neutrinos[0].abspid()-1)){
        MassZ12 = (dressedleptonsTotal[1].momentum()+dressedleptonsTotal[2].momentum()).mass();
        MassW0 = (dressedleptonsTotal[0].momentum()+neutrinos[0].momentum()).mass();
      }
      WeightZ1 = 1/(pow(MassZ01*MassZ01 - MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW1 = 1/(pow(MassW2*MassW2 - MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal1 = WeightZ1*WeightW1;
      M1 = -1*WeightTotal1;

      WeightZ2 = 1/(pow(MassZ02*MassZ02- MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW2 = 1/(pow(MassW1*MassW1- MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal2 = WeightZ2*WeightW2;
      M2 = -1*WeightTotal2;

      WeightZ3 = 1/(pow(MassZ12*MassZ12 - MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW3 = 1/(pow(MassW0*MassW0 - MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal3 = WeightZ3*WeightW3;
      M3 = -1*WeightTotal3;

      if( (M1 < M2 && M1 < M3) || (MassZ01 != 0 && MassW2 != 0 && MassZ02 == 0 && MassZ12 == 0) ){
        i = 0; j = 1; k = 2;
      }
      if((M2 < M1 && M2 < M3) || (MassZ02 != 0 && MassW1 != 0 && MassZ01 == 0 && MassZ12 == 0) ){
        i = 0; j = 2; k = 1;
      }
      if((M3 < M1 && M3 < M2) || (MassZ12 != 0 && MassW0 != 0 && MassZ01 == 0 && MassZ02 == 0) ){
        i = 1; j = 2; k = 0;
      }

      FourMomentum ZbosonTotal   = dressedleptonsTotal[i].momentum()+dressedleptonsTotal[j].momentum();
      if ( ZbosonTotal.mass() >= 66*GeV && ZbosonTotal.mass() <= 116*GeV )  _h_total->fill(13000, weight);

      //---end Total PS


      //---Fiducial PS: assign leptons to W and Z bosons using Resonant shape algorithm
      if (dressedleptons.size() < 3)  vetoEvent;

      int EventType = -1;
      int Nel = 0, Nmu = 0;

      for (const DressedLepton& l : dressedleptons) {
        if (l.abspid() == 11)  ++Nel;
        if (l.abspid() == 13)  ++Nmu;
      }

      if ( (Nel == 3)  && (Nmu==0) )  EventType = 3;
      if ( (Nel == 2)  && (Nmu==1) )  EventType = 2;
      if ( (Nel == 1)  && (Nmu==2) )  EventType = 1;
      if ( (Nel == 0)  && (Nmu==3) )  EventType = 0;

      int EventCharge = -dressedleptons[0].charge() * dressedleptons[1].charge() * dressedleptons[2].charge();

      MassZ01 = 0; MassZ02 = 0; MassZ12 = 0;
      MassW0 = 0;  MassW1 = 0;  MassW2 = 0;

      // try Z pair of leptons 01
      if (dressedleptons[0].pid() == -dressedleptons[1].pid()) {
        MassZ01 = (dressedleptons[0].momentum() + dressedleptons[1].momentum()).mass();
        MassW2 = (dressedleptons[2].momentum() + neutrinos[0].momentum()).mass();
      }
      // try Z pair of leptons 02
      if (dressedleptons[0].pid() == -dressedleptons[2].pid()) {
        MassZ02 = (dressedleptons[0].momentum() + dressedleptons[2].momentum()).mass();
        MassW1 = (dressedleptons[1].momentum() + neutrinos[0].momentum()).mass();
      }
      // try Z pair of leptons 12
      if (dressedleptons[1].pid() == -dressedleptons[2].pid()) {
        MassZ12 = (dressedleptons[1].momentum() + dressedleptons[2].momentum()).mass();
        MassW0 = (dressedleptons[0].momentum() + neutrinos[0].momentum()).mass();
      }
      WeightZ1 = 1/(pow(MassZ01*MassZ01 - MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW1 = 1/(pow(MassW2*MassW2 - MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal1 = WeightZ1*WeightW1;
      M1 = -1*WeightTotal1;

      WeightZ2 = 1/(pow(MassZ02*MassZ02- MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW2 = 1/(pow(MassW1*MassW1- MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal2 = WeightZ2*WeightW2;
      M2 = -1*WeightTotal2;

      WeightZ3 = 1/(pow(MassZ12*MassZ12 - MZ_PDG*MZ_PDG,2) + pow(MZ_PDG*GammaZ_PDG,2));
      WeightW3 = 1/(pow(MassW0*MassW0 - MW_PDG*MW_PDG,2) + pow(MW_PDG*GammaW_PDG,2));
      WeightTotal3 = WeightZ3*WeightW3;
      M3 = -1*WeightTotal3;

      if( (M1 < M2 && M1 < M3) || (MassZ01 != 0 && MassW2 != 0 && MassZ02 == 0 && MassZ12 == 0) ) {
        i = 0; j = 1; k = 2;
      }
      if((M2 < M1 && M2 < M3) || (MassZ02 != 0 && MassW1 != 0 && MassZ01 == 0 && MassZ12 == 0) ) {
        i = 0; j = 2; k = 1;
      }
      if((M3 < M1 && M3 < M2) || (MassZ12 != 0 && MassW0 != 0 && MassZ01 == 0 && MassZ02 == 0) ) {
        i = 1; j = 2; k = 0;
      }

      FourMomentum Zlepton1 = dressedleptons[i].momentum();
      FourMomentum Zlepton2 = dressedleptons[j].momentum();
      FourMomentum Wlepton  = dressedleptons[k].momentum();
      FourMomentum Zboson   = dressedleptons[i].momentum()+dressedleptons[j].momentum();

      double Wboson_mT = sqrt( 2 * Wlepton.pT() * neutrinos[0].pt() * (1 - cos(deltaPhi(Wlepton, neutrinos[0]))) );

      if (fabs(Zboson.mass()/GeV - MZ_PDG) >= 10.) vetoEvent;
      if (Wboson_mT <= 30*GeV)                     vetoEvent;
      if (Wlepton.pT() <= 20*GeV)                  vetoEvent;
      if (deltaR(Zlepton1, Zlepton2) < 0.2)        vetoEvent;
      if (deltaR(Zlepton1, Wlepton)  < 0.3)        vetoEvent;
      if (deltaR(Zlepton2, Wlepton)  < 0.3)        vetoEvent;

      if (EventType == 3)  _h_eee->fill(13000., weight);
      if (EventType == 2)  _h_mee->fill(13000., weight);
      if (EventType == 1)  _h_emm->fill(13000., weight);
      if (EventType == 0)  _h_mmm->fill(13000., weight);
      _h_fid->fill(13000, weight);

      if (EventCharge == 1) {
        if (EventType == 3)  _h_eee_Plus->fill(13000., weight);
        if (EventType == 2)  _h_mee_Plus->fill(13000., weight);
        if (EventType == 1)  _h_emm_Plus->fill(13000., weight);
        if (EventType == 0)  _h_mmm_Plus->fill(13000., weight);
        _h_fid_Plus->fill(13000, weight);
      } else {
        if (EventType == 3)  _h_eee_Minus->fill(13000., weight);
        if (EventType == 2)  _h_mee_Minus->fill(13000., weight);
        if (EventType == 1)  _h_emm_Minus->fill(13000., weight);
        if (EventType == 0)  _h_mmm_Minus->fill(13000., weight);
        _h_fid_Minus->fill(13000, weight);
      }

      if (jets.size() < 4)  _h_Njets->fill(jets.size(), weight);
      else  _h_Njets->fill(4, weight);

    }


    void finalize() {

      // Print summary info
      const double xs_pb(crossSection() / picobarn);
      const double xs_fb(crossSection() / femtobarn);
      const double sumw(sumOfWeights());
      const double sf_pb(xs_pb / sumw);
      const double sf_fb(xs_fb / sumw);

      const float totalBR= 4*0.1086*0.033658; // W and Z leptonic branching fractions

      scale(_h_fid,       sf_fb/4.);
      scale(_h_eee,       sf_fb);
      scale(_h_mee,       sf_fb);
      scale(_h_emm,       sf_fb);
      scale(_h_mmm,       sf_fb);
      scale(_h_fid_Plus,  sf_fb/4.);
      scale(_h_eee_Plus,  sf_fb);
      scale(_h_mee_Plus,  sf_fb);
      scale(_h_emm_Plus,  sf_fb);
      scale(_h_mmm_Plus,  sf_fb);
      scale(_h_fid_Minus, sf_fb/4.);
      scale(_h_eee_Minus, sf_fb);
      scale(_h_mee_Minus, sf_fb);
      scale(_h_emm_Minus, sf_fb);
      scale(_h_mmm_Minus, sf_fb);
      scale(_h_Njets, sf_fb/4.);
      scale(_h_total, sf_pb/totalBR);

    }

    //@}


  private:


    /// @name Histograms
    //@{
    Histo1DPtr _h_eee;
    Histo1DPtr _h_mee;
    Histo1DPtr _h_emm;
    Histo1DPtr _h_mmm;
    Histo1DPtr _h_fid;
    Histo1DPtr _h_eee_Plus;
    Histo1DPtr _h_mee_Plus;
    Histo1DPtr _h_emm_Plus;
    Histo1DPtr _h_mmm_Plus;
    Histo1DPtr _h_fid_Plus;
    Histo1DPtr _h_eee_Minus;
    Histo1DPtr _h_mee_Minus;
    Histo1DPtr _h_emm_Minus;
    Histo1DPtr _h_mmm_Minus;
    Histo1DPtr _h_fid_Minus;
    Histo1DPtr _h_total;
    Histo1DPtr _h_Njets;

    //@}

    double MZ_PDG = 91.1876;
    double MW_PDG = 83.385;
    double GammaZ_PDG = 2.4952;
    double GammaW_PDG = 2.085;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1469071);
}
