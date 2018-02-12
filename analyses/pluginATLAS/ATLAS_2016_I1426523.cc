// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief Measurement of the WZ production cross section at 8 TeV
  class ATLAS_2016_I1426523 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_I1426523);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Lepton cuts
      Cut FS_Zlept = Cuts::abseta < 2.5 && Cuts::pT > 15*GeV;

      const FinalState fs;
      Cut fs_z = Cuts::abseta < 2.5 && Cuts::pT > 15*GeV;
      Cut fs_j = Cuts::abseta < 4.5 && Cuts::pT > 25*GeV;

      // Get photons to dress leptons
      PromptFinalState photons(Cuts::abspid == PID::PHOTON);

      // Electrons and muons in Fiducial PS
      PromptFinalState leptons(fs_z && (Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON));
      leptons.acceptTauDecays(false);
      DressedLeptons dressedleptons(photons, leptons, 0.1, FS_Zlept, true);
      addProjection(dressedleptons, "DressedLeptons");

      // Electrons and muons in Total PS
      PromptFinalState leptons_total(Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON);
      leptons_total.acceptTauDecays(false);
      DressedLeptons dressedleptonsTotal(photons, leptons_total, 0.1, Cuts::open(), true);
      addProjection(dressedleptonsTotal, "DressedLeptonsTotal");

      // Neutrinos
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
      _h["eee"]        = bookHisto1D(1, 1, 1);
      _h["mee"]        = bookHisto1D(1, 1, 2);
      _h["emm"]        = bookHisto1D(1, 1, 3);
      _h["mmm"]        = bookHisto1D(1, 1, 4);
      _h["fid"]        = bookHisto1D(1, 1, 5);
      _h["eee_Plus"]   = bookHisto1D(2, 1, 1);
      _h["mee_Plus"]   = bookHisto1D(2, 1, 2);
      _h["emm_Plus"]   = bookHisto1D(2, 1, 3);
      _h["mmm_Plus"]   = bookHisto1D(2, 1, 4);
      _h["fid_Plus"]   = bookHisto1D(2, 1, 5);
      _h["eee_Minus"]  = bookHisto1D(3, 1, 1);
      _h["mee_Minus"]  = bookHisto1D(3, 1, 2);
      _h["emm_Minus"]  = bookHisto1D(3, 1, 3);
      _h["mmm_Minus"]  = bookHisto1D(3, 1, 4);
      _h["fid_Minus"]  = bookHisto1D(3, 1, 5);
      _h["total"]      = bookHisto1D(5, 1, 1);
      _h["Njets"]      = bookHisto1D(27, 1, 1);
      _h["Njets_norm"] = bookHisto1D(41, 1, 1);

      bookHandler("ZpT",	             12);
      bookHandler("ZpT_Plus",          13);
      bookHandler("ZpT_Minus",         14);
      bookHandler("WpT",	             15);
      bookHandler("WpT_Plus",          16);
      bookHandler("WpT_Minus",         17);
      bookHandler("mTWZ",              18);
      bookHandler("mTWZ_Plus",         19);
      bookHandler("mTWZ_Minus",        20);
      bookHandler("pTv",               21);
      bookHandler("pTv_Plus",          22);
      bookHandler("pTv_Minus",         23);
      bookHandler("Deltay",	           24);
      bookHandler("Deltay_Plus",       25);
      bookHandler("Deltay_Minus",      26);
      bookHandler("mjj",               28);
      bookHandler("Deltayjj",          29);
      bookHandler("ZpT_norm",          30);
      bookHandler("ZpT_Plus_norm",     31);
      bookHandler("ZpT_Minus_norm",    32);
      bookHandler("WpT_norm",          33);
      bookHandler("mTWZ_norm",         34);
      bookHandler("pTv_norm", 	       35);
      bookHandler("pTv_Plus_norm",	   36);
      bookHandler("pTv_Minus_norm",	   37);
      bookHandler("Deltay_norm",	     38);
      bookHandler("Deltay_Minus_norm", 39);
      bookHandler("Deltay_Plus_norm",  40);
      bookHandler("mjj_norm",          42);
      bookHandler("Deltayjj_norm",     43);
    }

    void bookHandler(const string& tag, size_t ID) {
      _s[tag] = bookScatter2D(ID, 1, 1);
      const string code1 = makeAxisCode(ID, 1, 1);
      const string code2 = makeAxisCode(ID, 1, 2);
      _h[tag] = bookHisto1D(code2, refData(code1));
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      /// @todo Do the event by event analysis here
      const double weight = event.weight();

      const vector<DressedLepton>& dressedleptons = apply<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
      const vector<DressedLepton>& dressedleptonsTotal = apply<DressedLeptons>(event, "DressedLeptonsTotal").dressedLeptons();
      const Particles& neutrinos = apply<PromptFinalState>(event, "Neutrinos").particlesByPt();
      Jets jets = apply<JetAlg>(event, "Jets").jetsByPt( (Cuts::abseta < 4.5) && (Cuts::pT > 25*GeV) );

      if ((dressedleptonsTotal.size()<3) || (neutrinos.size()<1)) vetoEvent;

      //---Total PS: assign leptons to W and Z bosons using Resonant shape algorithm
      // NB: This resonant shape algorithm assumes the Standard Model and can therefore
      // NOT be used for reinterpretation in terms of new-physics models.

      int i, j, k;
      double MassZ01 = 0., MassZ02 = 0., MassZ12 = 0.;
      double MassW0 = 0., MassW1 = 0., MassW2 = 0.;
      double WeightZ1, WeightZ2, WeightZ3;
      double WeightW1, WeightW2, WeightW3;
      double M1, M2, M3;
      double WeightTotal1, WeightTotal2, WeightTotal3;

      //try Z pair of leptons 01
      if ( (dressedleptonsTotal[0].pid()==-(dressedleptonsTotal[1].pid())) && (dressedleptonsTotal[2].abspid()==neutrinos[0].abspid()-1)){
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
      if( (M2 < M1 && M2 < M3) || (MassZ02 != 0 && MassW1 != 0 && MassZ01 == 0 && MassZ12 == 0) ){
        i = 0; j = 2; k = 1;
      }
      if( (M3 < M1 && M3 < M2) || (MassZ12 != 0 && MassW0 != 0 && MassZ01 == 0 && MassZ02 == 0) ){
        i = 1; j = 2; k = 0;
      }

      FourMomentum ZbosonTotal   = dressedleptonsTotal[i].momentum()+dressedleptonsTotal[j].momentum();

      if ( (ZbosonTotal.mass() >= 66*GeV) && (ZbosonTotal.mass() <= 116*GeV) ) _h["total"]->fill(8000, weight);

      //---end Total PS


      //---Fiducial PS: assign leptons to W and Z bosons using Resonant shape algorithm
      if (dressedleptons.size() < 3 || neutrinos.size() < 1)  vetoEvent;

      int EventType = -1;
      int Nel = 0, Nmu = 0;

      for (const DressedLepton& l : dressedleptons) {
        if (l.abspid() == 11)  ++Nel;
        if (l.abspid() == 13)  ++Nmu;
      }

      if ( Nel == 3  && Nmu==0 )  EventType = 3;
      if ( Nel == 2  && Nmu==1 )  EventType = 2;
      if ( Nel == 1  && Nmu==2 )  EventType = 1;
      if ( Nel == 0  && Nmu==3 )  EventType = 0;

      int EventCharge = -dressedleptons[0].charge()*dressedleptons[1].charge()*dressedleptons[2].charge();

      MassZ01 = 0; MassZ02 = 0; MassZ12 = 0;
      MassW0 = 0;  MassW1 = 0;  MassW2 = 0;

      //try Z pair of leptons 01
      if ( (dressedleptons[0].pid()==-(dressedleptons[1].pid())) && (dressedleptons[2].abspid()==neutrinos[0].abspid()-1)){
        MassZ01 = (dressedleptons[0].momentum()+dressedleptons[1].momentum()).mass();
        MassW2 = (dressedleptons[2].momentum()+neutrinos[0].momentum()).mass();
      }
      //try Z pair of leptons 02
      if ( (dressedleptons[0].pid()==-(dressedleptons[2].pid())) && (dressedleptons[1].abspid()==neutrinos[0].abspid()-1)){
        MassZ02 = (dressedleptons[0].momentum()+dressedleptons[2].momentum()).mass();
        MassW1 = (dressedleptons[1].momentum()+neutrinos[0].momentum()).mass();
      }
      //try Z pair of leptons 12
      if ( (dressedleptons[1].pid()==-(dressedleptons[2].pid())) && (dressedleptons[0].abspid()==neutrinos[0].abspid()-1)){
        MassZ12 = (dressedleptons[1].momentum()+dressedleptons[2].momentum()).mass();
        MassW0 = (dressedleptons[0].momentum()+neutrinos[0].momentum()).mass();
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
      if( (M2 < M1 && M2 < M3) || (MassZ02 != 0 && MassW1 != 0 && MassZ01 == 0 && MassZ12 == 0) ){
        i = 0; j = 2; k = 1;
      }
      if( (M3 < M1 && M3 < M2) || (MassZ12 != 0 && MassW0 != 0 && MassZ01 == 0 && MassZ02 == 0) ){
        i = 1; j = 2; k = 0;
      }

      FourMomentum Zlepton1 = dressedleptons[i].momentum();
      FourMomentum Zlepton2 = dressedleptons[j].momentum();
      FourMomentum Wlepton  = dressedleptons[k].momentum();
      FourMomentum Zboson   = dressedleptons[i].momentum()+dressedleptons[j].momentum();
      FourMomentum Wboson   = dressedleptons[k].momentum()+neutrinos[0].momentum();

      double Wboson_mT = sqrt( 2 * Wlepton.pT() * neutrinos[0].pt() * (1 - cos(deltaPhi(Wlepton, neutrinos[0]))) )/GeV;

      if (fabs(Zboson.mass()-MZ_PDG)>=10.)  vetoEvent;
      if (Wboson_mT<=30.)                   vetoEvent;
      if (Wlepton.pT()<=20.)                vetoEvent;
      if (deltaR(Zlepton1,Zlepton2) < 0.2)  vetoEvent;
      if (deltaR(Zlepton1,Wlepton)  < 0.3)  vetoEvent;
      if (deltaR(Zlepton2,Wlepton)  < 0.3)  vetoEvent;

      double WZ_pt = Zlepton1.pt() + Zlepton2.pt() + Wlepton.pt() + neutrinos[0].pt();
      double WZ_px = Zlepton1.px() + Zlepton2.px() + Wlepton.px() + neutrinos[0].px();
      double WZ_py = Zlepton1.py() + Zlepton2.py() + Wlepton.py() + neutrinos[0].py();
      double mTWZ = sqrt( pow(WZ_pt, 2) - ( pow(WZ_px, 2) + pow(WZ_py,2) ) )/GeV;

      double AbsDeltay = fabs(Zboson.rapidity()-Wlepton.rapidity());

      if (EventType == 3) _h["eee"]->fill(8000., weight);
      if (EventType == 2) _h["mee"]->fill(8000., weight);
      if (EventType == 1) _h["emm"]->fill(8000., weight);
      if (EventType == 0) _h["mmm"]->fill(8000., weight);
      _h["fid"]->fill(8000., weight);

      if (EventCharge == 1) {

        if (EventType == 3) _h["eee_Plus"]->fill(8000., weight);
        if (EventType == 2) _h["mee_Plus"]->fill(8000., weight);
        if (EventType == 1) _h["emm_Plus"]->fill(8000., weight);
        if (EventType == 0) _h["mmm_Plus"]->fill(8000., weight);
        _h["fid_Plus"]->fill(8000., weight);

        _h["Deltay_Plus"]->fill(AbsDeltay, weight);
        _h["Deltay_Plus_norm"]->fill(AbsDeltay, weight);
        fillWithOverflow("ZpT_Plus", Zboson.pT()/GeV, 220, weight);
        fillWithOverflow("WpT_Plus", Wboson.pT()/GeV, 220, weight);
        fillWithOverflow("mTWZ_Plus", mTWZ, 600, weight);
        fillWithOverflow("pTv_Plus", neutrinos[0].pt(), 90, weight);
        fillWithOverflow("ZpT_Plus_norm", Zboson.pT()/GeV, 220, weight);
        fillWithOverflow("pTv_Plus_norm", neutrinos[0].pt()/GeV, 90, weight);

      } else {

        if (EventType == 3) _h["eee_Minus"]->fill(8000., weight);
        if (EventType == 2) _h["mee_Minus"]->fill(8000., weight);
        if (EventType == 1) _h["emm_Minus"]->fill(8000., weight);
        if (EventType == 0) _h["mmm_Minus"]->fill(8000., weight);
        _h["fid_Minus"]->fill(8000., weight);

        _h["Deltay_Minus"]->fill(AbsDeltay, weight);
        _h["Deltay_Minus_norm"]->fill(AbsDeltay, weight);
        fillWithOverflow("ZpT_Minus", Zboson.pT()/GeV, 220, weight);
        fillWithOverflow("WpT_Minus", Wboson.pT()/GeV, 220, weight);
        fillWithOverflow("mTWZ_Minus", mTWZ, 600, weight);
        fillWithOverflow("pTv_Minus", neutrinos[0].pt()/GeV, 90, weight);
        fillWithOverflow("ZpT_Minus_norm", Zboson.pT()/GeV, 220, weight);
        fillWithOverflow("pTv_Minus_norm", neutrinos[0].pt()/GeV, 90, weight);

      }

      fillWithOverflow("ZpT", Zboson.pT()/GeV, 220, weight);
      fillWithOverflow("WpT", Wboson.pT()/GeV, 220, weight);
      fillWithOverflow("mTWZ", mTWZ, 600, weight);
      fillWithOverflow("pTv", neutrinos[0].pt()/GeV, 90, weight);

      _h["Deltay"]->fill(AbsDeltay, weight);

      fillWithOverflow("Njets", jets.size(), 5, weight);
      fillWithOverflow("Njets_norm", jets.size(), 5, weight);
      fillWithOverflow("ZpT_norm", Zboson.pT()/GeV, 220, weight);
      fillWithOverflow("WpT_norm", Wboson.pT()/GeV, 220, weight);
      fillWithOverflow("mTWZ_norm", mTWZ, 600, weight);
      fillWithOverflow("pTv_norm", neutrinos[0].pt()/GeV, 90, weight);

      _h["Deltay_norm"]->fill(AbsDeltay, weight);

      if (jets.size()>1) {
        double mjj = (jets[0].momentum()+jets[1].momentum()).mass()/GeV;
        fillWithOverflow("mjj",      mjj, 800, weight);
        fillWithOverflow("mjj_norm", mjj, 800, weight);
        double DeltaYjj = fabs(jets[0].rapidity()-jets[1].rapidity());
        fillWithOverflow("Deltayjj",      DeltaYjj, 5, weight);
        fillWithOverflow("Deltayjj_norm", DeltaYjj, 5, weight);
      }

    }


    void fillWithOverflow(const string& tag, const double value, const double overflow, const double weight){
      if (value < overflow) _h[tag]->fill(value,   weight);
      else _h[tag]->fill(overflow - 0.45,   weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double xs_pb(crossSection() / picobarn);
      const double xs_fb(crossSection() / femtobarn);
      const double sumw(sumOfWeights());
      MSG_INFO("Cross-Section/pb: " << xs_pb      );
      MSG_INFO("Cross-Section/fb: " << xs_fb      );
      MSG_INFO("Sum of weights  : " << sumw       );
      MSG_INFO("nEvents         : " << numEvents());

      const double sf_pb(xs_pb / sumw);
      const double sf_fb(xs_fb / sumw);

      MSG_INFO("sf_pb         : " << sf_pb);
      MSG_INFO("sf_fb         : " << sf_fb);

      float totalBR= 4*0.1086*0.033658; // W and Z leptonic branching fractions

      for (map<string, Histo1DPtr>::iterator it = _h.begin(); it != _h.end(); ++it) {
        if (it->first.find("total") != string::npos)        scale(it->second, sf_pb/totalBR);
        else if (it->first.find("norm") != string::npos)    normalize(it->second);
        else if (it->first.find("fid") != string::npos)     scale(it->second, sf_fb/4.);
        else if (it->first.find("Njets") != string::npos)   scale(it->second, sf_fb/4.);
        else if (it->first.find("ZpT") != string::npos)     scale(it->second, sf_fb/4.);
        else if (it->first.find("WpT") != string::npos)     scale(it->second, sf_fb/4.);
        else if (it->first.find("mTWZ") != string::npos)    scale(it->second, sf_fb/4.);
        else if (it->first.find("pTv") != string::npos)     scale(it->second, sf_fb/4.);
        else if (it->first.find("Deltay") != string::npos)  scale(it->second, sf_fb/4.);
        else if (it->first.find("mjj") != string::npos)     scale(it->second, sf_fb/4.);
        else                                                scale(it->second, sf_fb);
      }
      for (map<string, Scatter2DPtr>::iterator it = _s.begin(); it != _s.end(); ++it) {
        makeScatterWithoutDividingByBinwidth(it->first);
        removeAnalysisObject(_h[it->first]);
      }
    }

    void makeScatterWithoutDividingByBinwidth(const string& tag) {
      vector<Point2D> points;
      //size_t nBins = _dummy->numBins();
      for (const HistoBin1D &bin : _h[tag]->bins()) {
        double  x = bin.midpoint();
        double  y = bin.sumW();
        double ex = bin.xWidth()/2;
        double ey = sqrt(bin.sumW2());
        points.push_back(Point2D(x, y, ex, ey));
      }
      _s[tag]->addPoints(points);
    }


    //@}


  private:

    /// @name Histograms
    //@{

     map<string, Histo1DPtr> _h;
     map<string, Scatter2DPtr> _s;

     //@}

     double MZ_PDG = 91.1876;
     double MW_PDG = 83.385;
     double GammaZ_PDG = 2.4952;
     double GammaW_PDG = 2.085;

  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1426523);

}
