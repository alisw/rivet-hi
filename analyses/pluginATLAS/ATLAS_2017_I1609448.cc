// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  /// ATLAS pTmiss+jets cross-section ratios at 13 TeV
  class ATLAS_2017_I1609448 : public Analysis {
  public:

    /// Constructor
    ATLAS_2017_I1609448(string name="ATLAS_2017_I1609448")
      : Analysis(name)
    {
      _mode = 0; // using Z -> nunu channel by default
      setNeedsCrossSection(true);
    }


    struct HistoHandler {
      Histo1DPtr histo;
      Scatter2DPtr scatter;
      unsigned int d, x, y;

      void fill(double value, double weight) {
        histo->fill(value, weight);
      }
    };


    /// Initialize
    void init() {

      // Prompt photons
      PromptFinalState photon_fs(Cuts::abspid == PID::PHOTON && Cuts::abseta < 4.9);
      // Prompt electrons
      PromptFinalState el_fs(Cuts::abseta < 4.9 && Cuts::abspid == PID::ELECTRON);
      // Prompt muons
      PromptFinalState mu_fs(Cuts::abseta < 4.9 && Cuts::abspid == PID::MUON);

      // Dressed leptons
      Cut lep_cuts = Cuts::pT > 7*GeV && Cuts::abseta < 2.5;
      DressedLeptons dressed_leps(photon_fs, (_mode == 2 ? el_fs : mu_fs), 0.1, lep_cuts);
      declare(dressed_leps, "DressedLeptons");

      // In-acceptance leptons for lepton veto
      PromptFinalState veto_lep_fs(Cuts::abseta < 4.9 && (Cuts::abspid == PID::ELECTRON || Cuts::abspid == PID::MUON));
      veto_lep_fs.acceptTauDecays();
      veto_lep_fs.acceptMuonDecays();
      DressedLeptons veto_lep(photon_fs, veto_lep_fs, 0.1, lep_cuts);
      declare(veto_lep, "VetoLeptons");

      // MET
      VetoedFinalState met_fs(!(Cuts::abseta > 2.5 && Cuts::abspid == PID::MUON));
      if (_mode) met_fs.addVetoOnThisFinalState(dressed_leps);
      declare(MissingMomentum(met_fs), "MET");

      // Jet collection
      FastJets jets(FinalState(Cuts::abseta < 4.9), FastJets::ANTIKT, 0.4, JetAlg::NO_MUONS, JetAlg::NO_INVISIBLES);
      declare(jets, "Jets");

      _h["met_mono"] = bookHandler(1, 1, 2);
      _h["met_vbf" ] = bookHandler(2, 1, 2);
      _h["mjj_vbf" ] = bookHandler(3, 1, 2);
      _h["dphijj_vbf"] = bookHandler(4, 1, 2);
    }


    HistoHandler bookHandler(unsigned int id_d, unsigned int id_x, unsigned int id_y) {
      HistoHandler dummy;
      if (_mode < 2) {  // numerator mode
        const string histName = "_" + makeAxisCode(id_d, id_x, id_y);
        dummy.histo = bookHisto1D(histName, refData(id_d, id_x, id_y)); // hidden auxiliary output
        dummy.scatter = bookScatter2D(id_d, id_x, id_y - 1, true); // ratio
        dummy.d = id_d;
        dummy.x = id_x;
        dummy.y = id_y;
      } else {
        dummy.histo = bookHisto1D(id_d, id_x, 4); // denominator mode
      }
      return dummy;
    }


    bool isBetweenJets(const Jet& probe, const Jet& boundary1, const Jet& boundary2) {
      const double y_p = probe.rapidity();
      const double y_b1 = boundary1.rapidity();
      const double y_b2 = boundary2.rapidity();
      const double y_min = std::min(y_b1, y_b2);
      const double y_max = std::max(y_b1, y_b2);
      return (y_p > y_min && y_p < y_max);
    }


    int centralJetVeto(Jets& jets) {
      if (jets.size() < 2) return 0;
      const Jet bj1 = jets.at(0);
      const Jet bj2 = jets.at(1);

      // Start loop at the 3rd hardest pT jet
      int n_between = 0;
      for (size_t i = 2; i < jets.size(); ++i) {
        const Jet j = jets.at(i);
        if (isBetweenJets(j, bj1, bj2) && j.pT() > 25*GeV)  ++n_between;
      }
      return n_between;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      // Require 0 (Znunu) or 2 (Zll) dressed leptons
      bool isZll = bool(_mode);
      const vector<DressedLepton> &vetoLeptons = applyProjection<DressedLeptons>(event, "VetoLeptons").dressedLeptons();
      const vector<DressedLepton> &all_leps = applyProjection<DressedLeptons>(event, "DressedLeptons").dressedLeptons();
      if (!isZll && vetoLeptons.size())    vetoEvent;
      if ( isZll && all_leps.size() != 2)  vetoEvent;

      vector<DressedLepton> leptons;
      bool pass_Zll = true;
      if (isZll) {
        // Sort dressed leptons by pT
        if (all_leps[0].pt() > all_leps[1].pt()) {
          leptons.push_back(all_leps[0]);
          leptons.push_back(all_leps[1]);
        } else {
          leptons.push_back(all_leps[1]);
          leptons.push_back(all_leps[0]);
        }
        // Leading lepton pT cut
        pass_Zll &= leptons[0].pT() > 80*GeV;
        // Opposite-charge requirement
        pass_Zll &= threeCharge(leptons[0]) + threeCharge(leptons[1]) == 0;
        // Z-mass requirement
        const double Zmass = (leptons[0].mom() + leptons[1].mom()).mass();
        pass_Zll &= (Zmass >= 66*GeV && Zmass <= 116*GeV);
      }
      if (!pass_Zll)  vetoEvent;


      // Get jets and remove those within dR = 0.5 of a dressed lepton
      Jets jets = applyProjection<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::absrap < 4.4);
      for (const DressedLepton& lep : leptons)
        ifilter_discard(jets, deltaRLess(lep, 0.5));

      const size_t njets = jets.size();
      if (!njets)  vetoEvent;
      const int njets_gap = centralJetVeto(jets);

      double jpt1 = jets[0].pT();
      double jeta1 = jets[0].eta();
      double mjj = 0., jpt2 = 0., dphijj = 0.;
      if (njets >= 2) {
        mjj = (jets[0].momentum() + jets[1].momentum()).mass();
        jpt2 = jets[1].pT();
        dphijj = deltaPhi(jets[0], jets[1]);
      }

      // MET
      Vector3 met_vec = apply<MissingMomentum>(event, "MET").vectorMPT();
      double met = met_vec.mod();

      // Cut on deltaPhi between MET and first 4 jets, but only if jet pT > 30 GeV
      bool dphi_fail = false;
      for (size_t i = 0; i < jets.size() && i < 4; ++i) {
        dphi_fail |= (deltaPhi(jets[i], met_vec) < 0.4 && jets[i].pT() > 30*GeV);
      }

      const bool pass_met_dphi = met > 200*GeV && !dphi_fail;
      const bool pass_vbf = pass_met_dphi && mjj > 200*GeV && jpt1 > 80*GeV && jpt2 > 50*GeV && njets >= 2 && !njets_gap;
      const bool pass_mono = pass_met_dphi && jpt1 > 120*GeV && fabs(jeta1) < 2.4;
      if (pass_mono)  _h["met_mono"].fill(met, weight);
      if (pass_vbf) {
        _h["met_vbf"].fill(met/GeV, weight);
        _h["mjj_vbf"].fill(mjj/GeV, weight);
        _h["dphijj_vbf"].fill(dphijj, weight);
      }
    }


    /// Normalise, scale and otherwise manipulate histograms here
    void finalize() {
      const double sf(crossSection() / femtobarn / sumOfWeights());
      for (map<string, HistoHandler>::iterator hit = _h.begin(); hit != _h.end(); ++hit) {
        scale(hit->second.histo, sf);
        if (_mode < 2)  constructRmiss(hit->second);
      }
    }


    void constructRmiss(HistoHandler& handler) {
      // Load transfer function from reference data file
      const YODA::Scatter2D& rmiss = refData(handler.d, handler.x, handler.y);
      const YODA::Scatter2D& numer = refData(handler.d, handler.x, handler.y + 1);
      const YODA::Scatter2D& denom = refData(handler.d, handler.x, handler.y + 2);
      for (size_t i = 0; i < handler.scatter->numPoints(); ++i) {
        const Point2D& r = rmiss.point(i); // SM Rmiss
        const Point2D& n = numer.point(i); // SM numerator
        const Point2D& d = denom.point(i); // SM denominator
        const HistoBin1D& b = handler.histo->bin(i); // BSM
        double bsmy;
        try {
          bsmy = b.height();
        } catch (const Exception&) { // LowStatsError or WeightError
          bsmy = 0;
        }
        double bsmey;
        try {
          bsmey = b.heightErr();
        } catch (const Exception&) { // LowStatsError or WeightError
          bsmey = 0;
        }
        // Combined numerator
        double sm_plus_bsm = n.y() + bsmy;
        // Rmiss central value
        double rmiss_y = safediv(sm_plus_bsm, d.y());
        // Ratio error (Rmiss = SM_num/SM_denom + BSM/SM_denom ~ Rmiss_SM + BSM/SM_denom
        double rmiss_p = sqrt(r.yErrPlus()*r.yErrPlus()   + safediv(bsmey*bsmey, d.y()*d.y()));
        double rmiss_m = sqrt(r.yErrMinus()*r.yErrMinus() + safediv(bsmey*bsmey, d.y()*d.y()));
        // Set new values
        Point2D& p = handler.scatter->point(i); // (SM + BSM) Rmiss
        p.setY(rmiss_y);
        p.setYErrMinus(rmiss_m);
        p.setYErrPlus(rmiss_p);
      }
    }


  protected:

    // Analysis-mode switch
    size_t _mode;

    /// Histograms
    map<string, HistoHandler> _h;

  };



  /// ATLAS pTmiss+jets specialisation for Znunu channel
  class ATLAS_2017_I1609448_Znunu : public ATLAS_2017_I1609448 {
  public:
    ATLAS_2017_I1609448_Znunu()
      : ATLAS_2017_I1609448("ATLAS_2017_I1609448_Znunu")
    {
      _mode = 0;
    }
  };

  /// ATLAS pTmiss+jets specialisation for Zmumu channel
  class ATLAS_2017_I1609448_Zmumu : public ATLAS_2017_I1609448 {
  public:
    ATLAS_2017_I1609448_Zmumu()
      : ATLAS_2017_I1609448("ATLAS_2017_I1609448_Zmumu")
    {
      _mode = 1;
    }
  };

  /// ATLAS pTmiss+jets specialisation for Zee channel
  class ATLAS_2017_I1609448_Zee : public ATLAS_2017_I1609448 {
  public:
    ATLAS_2017_I1609448_Zee()
      : ATLAS_2017_I1609448("ATLAS_2017_I1609448_Zee")
    {
      _mode = 2;
    }
  };


  // Hooks for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1609448);
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1609448_Znunu);
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1609448_Zmumu);
  DECLARE_RIVET_PLUGIN(ATLAS_2017_I1609448_Zee);

}
