// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"

namespace Rivet {


  /// @brief Z/gamma cross section measurement at 8 TeV
  class ATLAS_2016_I1448301 : public Analysis {
  public:

    /// Constructor
    ATLAS_2016_I1448301(string name="ATLAS_2016_I1448301") : Analysis(name) {
      _mode = 0; // pick electron channel by default
      setNeedsCrossSection(true);
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Prompt photons
      const Cut photoncut = Cuts::abspid == PID::PHOTON && Cuts::pT > 15*GeV && Cuts::abseta < 2.37;
      PromptFinalState photon_fs(photoncut);
      declare(photon_fs, "Photons");

      // Prompt leptons
      const PromptFinalState barelepton_fs = _mode ? Cuts::abspid == PID::MUON : Cuts::abspid == PID::ELECTRON;

      // Dressed leptons
      const IdentifiedFinalState allphoton_fs(PID::PHOTON); // photons used for lepton dressing
      const Cut leptoncut = Cuts::pT > 25*GeV && Cuts::abseta < 2.47;
      const DressedLeptons dressedlepton_fs(allphoton_fs, barelepton_fs, 0.1, leptoncut, true); // use *all* photons for lepton dressing
      declare(dressedlepton_fs, "Leptons");

      // MET (prompt neutrinos)
      VetoedFinalState ivfs;
      ivfs.addVetoOnThisFinalState(VisibleFinalState());
      declare(PromptFinalState(ivfs), "MET");

      // Jets
      VetoedFinalState jet_fs;
      jet_fs.vetoNeutrinos();
      jet_fs.addVetoPairId(PID::MUON);
      const FastJets fastjets(jet_fs, FastJets::ANTIKT, 0.4);
      declare(fastjets, "Jets");


      // Histograms
      if (_mode == 2) {
        _h["vvg"]     = bookHisto1D( 2, 1, 1);
        _h["vvgg"]    = bookHisto1D( 4, 1, 1);
        _h["pT"]      = bookHisto1D( 7, 1, 1);
        _h["pT_0jet"] = bookHisto1D( 8, 1, 1);
      } else {
        const size_t ch = 1 + bool(_mode);
        _h["llg"]       = bookHisto1D( 1, 1, ch);
        _h["llg_comb"]  = bookHisto1D( 1, 1, 3);
        _h["llgg"]      = bookHisto1D( 3, 1, ch);
        _h["llgg_comb"] = bookHisto1D( 3, 1, 3);
        //
        _h["pT"]       = bookHisto1D( 5, 1, 1);
        _h["pT_0jet"]  = bookHisto1D( 6, 1, 1);
        _h["M"]        = bookHisto1D( 9, 1, 1);
        _h["M_0jet"]   = bookHisto1D(10, 1, 1);
        _h["Njets"]    = bookHisto1D(11, 1, 1);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const double weight = event.weight();

      // Get objects
      vector<DressedLepton> leptons = apply<DressedLeptons>(event, "Leptons").dressedLeptons();
      const Particles& photons = apply<PromptFinalState>(event, "Photons").particlesByPt();
      const Jets jets = apply<FastJets>(event, "Jets").jetsByPt();

      if (_mode == 2) {
        const FinalState& metfs = apply<PromptFinalState>(event, "MET");
        Vector3 met_vec;
        for (const Particle& p : metfs.particles()) met_vec += p.mom().perpVec();

        if (met_vec.mod() < 100*GeV) vetoEvent;
        if (photons.empty()) vetoEvent;

        if (photons.size() > 1) { // nu nu y y

          bool yy_veto = false;
          yy_veto |= photons[0].pT() < 22*GeV;
          yy_veto |= photons[1].pT() < 22*GeV;
          yy_veto |= met_vec.mod() < 110*GeV;
          const double yyPhi = (photons[0].momentum() + photons[1].momentum()).phi();
          yy_veto |= fabs(yyPhi - met_vec.phi()) < 2.62 || fabs(yyPhi - met_vec.phi()) > 3.66;
          yy_veto |= deltaR(photons[0], photons[1]) < 0.4;

          // Photon isolation calculated by jets, count jets
          Jet ph0_jet, ph1_jet;
          double min_dR_ph0_jet = 999., min_dR_ph1_jet = 999.;
          size_t njets = 0;
          for (const Jet& j : jets) {
            if (j.pT() > 30*GeV && j.abseta() < 4.5) {
              if (deltaR(j, photons[0]) > 0.3 && deltaR(j, photons[1]) > 0.3)  ++njets;
            }
            if (deltaR(j, photons[0]) < min_dR_ph0_jet) {
              min_dR_ph0_jet = deltaR(j, photons[0]);
              ph0_jet = j;
            }
            if (deltaR(j, photons[1]) < min_dR_ph1_jet) {
              min_dR_ph1_jet = deltaR(j, photons[1]);
              ph1_jet = j;
            }
          }
          double photon0iso = 0., photon1iso = 0.;
          if (min_dR_ph0_jet < 0.4)  photon0iso = ph0_jet.pT() - photons[0].pT();
          if (min_dR_ph1_jet < 0.4)  photon1iso = ph1_jet.pT() - photons[1].pT();
          yy_veto |= photon0iso/photons[0].pT() > 0.5;
          yy_veto |= photon1iso/photons[1].pT() > 0.5;

          if (!yy_veto) {
            _h["vvgg"]->fill(0.5, weight);
            if (!njets)  _h["vvgg"]->fill(1.5, weight);
          }
        } // end of nu nu y y section


        if (photons[0].pT() < 130*GeV)  vetoEvent;
        if (fabs(fabs(deltaPhi(photons[0], met_vec)) - 3.14) > 1.57)  vetoEvent;

        // Photon isolation calculated by jets, count jets
        Jet ph_jet;
        double min_dR_ph_jet = 999.;
        size_t njets = 0;
        for (const Jet& j : jets) {
          if (j.pT() > 30*GeV && j.abseta() < 4.5) {
            if (deltaR(j, photons[0]) > 0.3)  ++njets;
          }
          if (deltaR(j, photons[0]) < min_dR_ph_jet) {
            min_dR_ph_jet = deltaR(j, photons[0]);
            ph_jet = j;
          }
        }
        double photoniso = 0;
        if (min_dR_ph_jet < 0.4)  photoniso = ph_jet.pT() - photons[0].pT();
        if (photoniso/photons[0].pT() > 0.5)  vetoEvent;

        const double pTgamma = photons[0].pT()/GeV;
        _h["pT"]->fill(pTgamma, weight);
        _h["vvg"]->fill(0.5, weight);
        if (!njets) {
          _h["vvg"]->fill(1.5, weight);
          _h["pT_0jet"]->fill(pTgamma, weight);
        }
      } // end of nu nu y (y) section


      else {

        // Dilepton candidate
        if (leptons.size() < 2) vetoEvent;

        // Sort the dressed leptons by pt
        std::sort(leptons.begin(), leptons.end(), cmpMomByPt);

        vector<DressedLepton> lep_p, lep_m;
        for (const DressedLepton& lep : leptons) {
          if (lep.charge() > 0.)  lep_p.push_back(lep);
          if (lep.charge() < 0.)  lep_m.push_back(lep);
        }

        if (lep_p.empty() || lep_m.empty())  vetoEvent;
        if (lep_p[0].abspid() != lep_m[0].abspid())  vetoEvent;
        if ((lep_p[0].momentum() + lep_m[0].momentum()).mass() < 40*GeV)  vetoEvent;

        // Photon lepton overlap removal
        if (photons.empty())  vetoEvent;

        if (photons.size() > 1) {
          bool veto = false;
          veto |= deltaR(photons[0], lep_p[0]) < 0.4;
          veto |= deltaR(photons[0], lep_m[0]) < 0.4;
          veto |= deltaR(photons[1], lep_p[0]) < 0.4;
          veto |= deltaR(photons[1], lep_m[0]) < 0.4;
          veto |= deltaR(photons[0], photons[1]) < 0.4;

          Jet ph0_jet, ph1_jet;
          double min_dR_ph0_jet = 999., min_dR_ph1_jet=999.;
          int njets = 0;
          for (const Jet& j : jets){
            if (j.pT() > 30*GeV && j.abseta() < 4.5) {
              if (deltaR(j, lep_p[0]) > 0.3 && deltaR(j, lep_m[0]) > 0.3) {
                if (deltaR(j, photons[0]) > 0.3 && deltaR(j, photons[1]) > 0.3 )  ++njets;
              }
            }
            if (deltaR(j, photons[0]) < min_dR_ph0_jet) {
              min_dR_ph0_jet = deltaR(j, photons[0]);
              ph0_jet = j;
            }
            if (deltaR(j, photons[1]) < min_dR_ph1_jet) {
              min_dR_ph1_jet = deltaR(j, photons[1]);
              ph1_jet = j;
            }
          }
          double photon0iso = 0, photon1iso = 0;
          if (min_dR_ph0_jet < 0.4) photon0iso = ph0_jet.pT() - photons[0].pT();
          if (min_dR_ph1_jet < 0.4) photon1iso = ph1_jet.pT() - photons[1].pT();
          veto |= photon0iso/photons[0].pT() > 0.5;
          veto |= photon1iso/photons[1].pT() > 0.5;

          // Fill plots
          if (!veto) {
            _h["llgg"]->fill(0.5, weight);
            _h["llgg_comb"]->fill(0.5, weight);
            if (!njets) {
              _h["llgg"]->fill(1.5, weight);
              _h["llgg_comb"]->fill(1.5, weight);
            }
          }
        }

        if (deltaR(photons[0], lep_p[0]) < 0.7)  vetoEvent;
        if (deltaR(photons[0], lep_m[0]) < 0.7)  vetoEvent;

        // Photon isolation calculated by jets, count jets
        Jet ph_jet;
        double min_dR_ph_jet = 999.;
        size_t njets = 0;
        for (const Jet& j : jets) {
          if (j.pT() > 30*GeV && j.abseta() < 4.5) {
            if (deltaR(j, lep_p[0]) > 0.3 && deltaR(j, lep_m[0]) > 0.3 && deltaR(j, photons[0]) > 0.3)  ++njets;
          }
          if (deltaR(j, photons[0]) < min_dR_ph_jet) {
            min_dR_ph_jet = deltaR(j, photons[0]);
            ph_jet = j;
          }
        }
        double photoniso = 0;
        if (min_dR_ph_jet < 0.4)  photoniso = ph_jet.pT() - photons[0].pT();
        if (photoniso/photons[0].pT() > 0.5)  vetoEvent;


        // Fill plots
        const double pTgamma = photons[0].pT()/GeV;
        const double mllgamma = (lep_p[0].momentum() + lep_m[0].momentum() + photons[0].momentum()).mass()/GeV;

        _h["pT"]->fill(pTgamma,  weight);
        _h["M"]->fill(mllgamma, weight);
        _h["Njets"]->fill(njets < 3? njets : 3, weight);

        _h["llg"]->fill(0.5, weight);
        _h["llg_comb"]->fill(0.5, weight);

        if (!njets) {
          _h["pT_0jet"]->fill(pTgamma, weight);
          _h["M_0jet"]->fill(mllgamma, weight);
          _h["llg"]->fill(1.5, weight);
          _h["llg_comb"]->fill(1.5, weight);
        }
      } // end of _mode check
    } // end of analysis


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection()/femtobarn/sumOfWeights();
      for (const auto& kv : _h) scale(kv.second, sf);
    }

    //@}


  protected:

    // Data members like post-cuts event weight counters go here
    size_t _mode;

  private:

    /// Histograms
    map<string, Histo1DPtr> _h;

  };


  struct ATLAS_2016_I1448301_EL : public ATLAS_2016_I1448301 {
    ATLAS_2016_I1448301_EL() : ATLAS_2016_I1448301("ATLAS_2016_I1448301_EL") { _mode = 0; }
  };

  struct ATLAS_2016_I1448301_MU : public ATLAS_2016_I1448301 {
    ATLAS_2016_I1448301_MU() : ATLAS_2016_I1448301("ATLAS_2016_I1448301_MU") { _mode = 1; }
  };

  struct ATLAS_2016_I1448301_NU : public ATLAS_2016_I1448301 {
    ATLAS_2016_I1448301_NU() : ATLAS_2016_I1448301("ATLAS_2016_I1448301_NU") { _mode = 2; }
  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1448301);
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1448301_EL);
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1448301_MU);
  DECLARE_RIVET_PLUGIN(ATLAS_2016_I1448301_NU);

}
