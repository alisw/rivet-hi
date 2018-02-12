// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/WFinder.hh"

namespace Rivet {


  /// W + jets jet multiplicities and pT
  class ATLAS_2010_S8919674 : public Analysis {
  public:

    /// @name Constructors etc.
    //@{

    /// Constructor
    ATLAS_2010_S8919674()
      : Analysis("ATLAS_2010_S8919674")
    {    }

    //@}


  public:

    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Set up projections to find the electron and muon Ws
      FinalState fs;
      Cut cuts = (Cuts::abseta < 1.37 || Cuts::absetaIn(1.52, 2.47)) && Cuts::pT > 20*GeV;
      WFinder wfinder_e(fs, cuts, PID::ELECTRON, 0*GeV, 1000*GeV, 25*GeV);
      declare(wfinder_e, "W_e");
      WFinder wfinder_mu(fs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::MUON, 0*GeV, 1000*GeV, 25*GeV);
      declare(wfinder_mu, "W_mu");

      // Input for the jets: no neutrinos, no muons, and no electron which passed the electron cuts
      VetoedFinalState veto;
      veto.addVetoOnThisFinalState(wfinder_e);
      veto.addVetoOnThisFinalState(wfinder_mu);
      veto.addVetoPairId(PID::MUON);
      veto.vetoNeutrinos();
      FastJets jets(veto, FastJets::ANTIKT, 0.4);
      declare(jets, "jets");

      /// Book histograms
      _h_el_njet_inclusive = bookHisto1D(1,1,1);
      _h_mu_njet_inclusive = bookHisto1D(2,1,1);
      _h_el_pT_jet1 = bookHisto1D(5,1,1);
      _h_mu_pT_jet1 = bookHisto1D(6,1,1);
      _h_el_pT_jet2 = bookHisto1D(7,1,1);
      _h_mu_pT_jet2 = bookHisto1D(8,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      const Jets& jets = apply<FastJets>(event, "jets").jetsByPt(20.0*GeV);

      const WFinder& We = apply<WFinder>(event, "W_e");
      if (We.bosons().size() == 1) {
        const FourMomentum p_miss = We.constituentNeutrinos()[0];
        const FourMomentum p_lept = We.constituentLeptons()[0];
        if (p_miss.Et() > 25*GeV && We.mT() > 40*GeV) {
          Jets js;
          for (const Jet& j : jets) {
            if (j.abseta() < 2.8 && deltaR(p_lept, j.momentum()) > 0.5)
              js.push_back(j);
          }
          _h_el_njet_inclusive->fill(0, weight);
          if (js.size() >= 1) {
            _h_el_njet_inclusive->fill(1, weight);
            _h_el_pT_jet1->fill(js[0].pT(), weight);
          }
          if (js.size() >= 2) {
            _h_el_njet_inclusive->fill(2, weight);
            _h_el_pT_jet2->fill(js[1].pT(), weight);
          }
          if (js.size() >= 3) {
            _h_el_njet_inclusive->fill(3, weight);
          }
        }
      }

      const WFinder& Wm = apply<WFinder>(event, "W_mu");
      if (Wm.bosons().size() == 1) {
        const FourMomentum p_miss = Wm.constituentNeutrinos()[0];
        const FourMomentum p_lept = Wm.constituentLeptons()[0];
        if (p_miss.Et() > 25*GeV && Wm.mT() > 40*GeV) {
          Jets js;
          for (const Jet& j : jets) {
            if (j.abseta() < 2.8 && deltaR(p_lept, j.momentum()) > 0.5)
              js.push_back(j);
          }
          _h_mu_njet_inclusive->fill(0, weight);
          if (js.size() >= 1) {
            _h_mu_njet_inclusive->fill(1, weight);
            _h_mu_pT_jet1->fill(js[0].pT(), weight);
          }
          if (js.size() >= 2) {
            _h_mu_njet_inclusive->fill(2, weight);
            _h_mu_pT_jet2->fill(js[1].pT(), weight);
          }
          if (js.size() >= 3) {
            _h_mu_njet_inclusive->fill(3, weight);
          }
          if (js.size() >= 4) {
            _h_mu_njet_inclusive->fill(4, weight);
          }
        }
      }

    }


    /// Normalise histograms etc., after the run
    void finalize() {
      double normfac = crossSection()/sumOfWeights();
      scale(_h_el_njet_inclusive, normfac);
      scale(_h_mu_njet_inclusive, normfac);
      scale(_h_el_pT_jet1, normfac);
      scale(_h_mu_pT_jet1, normfac);
      scale(_h_el_pT_jet2, normfac);
      scale(_h_mu_pT_jet2, normfac);
    }

    //@}


  private:

    /// @name Histograms
    //@{

    Histo1DPtr _h_el_njet_inclusive;
    Histo1DPtr _h_mu_njet_inclusive;
    Histo1DPtr _h_el_pT_jet1;
    Histo1DPtr _h_mu_pT_jet1;
    Histo1DPtr _h_el_pT_jet2;
    Histo1DPtr _h_mu_pT_jet2;
    //@}

  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2010_S8919674);

}
