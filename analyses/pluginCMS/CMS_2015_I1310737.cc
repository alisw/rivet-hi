// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/VisibleFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ZFinder.hh"

namespace Rivet {


  class CMS_2015_I1310737 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2015_I1310737);


    /// Book histograms and initialise projections before the run
    void init() {

      FinalState fs; ///< @todo No cuts?
      VisibleFinalState visfs(fs);
      // Prompt leptons only
      PromptFinalState pfs(fs);

      ZFinder zeeFinder(pfs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::ELECTRON, 71.0*GeV, 111.0*GeV);
      declare(zeeFinder, "ZeeFinder");

      ZFinder zmumuFinder(pfs, Cuts::abseta < 2.4 && Cuts::pT > 20*GeV, PID::MUON, 71.0*GeV, 111.0*GeV);
      declare(zmumuFinder, "ZmumuFinder");

      VetoedFinalState jetConstits(visfs);
      jetConstits.addVetoOnThisFinalState(zeeFinder);
      jetConstits.addVetoOnThisFinalState(zmumuFinder);

      FastJets akt05Jets(jetConstits, FastJets::ANTIKT, 0.5);
      declare(akt05Jets, "AntiKt05Jets");


      _h_excmult_jets_tot = bookHisto1D(1, 1, 1);
      _h_incmult_jets_tot = bookHisto1D(2, 1, 1);
      _h_leading_jet_pt_tot = bookHisto1D(3, 1, 1);
      _h_second_jet_pt_tot = bookHisto1D(4, 1, 1);
      _h_third_jet_pt_tot = bookHisto1D(5, 1, 1);
      _h_fourth_jet_pt_tot = bookHisto1D(6, 1, 1);
      _h_leading_jet_eta_tot = bookHisto1D(7, 1, 1);
      _h_second_jet_eta_tot = bookHisto1D(8, 1, 1);
      _h_third_jet_eta_tot = bookHisto1D(9, 1, 1);
      _h_fourth_jet_eta_tot = bookHisto1D(10, 1, 1);
      _h_ht1_tot = bookHisto1D(11, 1, 1);
      _h_ht2_tot = bookHisto1D(12, 1, 1);
      _h_ht3_tot = bookHisto1D(13, 1, 1);
      _h_ht4_tot = bookHisto1D(14, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {;

      const ZFinder& zeeFS = apply<ZFinder>(event, "ZeeFinder");
      const ZFinder& zmumuFS = apply<ZFinder>(event, "ZmumuFinder");

      const Particles& zees = zeeFS.bosons();
      const Particles& zmumus = zmumuFS.bosons();

      // We did not find exactly one Z. No good.
      if (zees.size() + zmumus.size() != 1) {
        MSG_DEBUG("Did not find exactly one good Z candidate");
        vetoEvent;
      }

      // Find the (dressed!) leptons
      const Particles& dressedLeptons = zees.size() ? zeeFS.constituents() : zmumuFS.constituents();

      // Cluster jets
      // NB. Veto has already been applied on leptons and photons used for dressing
      const FastJets& fj = apply<FastJets>(event, "AntiKt05Jets");
      const Jets& jets = fj.jetsByPt(Cuts::abseta < 2.4 && Cuts::pT > 30*GeV);

      // Perform lepton-jet overlap and HT calculation
      double ht = 0;
      Jets goodjets;
      foreach (const Jet& j, jets) {
        // Decide if this jet is "good", i.e. isolated from the leptons
        /// @todo Nice use-case for any() and a C++11 lambda
        bool overlap = false;
        foreach (const Particle& l, dressedLeptons) {
          if (Rivet::deltaR(j, l) < 0.5) {
            overlap = true;
            break;
          }
        }

        // Fill HT and good-jets collection
        if (overlap) continue;
        goodjets.push_back(j);
        ht += j.pT();
      }

      // We don't care about events with no isolated jets
      if (goodjets.empty()) {
        MSG_DEBUG("No jets in event");
        vetoEvent;
      }


      /////////////////


      // Weight to be used for histo filling
      const double w = 0.5 * event.weight();

      // Fill jet number integral histograms
      _h_excmult_jets_tot->fill(goodjets.size(), w);
      /// @todo Could be better computed by toIntegral transform on exclusive histo
      for (size_t iJet = 1; iJet <= goodjets.size(); iJet++ )
        _h_incmult_jets_tot->fill(iJet, w);

      // Fill leading jet histograms
      const Jet& j1 = goodjets[0];
      _h_leading_jet_pt_tot->fill(j1.pT()/GeV, w);
      _h_leading_jet_eta_tot->fill(j1.abseta(), w);
      _h_ht1_tot->fill(ht/GeV, w);

      // Fill 2nd jet histograms
      if (goodjets.size() < 2) return;
      const Jet& j2 = goodjets[1];
      _h_second_jet_pt_tot->fill(j2.pT()/GeV, w);
      _h_second_jet_eta_tot->fill(j2.abseta(), w);
      _h_ht2_tot->fill(ht/GeV, w);

      // Fill 3rd jet histograms
      if (goodjets.size() < 3) return;
      const Jet& j3 = goodjets[2];
      _h_third_jet_pt_tot->fill(j3.pT()/GeV, w);
      _h_third_jet_eta_tot->fill(j3.abseta(), w);
      _h_ht3_tot->fill(ht/GeV, w);

      // Fill 4th jet histograms
      if (goodjets.size() < 4) return;
      const Jet& j4 = goodjets[3];
      _h_fourth_jet_pt_tot->fill(j4.pT()/GeV, w);
      _h_fourth_jet_eta_tot->fill(j4.abseta(), w);
      _h_ht4_tot->fill(ht/GeV, w);
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      const double norm = (sumOfWeights() != 0) ? crossSection()/sumOfWeights() : 1.0;

      MSG_INFO("Cross section = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << crossSection() << " pb");
      MSG_INFO("# Events      = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << numEvents() );
      MSG_INFO("SumW          = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(3) << sumOfWeights());
      MSG_INFO("Norm factor   = " << std::setfill(' ') << std::setw(14) << std::fixed << std::setprecision(6) << norm);

      scale(_h_excmult_jets_tot, norm );
      scale(_h_incmult_jets_tot, norm );
      scale(_h_leading_jet_pt_tot, norm );
      scale(_h_second_jet_pt_tot, norm );
      scale(_h_third_jet_pt_tot, norm );
      scale(_h_fourth_jet_pt_tot, norm );
      scale(_h_leading_jet_eta_tot, norm );
      scale(_h_second_jet_eta_tot, norm );
      scale(_h_third_jet_eta_tot, norm );
      scale(_h_fourth_jet_eta_tot, norm );
      scale(_h_ht1_tot, norm );
      scale(_h_ht2_tot, norm );
      scale(_h_ht3_tot, norm );
      scale(_h_ht4_tot, norm );
    }


  private:

    /// @name Histograms

    Histo1DPtr _h_excmult_jets_tot,  _h_incmult_jets_tot;
    Histo1DPtr _h_leading_jet_pt_tot, _h_second_jet_pt_tot, _h_third_jet_pt_tot, _h_fourth_jet_pt_tot;
    Histo1DPtr _h_leading_jet_eta_tot, _h_second_jet_eta_tot, _h_third_jet_eta_tot, _h_fourth_jet_eta_tot;
    Histo1DPtr _h_ht1_tot, _h_ht2_tot, _h_ht3_tot, _h_ht4_tot;

  };


  DECLARE_RIVET_PLUGIN(CMS_2015_I1310737);


}
