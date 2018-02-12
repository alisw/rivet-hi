#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/InvMassFinalState.hh"
#include "Rivet/Projections/MissingMomentum.hh"

namespace Rivet {


  class CMS_2016_I1473674 : public Analysis {
  public:

    // Minimal constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1473674);


    // Set up projections and book histograms
    void init() {

      // Complete final state
      FinalState fs;

      // Parton level top quarks
      declare(PartonicTops(PartonicTops::E_MU, false), "LeptonicPartonTops");
      declare(PartonicTops(PartonicTops::HADRONIC),    "HadronicPartonTops");

      // Projections for dressed electrons and muons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);
      //
      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      addProjection(electrons, "Electrons");
      DressedLeptons dressed_electrons(photons, electrons, 0.1);
      addProjection(dressed_electrons, "DressedElectrons");
      //
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      addProjection(muons, "Muons");
      DressedLeptons dressed_muons(photons, muons, 0.1);
      addProjection(dressed_muons, "DressedMuons");

      // Projection for jets
      VetoedFinalState fs_jets(FinalState(-MAXDOUBLE, MAXDOUBLE, 0*GeV));
      fs_jets.addVetoOnThisFinalState(dressed_muons);
      addProjection(FastJets(fs_jets, FastJets::ANTIKT, 0.5), "Jets");

      // Projections for MET
      addProjection(MissingMomentum(), "MET");


      // Booking of histograms
      _hist_met = bookHisto1D(5, 1, 1);
      _hist_ht  = bookHisto1D(6, 1, 1);
      _hist_st  = bookHisto1D(7, 1, 1);
      _hist_wpt = bookHisto1D(8, 1, 1);
    }


    /// Per-event analysis
    void analyze(const Event& event) {
      const double weight = event.weight();

      // Select ttbar -> lepton+jets at parton level, removing tau decays
      const Particles leptonicpartontops = apply<ParticleFinder>(event, "LeptonicPartonTops").particlesByPt();
      if (leptonicpartontops.size() != 1) vetoEvent;
      const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt();
      if (hadronicpartontops.size() != 1) vetoEvent;

      // Select ttbar -> lepton+jets at particle level
      const DressedLeptons& dressed_electrons = applyProjection<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& dressed_muons = applyProjection<DressedLeptons>(event, "DressedMuons");
      if (dressed_electrons.dressedLeptons().size() + dressed_muons.dressedLeptons().size() != 1) vetoEvent;
      const FourMomentum lepton = (dressed_electrons.dressedLeptons().empty() ? dressed_muons : dressed_electrons).dressedLeptons()[0];

      // MET
      const MissingMomentum& met = applyProjection<MissingMomentum>(event, "MET");
      _hist_met->fill(met.visibleMomentum().pT()/GeV, weight);

      // HT and ST
      const FastJets& jetpro = applyProjection<FastJets>(event, "Jets");
      const Jets jets = jetpro.jetsByPt(20*GeV);

      double ht = 0.0;
      for (const Jet& j : jets) {
        if (deltaR(j.momentum(), lepton) > 0.3) {
          ht += j.pT();
        }
      }

      double st = ht + lepton.pT() + met.visibleMomentum().pT();
      _hist_ht->fill(ht/GeV, weight);
      _hist_st->fill(st/GeV, weight);

      // WPT
      const FourMomentum w = lepton - met.visibleMomentum();
      _hist_wpt->fill(w.pT()/GeV, weight);
    }


    /// Normalize histograms
    void finalize() {
      normalize(_hist_met);
      normalize(_hist_ht);
      normalize(_hist_st);
      normalize(_hist_wpt);
    }

  private:
    Histo1DPtr _hist_met, _hist_ht, _hist_st, _hist_wpt;

  };


  DECLARE_RIVET_PLUGIN(CMS_2016_I1473674);

}
