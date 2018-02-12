// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

namespace Rivet {


  /// Leading jet mass for boosted top quarks at 8 TeV
  class CMS_2017_I1518399 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2017_I1518399);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Dressed leptons
      IdentifiedFinalState photons(PID::PHOTON);
      ChargedLeptons charged_leptons;
      PromptFinalState prompt_leptons(charged_leptons);
      Cut leptonCuts = Cuts::pT > 45*GeV && Cuts::abseta < 2.1;
      DressedLeptons dressed_leptons(photons, prompt_leptons, 0.1, leptonCuts);
      declare(dressed_leptons, "DressedLeptons");

      // Jets
      VetoedFinalState fs_jets;
      fs_jets.vetoNeutrinos();
      declare(FastJets(fs_jets, FastJets::CAM, 1.2), "JetsCA12");

      // Partonic top for decay channel definition
      declare(PartonicTops(PartonicTops::E_MU, false), "LeptonicTops");
      declare(PartonicTops(PartonicTops::HADRONIC), "HadronicTops");

      // Main histograms
      _hist_mass        = bookHisto1D("d01-x01-y01");
      _hist_mass_norm   = bookHisto1D("d02-x01-y01");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Decay mode check
      const Particles& leptonicTops = apply<PartonicTops>(event, "LeptonicTops").particlesByPt();
      const Particles& hadronicTops = apply<PartonicTops>(event, "HadronicTops").particlesByPt();
      if (leptonicTops.size() != 1 || hadronicTops.size() != 1) vetoEvent;

      // Get the leptons
      const DressedLeptons& dressed_leptons = apply<DressedLeptons>(event, "DressedLeptons");

      // Leading dressed lepton
      const vector<DressedLepton> leptons = dressed_leptons.dressedLeptons();
      if (leptons.empty()) vetoEvent;
      Particle lepton;
      for (const Particle& l : leptons) {
        if (l.pT() > lepton.pT()) lepton = l;
      }

      // Get the jets
      const Jets& psjetsCA12 = applyProjection<FastJets>(event, "JetsCA12").jetsByPt(Cuts::pT > 50*GeV);

      // Subtract the lepton four vector from a jet in case of overlap and clean jets
      Jets cleanedJets;
      for (Jet jet : psjetsCA12) {
        if (deltaR(jet, lepton) < 1.2 )
          jet = Jet(jet.momentum()-lepton.momentum(), jet.particles(), jet.tags());
        if (jet.abseta() < 2.5) cleanedJets.push_back(jet);
      }
      std::sort(cleanedJets.begin(), cleanedJets.end(), cmpMomByPt);

      // Jet pT cuts
      if (cleanedJets.size() < 2) vetoEvent;
      if (cleanedJets.at(0).pT() < 400*GeV) vetoEvent;
      if (cleanedJets.at(1).pT() < 150*GeV) vetoEvent;

      // Jet veto
      if (cleanedJets.size() > 2 && cleanedJets.at(2).pT() > 150*GeV) vetoEvent;

      // Small distance between 2nd jet and lepton
      if (deltaR(cleanedJets.at(1), lepton) > 1.2) vetoEvent;

      // m(jet1) > m(jet2 +lepton)
      FourMomentum secondJetLepton = cleanedJets.at(1).momentum() + lepton.momentum();
      if (cleanedJets.at(0).mass() < secondJetLepton.mass()) vetoEvent;

      // Fill histograms
      const double weight = event.weight();
      _hist_mass->fill(cleanedJets.at(0).mass(), weight);
      _hist_mass_norm->fill(cleanedJets.at(0).mass(), weight);
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double sf = crossSection() * 1000 / sumOfWeights();
      scale(_hist_mass, sf);
      normalize(_hist_mass_norm, 1.0, false);
    }

    //@}


  private:

    // Histograms
    Histo1DPtr _hist_mass, _hist_mass_norm;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2017_I1518399);


}
