// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  /// $Wt$ at 8 TeV
  class ATLAS_2015_I1397635 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2015_I1397635);


    void init() {

      // Eta ranges
      Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT >= 1.0*MeV;
      Cut eta_lep  = Cuts::abseta < 2.5;

      // All final state particles
      FinalState fs(eta_full);

      // Get photons to dress leptons
      IdentifiedFinalState photons(fs);
      photons.acceptIdPair(PID::PHOTON);

      // Projection to find the electrons
      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState electrons(el_id);
      electrons.acceptTauDecays(true);
      declare(electrons, "electrons");
      DressedLeptons dressedelectrons(photons, electrons, 0.1, eta_lep && Cuts::pT > 25*GeV, true);
      declare(dressedelectrons, "dressedelectrons");
      DressedLeptons ewdressedelectrons(photons, electrons, 0.1, eta_full, true);

      // Projection to find the muons
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(true);
      declare(muons, "muons");
      DressedLeptons dressedmuons(photons, muons, 0.1, eta_lep && Cuts::pT > 25*GeV, true);
      declare(dressedmuons, "dressedmuons");
      DressedLeptons ewdressedmuons(photons, muons, 0.1, eta_full, true);

      // Projection to find neutrinos and produce MET
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);
      declare(neutrinos, "neutrinos");


      // Jet clustering.
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(ewdressedelectrons);
      vfs.addVetoOnThisFinalState(ewdressedmuons);
      vfs.addVetoOnThisFinalState(neutrinos);
      FastJets jets(vfs,FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      declare(jets, "jets");

      _histo = bookHisto1D(1,1,1);
    }


    void analyze(const Event& event) {

      // Get the selected objects, using the projections.
      vector<DressedLepton> electrons = apply<DressedLeptons>(event, "dressedelectrons").dressedLeptons();
      vector<DressedLepton> muons = apply<DressedLeptons>(event, "dressedmuons").dressedLeptons();
      // also make basic event selection cuts for leptons
      if (electrons.empty() && muons.empty())  vetoEvent;
      if (electrons.size() + muons.size() != 2) vetoEvent;

      // next selection cuts for jets
      const Jets jets = apply<FastJets>(event, "jets").jets(Cuts::pT>20*GeV && Cuts::abseta < 2.5, cmpMomByPt);
      if (jets.size() != 1) vetoEvent;

      // and selection cuts for b-tagging
      Jets bjets;
      // Check overlap of jets/leptons.
      foreach (Jet jet, jets) {
        // if dR(el,jet) < 0.4 skip the event
        foreach (DressedLepton el, electrons) {
          if (deltaR(jet, el) < 0.4)  vetoEvent;
        }
        // if dR(mu,jet) < 0.4 skip the event
        foreach (DressedLepton mu, muons) {
          if (deltaR(jet, mu) < 0.4)  vetoEvent;
        }
        // Count the number of b-tags
        // We have to check that the ghost-matched B hadrons have pT > 5 GeV
        // By default jet.bTags() returns all B hadrons without cuts
        bool btagged = jet.bTags(Cuts::pT >= 5*GeV).size();
        if (btagged)  bjets += jet;
      }
      if (bjets.size() != 1) vetoEvent;

      // Now evaluate MET selection
      // Get the neutrinos from the event record (they have pT > 0.0 and |eta| < 4.5 at this stage
      const Particles& neutrinos = apply<PromptFinalState>(event, "neutrinos").particlesByPt();
      FourMomentum met;
      foreach (const Particle& nu, neutrinos)  met += nu.momentum();
      if (met.pT() <= 20*GeV)  vetoEvent;

      // Make the plot
      _histo->fill(1, event.weight());
    }


    // Normalise histograms etc., after the run
    void finalize() {
      scale(_histo, crossSection() / femtobarn / sumOfWeights());
    }


  private:

    Histo1DPtr _histo;

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1397635);

}
