// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  // ATLAS $tt+b(b)$ at 8~TeV
  class ATLAS_2015_I1390114 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2015_I1390114);


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

      // Jet clustering.
      VetoedFinalState vfs;
      vfs.addVetoOnThisFinalState(ewdressedelectrons);
      vfs.addVetoOnThisFinalState(ewdressedmuons);
      vfs.addVetoOnThisFinalState(neutrinos);
      FastJets jets(vfs,FastJets::ANTIKT, 0.4);
      jets.useInvisibles();
      declare(jets, "jets");

      _histo = bookHisto1D(1,1,1);
      _ratio = bookScatter2D(2,1,1, true);
      _aux   = bookHisto1D("_aux", 1, 0.5, 1.5);
    }


    void analyze(const Event& event) {

      const double weight = event.weight();

      // Get the selected objects, using the projections.
      vector<DressedLepton> electrons = apply<DressedLeptons>(event, "dressedelectrons").dressedLeptons();
      vector<DressedLepton> muons = apply<DressedLeptons>(event, "dressedmuons").dressedLeptons();
      if (electrons.empty() && muons.empty())  vetoEvent;

      Jets jets = apply<FastJets>(event, "jets").jets((Cuts::pT>20*GeV) && (Cuts::abseta < 2.5), cmpMomByPt);
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

      // Evaluate basic event selection
      bool pass_1lep = (electrons.size() == 1 && muons.empty()) || (muons.size() == 1 && electrons.empty());
      if (pass_1lep && bjets.size() >= 3 && jets.size() >= 5)  _histo->fill(1, weight);

      if (muons.size() == 1 && electrons.size() == 1 && bjets.size() >= 3) {
        if (muons[0].charge() * electrons[0].charge() < 0.0)  _histo->fill(2, weight);
      }

      DressedLepton *lep1 = NULL, *lep2 = NULL;
      bool zveto = false;
      if (electrons.size() == 2 && muons.empty()) {
        lep1 = &electrons[0];
        lep2 = &electrons[1];
      }
      else if (muons.size() == 2 && electrons.empty()) {
        lep1 = &muons[0];
        lep2 = &muons[1];
      }
      else if (electrons.size() == 1 && muons.size() == 1) {
        lep1 = &electrons[0];
        lep2 = &muons[0];
        zveto = true;
      }
      if (lep1 && lep2 && lep1->charge() * lep2->charge() < 0.0 && bjets.size() >= 2) {
        double mass = (lep1->momentum() + lep2->momentum()).mass();
        bool pass_2lep = mass > 15*GeV && (zveto || !(mass > 81*GeV && mass < 101*GeV));
        if ( pass_2lep && bjets.size() >= 4) {
          _histo->fill(3, weight);
          _histo->fill(4, weight);
        }
        if ( pass_2lep && jets.size() >= 4) {
          _aux->fill(1, weight);
        }
      }
    }


    void finalize() {
      const double sf(crossSection()/sumOfWeights()/femtobarn);
      scale(_histo, sf);
      scale(_aux,  sf);

      // construct ratio
      const double  n = _histo->bin(3).sumW();
      const double dN = _histo->bin(3).sumW2();
      const double  d = _aux->bin(0).sumW();
      const double dD = _aux->bin(0).sumW2();
      const double  r = safediv(n, d);
      double e = sqrt( safediv(r * (1 - r), d) );
      if ( _aux->effNumEntries() != _aux->numEntries() ) {
        // use F. James's approximation for weighted events:
        e = sqrt( safediv((1 - 2 * r) * dN + r * r * dD, d * d) );
      }
      _ratio->point(0).setY(100.0 * r, 100.0 * e); // convert into percentage
    }


  private:

    Histo1DPtr _histo, _aux;
    Scatter2DPtr _ratio;

  };


  // Hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1390114);

}
