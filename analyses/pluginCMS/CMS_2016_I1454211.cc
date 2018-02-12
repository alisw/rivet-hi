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


  /// Boosted ttbar in pp collisions at sqrtS = 8 TeV
  /// @todo Use persistent weight counters
  class CMS_2016_I1454211 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_I1454211);


    // Set up projections and book histograms
    void init() {

      // Complete final state
      FinalState fs;

      // Partonic tops
      declare(PartonicTops(PartonicTops::ELECTRON, false), "ElectronPartonTops");
      declare(PartonicTops(PartonicTops::MUON, false),     "MuonPartonTops");
      declare(PartonicTops(PartonicTops::HADRONIC),        "HadronicPartonTops");

      // Projection for electrons and muons
      IdentifiedFinalState photons(fs, PID::PHOTON);

      const Cut leptonCuts = Cuts::pt > 45*GeV && Cuts::abseta < 2.1;

      IdentifiedFinalState el_id(fs, {{PID::ELECTRON, -PID::ELECTRON}});
      PromptFinalState electrons(el_id);
      DressedLeptons dressed_electrons(photons, electrons, 0.1, leptonCuts);
      declare(dressed_electrons, "DressedElectrons");

      IdentifiedFinalState mu_id(fs, {{PID::MUON, -PID::MUON}});
      PromptFinalState muons(mu_id);
      DressedLeptons dressed_muons(photons, muons, 0.1, leptonCuts);
      declare(dressed_muons, "DressedMuons");

      // Projection for jets
      VetoedFinalState fs_jets(fs);
      fs_jets.addVetoOnThisFinalState(dressed_muons);
      fs_jets.addVetoOnThisFinalState(dressed_electrons);
      fs_jets.vetoNeutrinos();
      declare(FastJets(fs_jets, FastJets::ANTIKT, 0.5), "ak5jets");
      declare(FastJets(fs_jets, FastJets::CAM, 0.8), "ca8jets");

      _hEl_topPt_parton          = bookHisto1D("d01-x01-y01"); // dsigma/dpt(top quark), el ch
      _hEl_topPt_particle        = bookHisto1D("d02-x01-y01"); // dsigma/dpt(top jet), el ch
      _hEl_topY_parton           = bookHisto1D("d03-x01-y01"); // dsigma/dy(top quark), el ch
      _hEl_topY_particle         = bookHisto1D("d04-x01-y01"); // dsigma/dy(top jet), el ch
      _hMu_topPt_parton          = bookHisto1D("d05-x01-y01"); // dsigma/dpt(top quark), mu ch
      _hMu_topPt_particle        = bookHisto1D("d06-x01-y01"); // dsigma/dpt(top jet), mu ch
      _hMu_topY_parton           = bookHisto1D("d07-x01-y01"); // dsigma/dy(top quark), mu ch
      _hMu_topY_particle         = bookHisto1D("d08-x01-y01"); // dsigma/dy(top jet), mu ch
      _hComb_topPt_parton        = bookHisto1D("d09-x01-y01"); // dsigma/dpt(top quark), comb ch
      _hComb_topPt_particle      = bookHisto1D("d10-x01-y01"); // dsigma/dpt(top jet), comb ch
      _hComb_topY_parton         = bookHisto1D("d11-x01-y01"); // dsigma/dy(top quark), comb ch
      _hComb_topY_particle       = bookHisto1D("d12-x01-y01"); // dsigma/dy(top jet), comb ch

      _hEl_topPt_parton_norm     = bookHisto1D("d13-x01-y01"); // 1/sigma dsigma/dpt(top quark), el ch
      _hEl_topPt_particle_norm   = bookHisto1D("d14-x01-y01"); // 1/sigma dsigma/dpt(top jet), el ch
      _hEl_topY_parton_norm      = bookHisto1D("d15-x01-y01"); // 1/sigma dsigma/dy(top quark), el ch
      _hEl_topY_particle_norm    = bookHisto1D("d16-x01-y01"); // 1/sigma dsigma/dy(top jet), el ch
      _hMu_topPt_parton_norm     = bookHisto1D("d17-x01-y01"); // 1/sigma dsigma/dpt(top quark), mu ch
      _hMu_topPt_particle_norm   = bookHisto1D("d18-x01-y01"); // 1/sigma dsigma/dpt(top jet), mu ch
      _hMu_topY_parton_norm      = bookHisto1D("d19-x01-y01"); // 1/sigma dsigma/dy(top quark), mu ch
      _hMu_topY_particle_norm    = bookHisto1D("d20-x01-y01"); // 1/sigma dsigma/dy(top jet), mu ch
      _hComb_topPt_parton_norm   = bookHisto1D("d21-x01-y01"); // 1/sigma dsigma/dpt(top quark), comb ch
      _hComb_topPt_particle_norm = bookHisto1D("d22-x01-y01"); // 1/sigma dsigma/dpt(top jet), comb ch
      _hComb_topY_parton_norm    = bookHisto1D("d23-x01-y01"); // 1/sigma dsigma/dy(top quark), comb ch
      _hComb_topY_particle_norm  = bookHisto1D("d24-x01-y01"); // 1/sigma dsigma/dy(top jet), comb ch

      _hMu_cutflow = bookHisto1D("mu_cutflow", 7, -0.5, 6.5);
      _hEl_cutflow = bookHisto1D("el_cutflow", 7, -0.5, 6.5);
    }


    // per event analysis
    void analyze(const Event& event) {

       // Total-events cutflow entries
      _hMu_cutflow->fill(0.);
      _hEl_cutflow->fill(0.);

      // Do parton-level selection and channel determination
      int partonCh = 0; //0 non-semi-lep, 1 muon, 2 electron
      const Particles muonpartontops = apply<ParticleFinder>(event, "MuonPartonTops").particlesByPt();
      const Particles electronpartontops = apply<ParticleFinder>(event, "ElectronPartonTops").particlesByPt();
      if (electronpartontops.size() == 0 && muonpartontops.size() == 1) partonCh = 1;
      else if (electronpartontops.size() == 1 && muonpartontops.size() == 0) partonCh = 2;
      else vetoEvent;
      const Particles hadronicpartontops = apply<ParticleFinder>(event, "HadronicPartonTops").particlesByPt();
      if (hadronicpartontops.size() != 1) vetoEvent;

      if (partonCh == 1) _hMu_cutflow->fill(1.); // muon at parton level
      if (partonCh == 2) _hEl_cutflow->fill(1.); // electron at parton level

      // Get hadronic parton-level top
      const FourMomentum& partonTopP4 = hadronicpartontops.front();

      // Do particle-level selection and channel determination
      const DressedLeptons& dressed_electrons = apply<DressedLeptons>(event, "DressedElectrons");
      const DressedLeptons& dressed_muons = apply<DressedLeptons>(event, "DressedMuons");

      bool passParticleLep = false, passParticleTop = false;
      FourMomentum lepton, particleTopP4;
      if (partonCh == 1 && dressed_muons.dressedLeptons().size() == 1 && dressed_electrons.dressedLeptons().size() == 0) {
        passParticleLep = true;
        _hMu_cutflow->fill(3.); //muon at particle level
        lepton = dressed_muons.dressedLeptons()[0].momentum();
      }
      if (partonCh == 2 && dressed_muons.dressedLeptons().size() == 0 && dressed_electrons.dressedLeptons().size() == 1) {
        passParticleLep = true;
        _hEl_cutflow->fill(3.); //electron at particle level
        lepton = dressed_electrons.dressedLeptons()[0].momentum();
      }

      if (passParticleLep) {

        // Jet cuts
        Cut jetCuts = Cuts::pt > 30*GeV && Cuts::abseta < 2.4;
        Jets genBjets, genTjets;
        int nGenBjets = 0, nGenTjets = 0;

        const FastJets& AK5jets = apply<FastJets>(event, "ak5jets");
        for (const Jet& jet : AK5jets.jetsByPt(jetCuts)) {
          if (deltaR(jet, lepton) > M_PI / 2.0) continue;
          if (deltaR(jet, lepton) < 0.1) continue;
          genBjets.push_back(jet);
          nGenBjets += 1;
        }

        const FastJets& CA8jets = apply<FastJets>(event, "ca8jets");
        for (const Jet& jet : CA8jets.jetsByPt(jetCuts)) {
          if (deltaR(jet, lepton) < M_PI / 2.0) continue;
          if (jet.mass() < 140*GeV) continue;
          if (jet.mass() > 250*GeV) continue;
          genTjets.push_back(jet);
          nGenTjets += 1;
        }

        if (nGenBjets >=1) {
          if (partonCh == 1) _hMu_cutflow->fill(4.); // muon at parton level
          if (partonCh == 2) _hEl_cutflow->fill(4.); // electron at parton level
          if (nGenTjets >= 1) {
            passParticleTop = true;
            if (partonCh == 1) _hMu_cutflow->fill(5.); // muon at parton level
            if (partonCh == 2) _hEl_cutflow->fill(5.); // electron at parton level
            particleTopP4 = genTjets[0];
          }
        }
      }

      const double weight = event.weight();
      if (partonCh == 1) {
        _nMu += weight;
        _hMu_topPt_parton->fill(partonTopP4.pT()/GeV, weight);
        _hMu_topPt_parton_norm->fill(partonTopP4.pT()/GeV, weight);
        _hComb_topPt_parton->fill(partonTopP4.pT()/GeV, weight);
        _hComb_topPt_parton_norm->fill(partonTopP4.pT()/GeV, weight);

        if (partonTopP4.pT() >= 400*GeV) {
          _nPassParton_mu += weight;
          _hMu_cutflow->fill(2.);
          _hMu_topY_parton->fill(partonTopP4.rapidity(), weight);
          _hMu_topY_parton_norm->fill(partonTopP4.rapidity(), weight);
          _hComb_topY_parton->fill(partonTopP4.rapidity(), weight);
          _hComb_topY_parton_norm->fill(partonTopP4.rapidity(), weight);
        }

        if (passParticleTop) {
          _hMu_topPt_particle->fill(particleTopP4.pT()/GeV, weight);
          _hMu_topPt_particle_norm->fill(particleTopP4.pT()/GeV, weight);
          _hComb_topPt_particle->fill(particleTopP4.pT()/GeV, weight);
          _hComb_topPt_particle_norm->fill(particleTopP4.pT()/GeV, weight);

          if (particleTopP4.pT() >= 400*GeV) {
            _nPassParticle_mu += weight;
            _hMu_cutflow->fill(6.);
            _hMu_topY_particle->fill(particleTopP4.rapidity(), weight);
            _hMu_topY_particle_norm->fill(particleTopP4.rapidity(), weight);
            _hComb_topY_particle->fill(particleTopP4.rapidity(), weight);
            _hComb_topY_particle_norm->fill(particleTopP4.rapidity(), weight);
          }
        }
      }

      if (partonCh == 2){
        _nEl += weight;
        _hEl_topPt_parton->fill(partonTopP4.pT()/GeV, weight);
        _hEl_topPt_parton_norm->fill(partonTopP4.pT()/GeV, weight);
        _hComb_topPt_parton->fill(partonTopP4.pT()/GeV, weight);
        _hComb_topPt_parton_norm->fill(partonTopP4.pT()/GeV, weight);

        if (partonTopP4.pT() >= 400*GeV) {
          _nPassParton_el += weight;
          _hEl_cutflow->fill(2.);
          _hEl_topY_parton->fill(partonTopP4.rapidity(), weight);
          _hEl_topY_parton_norm->fill(partonTopP4.rapidity(), weight);
          _hComb_topY_parton->fill(partonTopP4.rapidity(), weight);
          _hComb_topY_parton_norm->fill(partonTopP4.rapidity(), weight);
        }

        if (passParticleTop) {
          _hEl_topPt_particle->fill(particleTopP4.pT()/GeV, weight);
          _hEl_topPt_particle_norm->fill(particleTopP4.pT()/GeV, weight);
          _hComb_topPt_particle->fill(particleTopP4.pT()/GeV, weight);
          _hComb_topPt_particle_norm->fill(particleTopP4.pT()/GeV, weight);

          if (particleTopP4.pT() >= 400*GeV) {
            _nPassParticle_el += weight;
            _hEl_cutflow->fill(6.);
            _hEl_topY_particle->fill(particleTopP4.rapidity(), weight);
            _hEl_topY_particle_norm->fill(particleTopP4.rapidity(), weight);
            _hComb_topY_particle->fill(particleTopP4.rapidity(), weight);
            _hComb_topY_particle_norm->fill(particleTopP4.rapidity(), weight);
          }
        }
      }
    }


    void finalize() {

      normalize({_hMu_topPt_parton_norm, _hMu_topY_parton_norm, _hEl_topPt_parton_norm,
            _hEl_topY_parton_norm, _hComb_topPt_parton_norm, _hComb_topY_parton_norm}, 1.0, false);
      normalize({_hMu_topPt_particle_norm, _hMu_topY_particle_norm, _hEl_topPt_particle_norm,
            _hEl_topY_particle_norm, _hComb_topPt_particle_norm, _hComb_topY_particle_norm}, 1.0, false);

      scale(_hMu_topPt_particle, crossSection()/femtobarn / sumOfWeights());
      scale(_hEl_topPt_particle, crossSection()/femtobarn / sumOfWeights());
      scale(_hMu_topY_particle, crossSection()/femtobarn / sumOfWeights());
      scale(_hEl_topY_particle, crossSection()/femtobarn / sumOfWeights());
      scale(_hComb_topPt_particle, crossSection()/femtobarn / sumOfWeights());
      scale(_hComb_topY_particle, crossSection()/femtobarn / sumOfWeights());

      scale(_hMu_topPt_parton, crossSection()/femtobarn / sumOfWeights());
      scale(_hEl_topPt_parton, crossSection()/femtobarn / sumOfWeights());
      scale(_hMu_topY_parton, crossSection()/femtobarn / sumOfWeights());
      scale(_hEl_topY_parton, crossSection()/femtobarn / sumOfWeights());
      scale(_hComb_topPt_parton, crossSection()/femtobarn / sumOfWeights());
      scale(_hComb_topY_parton, crossSection()/femtobarn / sumOfWeights());

    }


  private:

    Histo1DPtr _hMu_topPt_parton, _hMu_topY_parton, _hEl_topPt_parton, _hEl_topY_parton, _hComb_topPt_parton, _hComb_topY_parton;
    Histo1DPtr _hMu_topPt_particle, _hMu_topY_particle, _hEl_topPt_particle, _hEl_topY_particle, _hComb_topPt_particle, _hComb_topY_particle;
    Histo1DPtr _hMu_topPt_parton_norm, _hMu_topY_parton_norm, _hEl_topPt_parton_norm, _hEl_topY_parton_norm, _hComb_topPt_parton_norm, _hComb_topY_parton_norm;
    Histo1DPtr _hMu_topPt_particle_norm, _hMu_topY_particle_norm, _hEl_topPt_particle_norm, _hEl_topY_particle_norm, _hComb_topPt_particle_norm, _hComb_topY_particle_norm;
    Histo1DPtr _hMu_cutflow, _hEl_cutflow;

    double _nMu = 0., _nEl = 0.;
    double _nPassParton_mu = 0.,_nPassParton_el = 0.;
    double _nPassParticle_mu = 0., _nPassParticle_el = 0.;

  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_I1454211);

}
