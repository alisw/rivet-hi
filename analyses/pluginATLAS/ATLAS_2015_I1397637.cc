// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/UnstableFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {


  class ATLAS_2015_I1397637 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2015_I1397637);


    /// Book projections and histograms
    void init() {

      // Base final state definition
      const FinalState fs(Cuts::abseta < 4.5);

      // Neutrinos for MET
      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);
      declare(neutrinos, "neutrinos");

      // Get photons used to dress leptons
      IdentifiedFinalState photons(fs);
      photons.acceptId(PID::PHOTON);

      // Use all bare muons as input to the DressedMuons projection
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState bare_mu(mu_id);
      bare_mu.acceptTauDecays(true);

      // Use all bare electrons as input to the DressedElectrons projection
      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState bare_el(el_id);
      bare_el.acceptTauDecays(true);

      // Use all bare leptons including taus for single-lepton filter
      IdentifiedFinalState lep_id(fs);
      lep_id.acceptIdPair(PID::MUON);
      lep_id.acceptIdPair(PID::ELECTRON);
      PromptFinalState bare_lep(lep_id);
      declare(bare_lep, "bare_lep");

      // Tau finding
      /// @todo Use TauFinder
      UnstableFinalState ufs;
      IdentifiedFinalState tau_id(ufs);
      tau_id.acceptIdPair(PID::TAU);
      PromptFinalState bare_tau(tau_id);
      declare(bare_tau, "bare_tau");

      // Muons and electrons must have |eta| < 2.5
      Cut eta_ranges = Cuts::abseta < 2.5;

      // Get dressed muons and the good muons (pt>25GeV)
      DressedLeptons all_dressed_mu(photons, bare_mu, 0.1, eta_ranges, true);
      DressedLeptons dressed_mu(photons, bare_mu, 0.1, eta_ranges && Cuts::pT > 25*GeV, true);
      declare(dressed_mu, "muons");

      // Get dressed electrons and the good electrons (pt>25GeV)
      DressedLeptons all_dressed_el(photons, bare_el, 0.1, eta_ranges, true);
      DressedLeptons dressed_el(photons, bare_el, 0.1, eta_ranges && Cuts::pT > 25*GeV, true);
      declare(dressed_el, "electrons");

      // Jet clustering
      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(all_dressed_el);
      vfs.addVetoOnThisFinalState(all_dressed_mu);
      vfs.addVetoOnThisFinalState(neutrinos);

      // Small-R jets
      /// @todo Use extra constructor args
      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(JetAlg::ALL_INVISIBLES);
      jets.useMuons(JetAlg::DECAY_MUONS);
      declare(jets, "jets");

      // Large-R jets
      /// @todo Use extra constructor args
      FastJets large_jets(vfs, FastJets::ANTIKT, 1.0);
      large_jets.useInvisibles(JetAlg::ALL_INVISIBLES);
      large_jets.useMuons(JetAlg::DECAY_MUONS);
      declare(large_jets, "fat_jets");


      /// Book histogram
      _h_pttop = bookHisto1D(1, 1, 1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Single lepton filter on bare leptons with no cuts
      const Particles& bare_lep = apply<PromptFinalState>(event, "bare_lep").particles();
      const Particles& bare_tau = apply<PromptFinalState>(event, "bare_tau").particles();
      if (bare_lep.size() + bare_tau.size() != 1) vetoEvent;

      // Electrons and muons
      const vector<DressedLepton>& electrons = apply<DressedLeptons>(event, "electrons").dressedLeptons();
      const vector<DressedLepton>& muons = apply<DressedLeptons>(event, "muons").dressedLeptons();
      if (electrons.size() + muons.size() != 1) vetoEvent;
      const DressedLepton& lepton = muons.empty() ? electrons[0] : muons[0];

      // Get the neutrinos from the event record (they have pT > 0.0 and |eta| < 4.5 at this stage
      const Particles& neutrinos = apply<PromptFinalState>(event, "neutrinos").particlesByPt();
      FourMomentum met;
      for (const Particle& nu : neutrinos) met += nu.momentum();
      if (met.pT() < 20*GeV) vetoEvent;


      // Thin jets and trimmed fat jets
      /// @todo Use Rivet built-in FJ trimming support
      const Jets& jets  = apply<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      const PseudoJets& fat_pjets = apply<FastJets>(event, "fat_jets").pseudoJetsByPt();
      const double Rfilt = 0.3, ptFrac_min = 0.05; ///< @todo Need to be careful about the units for the pT cut passed to FJ?
      PseudoJets trimmed_fat_pjets;
      fastjet::Filter trimmer(fastjet::JetDefinition(fastjet::kt_algorithm, Rfilt), fastjet::SelectorPtFractionMin(ptFrac_min));
      for (const PseudoJet& pjet : fat_pjets) trimmed_fat_pjets += trimmer(pjet);
      trimmed_fat_pjets = fastjet::sorted_by_pt(trimmed_fat_pjets);

      // Jet reclustering
      // Use a kT cluster sequence to recluster the trimmed jets so that a d12 can then be obtained from the reclustered jet
      vector<double> splittingScales;
      for (const PseudoJet& tpjet : trimmed_fat_pjets) {
        const PseudoJets tpjet_constits = tpjet.constituents();
        const fastjet::ClusterSequence kt_cs(tpjet_constits, fastjet::JetDefinition(fastjet::kt_algorithm, 1.5, fastjet::E_scheme, fastjet::Best));
        const PseudoJets kt_jets = kt_cs.inclusive_jets();
        const double d12 = 1.5 * sqrt(kt_jets[0].exclusive_subdmerge(1));
        splittingScales += d12;
      }
      Jets trimmed_fat_jets;
      for (size_t i = 0; i < trimmed_fat_pjets.size(); ++i) {
        const Jet tj = trimmed_fat_pjets[i];
        if (tj.mass()          <= 100*GeV) continue;
        if (tj.pT()            <= 300*GeV) continue;
        if (splittingScales[i] <=  40*GeV) continue;
        if (tj.abseta()        >=     2.0) continue;
        trimmed_fat_jets += tj;
      }
      if (trimmed_fat_jets.empty()) vetoEvent;

      // Jet b-tagging
      Jets bjets, non_bjets;
      for (const Jet& jet : jets)
        (jet.bTagged() ? bjets : non_bjets) += jet;
      if (bjets.empty()) vetoEvent;


      // Boosted selection: lepton/jet overlap
      const double transmass = sqrt( 2 * lepton.pT() * met.pT() * (1 - cos(deltaPhi(lepton, met))) );
      if (transmass + met.pt() <= 60*GeV) vetoEvent;
      int lepJetIndex = -1;
      for (size_t i = 0; i < jets.size(); ++i) {
        const Jet& jet = jets[i];
        if (deltaR(jet, lepton) < 1.5) {
          lepJetIndex = i;
          break;
        }
      }
      if (lepJetIndex < 0) vetoEvent;
      const Jet& ljet = jets[lepJetIndex];

      // Boosted selection: lepton-jet/fat-jet matching
      int fatJetIndex = -1;
      for (size_t j = 0; j < trimmed_fat_jets.size(); ++j) {
        const Jet& fjet = trimmed_fat_jets[j];
        const double dR_fatjet = deltaR(ljet, fjet);
        const double dPhi_fatjet = deltaPhi(lepton, fjet);
        if (dR_fatjet > 1.5 && dPhi_fatjet > 2.3) {
          fatJetIndex = j;
          break;
        }
      }
      if (fatJetIndex < 0) vetoEvent;
      const Jet& fjet = trimmed_fat_jets[fatJetIndex];

      // Boosted selection: b-tag matching
      const bool lepbtag = ljet.bTagged();
      bool hadbtag = false;
      for (const Jet& bjet : bjets) {
        hadbtag |= (deltaR(fjet, bjet) < 1.0);
      }

      // Fill histo if selection passed
      if (hadbtag || lepbtag) _h_pttop->fill(fjet.pT()/GeV, event.weight());
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      scale(_h_pttop, crossSection()/femtobarn / sumOfWeights());
    }


  private:

    Histo1DPtr _h_pttop;

  };


  DECLARE_RIVET_PLUGIN(ATLAS_2015_I1397637);

}
