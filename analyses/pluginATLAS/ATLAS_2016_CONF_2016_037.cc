// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/Sphericity.hh"
#include "Rivet/Projections/SmearedParticles.hh"
#include "Rivet/Projections/SmearedJets.hh"
#include "Rivet/Projections/SmearedMET.hh"
#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  /// @brief ATLAS 2016 2 -SS-lepton / 3-lepton SUSY search, from 13.2/fb CONF note
  class ATLAS_2016_CONF_2016_037 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_CONF_2016_037);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState calofs(Cuts::abseta < 4.9);
      declare(calofs, "Clusters");
      FastJets fj(calofs, FastJets::ANTIKT, 0.4);
      declare(fj, "TruthJets");
      declare(SmearedJets(fj, JET_SMEAR_ATLAS_RUN2, [](const Jet& j) {
            if (j.abseta() > 2.5) return 0.;
            return j.bTagged(Cuts::pT > 5*GeV) ? 0.70 :
              j.cTagged(Cuts::pT > 5*GeV) ? 1/12. :
              j.tauTagged(Cuts::pT > 5*GeV) ? 1/54. : 1/380.; }), "Jets");

      MissingMomentum mm(calofs);
      declare(mm, "TruthMET");
      declare(SmearedMET(mm, MET_SMEAR_ATLAS_RUN2), "MET");

      FinalState es(Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.47 && !Cuts::absetaIn(1.37, 1.52) && Cuts::pT > 10*GeV);
      declare(es, "TruthElectrons");
      declare(SmearedParticles(es, ELECTRON_EFF_ATLAS_RUN2, ELECTRON_SMEAR_ATLAS_RUN2), "Electrons");

      FinalState mus(Cuts::abspid == PID::MUON && Cuts::abseta < 2.5 && Cuts::pT > 10*GeV);
      declare(mus, "TruthMuons");
      declare(SmearedParticles(mus, MUON_EFF_ATLAS_RUN2, MUON_SMEAR_ATLAS_RUN2), "Muons");

      ChargedFinalState cfs(Cuts::abseta < 2.5);
      declare(cfs, "TruthTracks");
      declare(SmearedParticles(cfs, TRK_EFF_ATLAS_RUN2), "Tracks");


      // Book histograms/counters
      _h_3l1 = bookCounter("SR3l1");
      _h_3l2 = bookCounter("SR3l2");
      _h_0b1 = bookCounter("SR0b1");
      _h_0b2 = bookCounter("SR0b2");
      _h_1b = bookCounter("SR1b");
      _h_3b = bookCounter("SR3b");
      _h_1bDD = bookCounter("SR1bDD");
      _h_3bDD = bookCounter("SR3bDD");
      _h_1bGG = bookCounter("SR1bGG");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get baseline electrons, muons, and jets
      Particles elecs = apply<ParticleFinder>(event, "Electrons").particlesByPt();
      Particles muons = apply<ParticleFinder>(event, "Muons").particlesByPt();
      Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.8);
      const Jets bjets = filter_select(jets, [&](const Jet& j) { return j.bTagged(Cuts::pT > 5*GeV); });


      // Jet/electron/muon overlap removal and selection
      // Remove any electron or muon within dR = 0.2 of a b-tagged jet
      for (const Jet& bj : bjets) {
        ifilter_discard(elecs, deltaRLess(bj, 0.2, RAPIDITY));
        ifilter_discard(muons, deltaRLess(bj, 0.2, RAPIDITY));
      }
      // Remove any untagged jet within dR = 0.2 of an electron or muon
      for (const Particle& e : elecs)
        ifilter_discard(jets, deltaRLess(e, 0.2, RAPIDITY));
      for (const Particle& m : muons)
        ifilter_discard(jets, deltaRLess(m, 0.2, RAPIDITY));
      // Remove any untagged low-multiplicity/muon-dominated jet within dR = 0.4 of a muon
      for (const Particle& m : muons)
        ifilter_discard(jets, [&](const Jet& j) {
            if (deltaR(m, j, RAPIDITY) > 0.4) return false;
            const Particles trks = j.particles(Cuts::abscharge != 0);
            if (trks.size() < 3) return true;
            return m.pT()/j.pT() > 0.5 && m.pT()/sum(trks, pT, 0.0) > 0.7;
          });
      // Remove any electron or muon near a remaining jet, with a shrinking cone
      const auto lcone_iso_fn = [&](const Particle& l) {
        const double dr = min(0.4, 0.04 + 10*GeV/l.pT());
        return any(jets, deltaRLess(l, dr, RAPIDITY));
      };
      ifilter_discard(elecs, lcone_iso_fn);
      ifilter_discard(muons, lcone_iso_fn);
      // Track-sharing e,mu also filtered, but that decision can't be made here
      const Jets& sigjets = jets;
      const Jets& sigbjets = bjets;


      // Lepton isolation
      Particles sigelecs = filter_select(elecs, Cuts::abseta < 2);
      Particles sigmuons = muons;
      ifilter_select(sigelecs, ParticleEffFilter(ELECTRON_IDEFF_ATLAS_RUN2_MEDIUM));
      const Particles trks = apply<ParticleFinder>(event, "Tracks").particles();
      const Particles clus = apply<ParticleFinder>(event, "Clusters").particles();
      ifilter_discard(sigelecs, [&](const Particle& e) {
          const double R = min(0.2, 10*GeV/e.pT());
          double ptsum = -e.pT(), etsum = -e.Et();
          for (const Particle& t : trks)
            if (deltaR(t,e) < R) ptsum += t.pT();
          for (const Particle& c : clus)
            if (deltaR(c,e) < 0.2) etsum += c.pT(); ///< @todo Bit vague about "energy"
          return ptsum / e.pT() > 0.06 || etsum / e.pT() > 0.06;
        });
      ifilter_discard(sigmuons, [&](const Particle& m) {
          const double R = min(0.3, 10*GeV/m.pT());
          double ptsum = -m.pT();
          for (const Particle& t : trks)
            if (deltaR(t,m) < R) ptsum += t.pT();
          return ptsum / m.pT() > 0.06;
        });
      /// @todo Note is vague about whether "signal lepton" defn includes pT > 20?
      ifilter_discard(sigelecs, Cuts::pT > 20*GeV);
      ifilter_discard(sigmuons, Cuts::pT > 20*GeV);


      // MET calculation (NB. done generically, with smearing, rather than via explicit physics objects)
      const Vector3 vmet = -apply<SmearedMET>(event, "MET").vectorEt();
      const double etmiss = vmet.mod();


      //////////////////


      // Event selection cuts
      const Particles sigleptons = sigelecs + sigmuons;
      if (sigleptons.size() < 2) vetoEvent;
      if (sigleptons.size() == 2 && sigleptons[0].charge() != sigleptons[1].charge()) vetoEvent;

      // Jet sub-selections and meff calculation
      const Jets sigjets25 = filter_select(sigjets, Cuts::pT > 25*GeV);
      const Jets sigjets40 = filter_select(sigjets25, Cuts::pT > 40*GeV);
      const Jets sigjets50 = filter_select(sigjets40, Cuts::pT > 50*GeV);
      /// @todo Is meff specific to the jet pT cut?
      const double meff = sum(sigjets, pT, 0.0) + sum(sigleptons, pT, 0.0);

      // Fill counters
      const double w = event.weight();
      if (sigleptons.size() >= 3 && sigbjets.empty() && sigjets40.size() >= 4 && etmiss > 150*GeV) _h_3l1->fill(w);
      if (sigleptons.size() >= 3 && sigbjets.empty() && sigjets40.size() >= 4 && etmiss > 200*GeV && meff > 1500*GeV) _h_3l2->fill(w);
      if (sigleptons.size() >= 2 && sigbjets.empty() && sigjets25.size() >= 6 && etmiss > 150*GeV && meff > 500*GeV) _h_0b1->fill(w);
      if (sigleptons.size() >= 2 && sigbjets.empty() && sigjets40.size() >= 6 && etmiss > 150*GeV && meff > 900*GeV) _h_0b2->fill(w);
      if (sigleptons.size() >= 2 && sigbjets.size() >= 1 && sigjets25.size() >= 6 && etmiss > 200*GeV && meff > 650*GeV) _h_1b->fill(w);
      if (sigleptons.size() >= 2 && sigbjets.size() >= 3 && sigjets25.size() >= 6 && etmiss > 150*GeV && meff > 600*GeV) _h_3b->fill(w);
      if (filter_select(sigleptons, Cuts::charge < 0).size() >= 2) {
        if (sigleptons.size() >= 2 && sigbjets.size() >= 1 && sigjets50.size() >= 6 && meff > 1200*GeV) _h_1bDD->fill(w);
        if (sigleptons.size() >= 2 && sigbjets.size() >= 3 && sigjets50.size() >= 6 && meff > 1000*GeV) _h_3bDD->fill(w);
        if (sigleptons.size() >= 2 && sigbjets.size() >= 1 && sigjets50.size() >= 6 && meff > 1800*GeV) _h_1bGG->fill(w);
      }

    }


    /// Normalise counters after the run
    void finalize() {

      const double sf = 13.2*crossSection()/femtobarn/sumOfWeights();
      scale({_h_3l1, _h_3l2, _h_0b1, _h_0b2, _h_1b, _h_3b, _h_1bDD, _h_3bDD, _h_1bGG}, sf);

    }

    //@}


  private:

    /// @name Histograms
    //@{
    CounterPtr _h_3l1, _h_3l2, _h_0b1, _h_0b2, _h_1b, _h_3b, _h_1bDD, _h_3bDD, _h_1bGG;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_CONF_2016_037);


}
