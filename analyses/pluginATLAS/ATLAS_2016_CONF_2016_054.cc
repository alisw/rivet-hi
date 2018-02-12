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


  /// @brief ATLAS 2016 1-lepton SUSY search, from 14.8/fb CONF note
  class ATLAS_2016_CONF_2016_054 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_CONF_2016_054);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState calofs(Cuts::abseta < 4.9);
      FastJets fj(calofs, FastJets::ANTIKT, 0.4);
      declare(fj, "TruthJets");
      declare(SmearedJets(fj, JET_SMEAR_ATLAS_RUN2, [](const Jet& j) {
            if (j.abseta() > 2.5) return 0.;
            return j.bTagged(Cuts::pT > 5*GeV) ? 0.77 :
              j.cTagged(Cuts::pT > 5*GeV) ? 1/6.2 : 1/134.; }), "Jets");

      MissingMomentum mm(calofs);
      declare(mm, "TruthMET");
      declare(SmearedMET(mm, MET_SMEAR_ATLAS_RUN2), "MET");

      FinalState es(Cuts::abseta < 2.47 && Cuts::pT > 7*GeV && Cuts::abspid == PID::ELECTRON);
      declare(es, "TruthElectrons");
      declare(SmearedParticles(es, ELECTRON_EFF_ATLAS_RUN2, ELECTRON_SMEAR_ATLAS_RUN2), "Electrons");

      FinalState mus(Cuts::abseta < 2.5 && Cuts::pT > 6*GeV && Cuts::abspid == PID::MUON);
      declare(mus, "TruthMuons");
      declare(SmearedParticles(mus, MUON_EFF_ATLAS_RUN2, MUON_SMEAR_ATLAS_RUN2), "Muons");


      // Book histograms/counters
      _h_gg2j = bookCounter("GG-2j");
      _h_gg6j0 = bookCounter("GG-6j-0bulk");
      _h_gg6j1 = bookCounter("GG-6j-1highmass");
      _h_gg4j0 = bookCounter("GG-4j-0lowx");
      _h_gg4j1 = bookCounter("GG-4j-1lowxbveto");
      _h_gg4j2 = bookCounter("GG-4j-2highx");
      _h_ss4j0 = bookCounter("SS-4j-0x12");
      _h_ss4j1 = bookCounter("SS-4j-1lowx");
      _h_ss5j0 = bookCounter("SS-5j-0x12");
      _h_ss5j1 = bookCounter("SS-5j-1highx");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get baseline electrons, muons, and jets
      Particles elecs = apply<ParticleFinder>(event, "Electrons").particles();
      Particles muons = apply<ParticleFinder>(event, "Muons").particles();
      Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 4.5);

      // Jet/electron/muons overlap removal and selection
      // Remove any jet within dR = 0.2 of an electron
      for (const Particle& e : elecs)
        ifilter_discard(jets, deltaRLess(e, 0.2, RAPIDITY));
      // Remove any electron within dR = 0.01 of a muon
      for (const Particle& m : muons)
        ifilter_discard(elecs, deltaRLess(m, 0.01, RAPIDITY));
      // Assemble b-jets collection, and remove muons within dR = 0.2 of a b-tagged jet
      Jets bjets;
      for (const Jet& j : jets) {
        if (j.abseta() < 2.5 && j.pT() > 30*GeV && j.bTagged(Cuts::pT > 5*GeV)) {
          bjets += j;
          ifilter_discard(muons, deltaRLess(j, 0.2, RAPIDITY));
        }
      }
      // Remove any jet within dR = 0.2 of a muon if track conditions are met
      for (const Particle& m : muons)
        ifilter_discard(jets, [&](const Jet& j){
            if (deltaR(j,m) > 0.2) return false;
            /// @todo Add track efficiency random filtering
            const Particles trks = j.particles(Cuts::abscharge > 0 && Cuts::pT > 0.5*GeV);
            return trks.size() < 4 || m.pT()/j.pT() > 0.7;
          });
      // Remove any muon within dR = 0.2 of a remaining jet if the same track conditions are *not* met
      /// @todo There must be nicer way to do complementary removal...
      for (const Jet& j : jets) {
        /// @todo Add track efficiency random filtering
        const size_t ntrks = j.particles(Cuts::abscharge > 0 && Cuts::pT > 0.5*GeV).size();
        ifilter_discard(muons, [&](const Particle& m){
            if (deltaR(j,m) > 0.2) return false;
            return ntrks > 3 && m.pT()/j.pT() < 0.7;
          });
      }
      // Remove any muon with dR close to a remaining jet, via a functional form
      for (const Jet& j : jets)
        ifilter_discard(muons, [&](const Particle& m) { return deltaR(m,j, RAPIDITY) < min(0.4, 0.04 + 10*GeV/m.pT()); });


      // Signal jet selection
      const Jets sigjets = filter_select(jets, Cuts::pT > 30*GeV && Cuts::abseta < 2.8);
      const Jets sigbjets = bjets;

      // "Gradient-loose" signal lepton selection
      const ParticleEffFilter grad_loose_filter([](const Particle& e) { return e.pT() > 60*GeV ? 0.98 : 0.95; });
      Particles sigelecs = filter_select(elecs, grad_loose_filter);
      Particles sigmuons = filter_select(muons, grad_loose_filter);
      // Tight electron selection (NB. assuming independent eff to gradient-loose... hmm)
      ifilter_select(sigelecs, ParticleEffFilter(ELECTRON_IDEFF_ATLAS_RUN2_TIGHT));


      // MET calculation (NB. done generically, with smearing, rather than via explicit physics objects)
      const Vector3 vmet = -apply<SmearedMET>(event, "MET").vectorEt();
      const double etmiss = vmet.mod();


      //////////////////


      // Event selection cuts
      if (sigelecs.size() + sigmuons.size() != 1) vetoEvent;
      const Particle siglepton = sigelecs.empty() ? sigmuons.front() : sigelecs.front();

      // mT and m_eff
      const double mT = sqrt(2*siglepton.pT()*etmiss*(1-cos(deltaPhi(siglepton,vmet))));
      const double meff = siglepton.pT() + sum(sigjets, pT, 0.0) + etmiss;

      // Aplanarities
      Sphericity sph;
      vector<FourMomentum> vecs;
      transform(sigjets, vecs, mom);
      sph.calc(vecs);
      const double jet_aplanarity = sph.aplanarity();
      vecs += siglepton.mom();
      sph.calc(vecs);
      const double lepton_aplanarity = sph.aplanarity();


      //////////////////


      // Fill counters
      const double w = event.weight();
      // GG
      if (siglepton.pT() < 35*GeV && sigjets.size() >= 2 &&
          sigjets[0].pT() > 200*GeV && sigjets[1].pT() > 30*GeV &&
          mT > 100*GeV && etmiss > 460*GeV && etmiss/meff > 0.35) _h_gg2j->fill(w);
      if (siglepton.pT() > 35*GeV && sigjets.size() >= 6 &&
          sigjets[0].pT() > 125*GeV && sigjets[5].pT() > 30*GeV &&
          mT > 225*GeV && etmiss > 250*GeV && meff > 1000*GeV && etmiss/meff > 0.2 &&
          jet_aplanarity > 0.04) _h_gg6j0->fill(w);
      if (siglepton.pT() > 35*GeV && sigjets.size() >= 6 &&
          sigjets[0].pT() > 125*GeV && sigjets[5].pT() > 30*GeV &&
          mT > 225*GeV && etmiss > 250*GeV && meff > 2000*GeV && etmiss/meff > 0.1 &&
          jet_aplanarity > 0.04) _h_gg6j1->fill(w);
      if (sigjets.size() >= 4 && sigjets[3].pT() > 100*GeV &&
          mT > 125*GeV && etmiss > 250*GeV && meff > 2000*GeV && jet_aplanarity > 0.06) _h_gg4j0->fill(w);
      if (sigjets.size() >= 4 && sigjets[3].pT() > 100*GeV && sigbjets.empty() &&
          mT > 125*GeV && etmiss > 250*GeV && meff > 2000*GeV && jet_aplanarity > 0.03) _h_gg4j1->fill(w);
      if (siglepton.pT() > 35*GeV &&
          sigjets.size() >= 4 && sigjets[0].pT() > 400*GeV && inRange(sigjets[3].pT(), 30*GeV, 100*GeV) &&
          mT > 475*GeV && etmiss > 250*GeV && meff > 1600*GeV && etmiss/meff > 0.3) _h_gg4j2->fill(w);
      // SS
      if (siglepton.pT() > 35*GeV && sigjets.size() >= 4 && sigjets[3].pT() > 50*GeV &&
          mT > 175*GeV && etmiss > 300*GeV && meff > 1200*GeV && lepton_aplanarity > 0.08) _h_ss4j0->fill(w);
      if (siglepton.pT() > 35*GeV && sigjets.size() >= 5 && sigjets[4].pT() > 50*GeV && sigbjets.empty() &&
          mT > 175*GeV && etmiss > 300*GeV && etmiss/meff > 0.2) _h_ss5j0->fill(w);
      if (siglepton.pT() > 35*GeV && sigjets.size() >= 4 && sigjets[0].pT() > 250*GeV && sigjets[3].pT() > 30*GeV &&
          inRange(mT, 150*GeV, 400*GeV) && etmiss > 250*GeV && lepton_aplanarity > 0.03) _h_ss4j1->fill(w);
      if (siglepton.pT() > 35*GeV && sigjets.size() >= 5 && sigjets[4].pT() > 30*GeV &&
          mT > 400*GeV && etmiss > 400*GeV && lepton_aplanarity > 0.03) _h_ss5j1->fill(w);

    }


    /// Normalise counters after the run
    void finalize() {

      const double sf = 14.8*crossSection()/femtobarn/sumOfWeights();
      scale({_h_gg2j, _h_gg6j0, _h_gg6j1, _h_gg4j0, _h_gg4j1, _h_gg4j2}, sf);
      scale({_h_ss4j0, _h_ss4j1, _h_ss5j0,_h_ss5j1}, sf);

    }

    //@}


  private:

    /// @name Histograms
    //@{
    CounterPtr _h_gg2j, _h_gg6j0, _h_gg6j1, _h_gg4j0, _h_gg4j1, _h_gg4j2;
    CounterPtr _h_ss4j0, _h_ss4j1, _h_ss5j0,_h_ss5j1;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_CONF_2016_054);


}
