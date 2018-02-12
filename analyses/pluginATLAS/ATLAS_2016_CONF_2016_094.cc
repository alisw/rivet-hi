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


  /// @brief ATLAS 2016 1-lepton + many jets SUSY search, from 14.8/fb CONF note
  class ATLAS_2016_CONF_2016_094 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(ATLAS_2016_CONF_2016_094);


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
            return j.bTagged(Cuts::pT > 5*GeV) ? 0.80 :
              j.cTagged(Cuts::pT > 5*GeV) ? 1/6. : 1/106.; }), "Jets");

      // MissingMomentum mm(calofs);
      // declare(mm, "TruthMET");
      // declare(SmearedMET(mm, MET_SMEAR_ATLAS_RUN2), "MET");

      FinalState es(Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.47 && !Cuts::absetaIn(1.37, 1.52) && Cuts::pT > 10*GeV);
      declare(es, "TruthElectrons");
      declare(SmearedParticles(es, ELECTRON_EFF_ATLAS_RUN2, ELECTRON_SMEAR_ATLAS_RUN2), "Electrons");

      FinalState mus(Cuts::abspid == PID::MUON && Cuts::abseta < 2.4 && Cuts::pT > 10*GeV);
      declare(mus, "TruthMuons");
      declare(SmearedParticles(mus, MUON_EFF_ATLAS_RUN2, MUON_SMEAR_ATLAS_RUN2), "Muons");


      // Book histograms/counters
      _h_08j40_0b = bookCounter("08j40_0b");
      _h_09j40_0b = bookCounter("09j40_0b");
      _h_10j40_0b = bookCounter("10j40_0b");
      _h_08j40_3b = bookCounter("08j40_3b");
      _h_09j40_3b = bookCounter("09j40_3b");
      _h_10j40_3b = bookCounter("10j40_3b");
      _h_08j60_0b = bookCounter("08j60_0b");
      _h_09j60_0b = bookCounter("09j60_0b");
      _h_10j60_0b = bookCounter("10j60_0b");
      _h_08j60_3b = bookCounter("08j60_3b");
      _h_09j60_3b = bookCounter("09j60_3b");
      _h_10j60_3b = bookCounter("10j60_3b");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get baseline electrons, muons, and jets
      // NB. for electrons, we don't apply the loose ID here, since we don't want to double-count effs with later use of tight ID
      Particles elecs = apply<ParticleFinder>(event, "Electrons").particles();
      Particles muons = apply<ParticleFinder>(event, "Muons").particles();
      Jets jets = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 20*GeV && Cuts::abseta < 2.4);
      ifilter_select(jets, JetEffFilter([](const Jet& j) { return j.pT() > 60*GeV ? 1.0 : 0.94; }));


      // Jet/electron/muon overlap removal and selection
      // Remove any untagged jet within dR = 0.2 of an electron
      for (const Particle& e : elecs)
        ifilter_discard(jets, [&](const Jet& j) { return !j.bTagged(Cuts::pT > 5*GeV) && deltaR(e, j, RAPIDITY) < 0.2; });
      // Remove any untagged low-multiplicity/muon-dominated jet within dR = 0.4 of a muon
      for (const Particle& m : muons)
        ifilter_discard(jets, [&](const Jet& j) {
            if (!j.bTagged(Cuts::pT > 5*GeV)) return false; /// @note A different b-tag working point, 85%, was actually used here *sigh*
            if (deltaR(m, j, RAPIDITY) > 0.4) return false;
            if (j.particles(Cuts::abscharge != 0).size() < 3) return true;
            return m.pT()/j.pT() > 0.5;
          });
      // Removing leptons within dR = 0.4 of remaining jets
      for (const Jet& j : jets) {
        ifilter_discard(elecs, deltaRLess(j, 0.4, RAPIDITY));
        ifilter_discard(muons, deltaRLess(j, 0.4, RAPIDITY));
      }

      // Signal jet and lepton selection
      const Jets sigjets40 = filter_select(jets, Cuts::pT > 40*GeV);
      const Jets sigjets60 = filter_select(sigjets40, Cuts::pT > 60*GeV);
      const Jets sigbjets40 = filter_select(sigjets40, [](const Jet& j) { return j.bTagged(Cuts::pT > 5*GeV); });
      const Jets sigbjets60 = filter_select(sigjets60, [](const Jet& j) { return j.bTagged(Cuts::pT > 5*GeV); });
      const Particles sigmuons = filter_select(muons, Cuts::pT > 35*GeV);
      Particles sigelecs = filter_select(elecs, Cuts::pT > 35*GeV);
      ifilter_select(sigelecs, ParticleEffFilter(ELECTRON_IDEFF_ATLAS_RUN2_TIGHT));


      //////////////////


      // Event selection cuts
      if (sigelecs.size() + sigmuons.size() != 1) vetoEvent;
      const Particle siglepton = sigelecs.empty() ? sigmuons.front() : sigelecs.front();

      /// @note The note describes Nj = 5, 6, 7, 8, 9, >= 10 and Nb = 0, 1, 2, 3, >= 4 = 30 2D bins
      ///  for each jet cut... but only provides data for six Nj = >= 8, 9, 10, Nb = 0, >= 3 bins.
      /// We just implement the latter for now.

      // Fill counters
      const double w = event.weight();
      if (sigjets40.size() >= 8  && sigbjets40.empty()) _h_08j40_0b->fill(w);
      if (sigjets40.size() >= 9  && sigbjets40.empty()) _h_09j40_0b->fill(w);
      if (sigjets40.size() >= 10 && sigbjets40.empty()) _h_10j40_0b->fill(w);
      if (sigjets40.size() >= 8  && sigbjets40.size() >= 3) _h_08j40_3b->fill(w);
      if (sigjets40.size() >= 9  && sigbjets40.size() >= 3) _h_09j40_3b->fill(w);
      if (sigjets40.size() >= 10 && sigbjets40.size() >= 3) _h_10j40_3b->fill(w);

      if (sigjets60.size() >= 8  && sigbjets60.empty()) _h_08j60_0b->fill(w);
      if (sigjets60.size() >= 9  && sigbjets60.empty()) _h_09j60_0b->fill(w);
      if (sigjets60.size() >= 10 && sigbjets60.empty()) _h_10j60_0b->fill(w);
      if (sigjets60.size() >= 8  && sigbjets60.size() >= 3) _h_08j60_3b->fill(w);
      if (sigjets60.size() >= 9  && sigbjets60.size() >= 3) _h_09j60_3b->fill(w);
      if (sigjets60.size() >= 10 && sigbjets60.size() >= 3) _h_10j60_3b->fill(w);

    }


    /// Normalise counters after the run
    void finalize() {

      const double sf = 14.8*crossSection()/femtobarn/sumOfWeights();
      scale({_h_08j40_0b, _h_09j40_0b, _h_10j40_0b, _h_08j40_3b, _h_09j40_3b, _h_10j40_3b}, sf);
      scale({_h_08j60_0b, _h_09j60_0b, _h_10j60_0b, _h_08j60_3b, _h_09j60_3b, _h_10j60_3b}, sf);

    }

    //@}


  private:

    /// @name Histograms
    //@{
    CounterPtr _h_08j40_0b, _h_09j40_0b, _h_10j40_0b, _h_08j40_3b, _h_09j40_3b, _h_10j40_3b;
    CounterPtr _h_08j60_0b, _h_09j60_0b, _h_10j60_0b, _h_08j60_3b, _h_09j60_3b, _h_10j60_3b;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_2016_CONF_2016_094);


}
