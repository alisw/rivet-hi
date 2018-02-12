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


  /// @brief CMS 2016 0-lepton SUSY search, from 13/fb PAS note
  class CMS_2016_PAS_SUS_16_14 : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(CMS_2016_PAS_SUS_16_14);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      FinalState calofs(Cuts::abseta < 5.0);
      FastJets fj(calofs, FastJets::ANTIKT, 0.4);
      declare(fj, "TruthJets");
      declare(SmearedJets(fj, JET_SMEAR_CMS_RUN2, [](const Jet& j) {
            if (j.abseta() > 2.5) return 0.;
            return j.bTagged() ? 0.55 : j.cTagged() ? 0.12 : 0.016; }), "Jets");

      FinalState es(Cuts::abspid == PID::ELECTRON && Cuts::abseta < 2.5);
      declare(es, "TruthElectrons");
      declare(SmearedParticles(es, ELECTRON_EFF_CMS_RUN2, ELECTRON_SMEAR_CMS_RUN2), "Electrons");

      FinalState mus(Cuts::abspid == PID::MUON && Cuts::abseta < 2.4);
      declare(mus, "TruthMuons");
      declare(SmearedParticles(mus, MUON_EFF_CMS_RUN2, MUON_SMEAR_CMS_RUN2), "Muons");

      FinalState isofs(Cuts::abseta < 3.0 && Cuts::abspid != PID::ELECTRON && Cuts::abspid != PID::MUON);
      declare(isofs, "IsoFS");
      FinalState cfs(Cuts::abseta < 2.5 && Cuts::abscharge != 0);
      declare(cfs, "TruthTracks");
      declare(SmearedParticles(cfs, TRK_EFF_CMS_RUN2), "Tracks");

      // Book histograms/counters
      _h_srcounts.resize(160);
      for (size_t ij = 0; ij < 4; ++ij) {
        for (size_t ib = 0; ib < 4; ++ib) {
          for (size_t ih = 0; ih < 10; ++ih) {
            const size_t i = 40*ij + 10*ib + ih;
            _h_srcounts[i] = bookCounter(toString(2*ij+3) + "j-" + toString(ib) + "b-" + toString(ih));
          }
        }
      }
      _h_srcountsagg.resize(12);
      for (size_t ia = 0; ia < 12; ++ia) {
        _h_srcountsagg[ia] = bookCounter("agg-" + toString(ia));
      }

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      // Get jets and require Nj >= 3
      const Jets jets24 = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 2.4);
      if (jets24.size() < 3) vetoEvent;

      // HT cut
      vector<double> jetpts24; transform(jets24, jetpts24, pT);
      const double ht = sum(jetpts24, 0.0);
      if (ht < 300*GeV) vetoEvent;

      // HTmiss cut
      const Jets jets50 = apply<JetAlg>(event, "Jets").jetsByPt(Cuts::pT > 30*GeV && Cuts::abseta < 5.0);
      const FourMomentum htmissvec = -sum(jets24, mom, FourMomentum());
      const double htmiss = htmissvec.pT();
      if (htmissvec.pT() < 300*GeV) vetoEvent;


      // Get baseline electrons & muons
      Particles elecs = apply<ParticleFinder>(event, "Electrons").particles(Cuts::pT > 10*GeV);
      Particles muons = apply<ParticleFinder>(event, "Muons").particles(Cuts::pT > 10*GeV);

      // Electron/muon isolation
      const Particles calofs = apply<ParticleFinder>(event, "IsoFS").particles();
      ifilter_discard(elecs, [&](const Particle& e) {
          const double R = max(0.05, min(0.2, 10*GeV/e.pT()));
          double ptsum = -e.pT();
          for (const Particle& p : calofs)
            if (deltaR(p,e) < R) ptsum += p.pT();
          return ptsum / e.pT() > 0.1;
        });
      ifilter_discard(muons, [&](const Particle& m) {
          const double R = max(0.05, min(0.2, 10*GeV/m.pT()));
          double ptsum = -m.pT();
          for (const Particle& p : calofs)
            if (deltaR(p,m) < R) ptsum += p.pT();
          return ptsum / m.pT() > 0.2;
        });

      // Veto the event if there are any remaining baseline leptons
      if (!elecs.empty()) vetoEvent;
      if (!muons.empty()) vetoEvent;


      // Get isolated tracks
      Particles trks25 = apply<ParticleFinder>(event, "Tracks").particles();
      ifilter_discard(trks25, [&](const Particle& t) {
          double ptsum = -t.pT();
          for (const Particle& p : trks25)
            if (deltaR(p,t) < 0.3) ptsum += p.pT();
          return ptsum/t.pT() > ((t.abspid() == PID::ELECTRON || t.abspid() == PID::MUON) ? 0.2 : 0.1);
        });
      const Particles trks = filter_select(trks25, Cuts::abseta < 2.4);

      // Isolated track pT, pTmiss and mT cut
      // mT^2 = m1^2 + m2^2 + 2(ET1 ET2 - pT1 . pT2))
      // => mT0^2 = 2(ET1 |pT2| - pT1 . pT2)) for m1, m2 -> 0
      FourMomentum ptmissvec = htmissvec; ///< @todo Can we do better? No e,mu left...
      const double ptmiss = ptmissvec.pT();
      for (const Particle& t : trks) {
        const double ptcut = (t.abspid() == PID::ELECTRON || t.abspid() == PID::MUON) ? 5*GeV : 10*GeV;
        const double mT = sqrt( t.mass2() + 2*(t.Et()*ptmiss - t.pT()*ptmiss*cos(deltaPhi(t,ptmissvec))) );
        if (mT < 100*GeV && t.pT() < ptcut) vetoEvent;
      }

      // Lead jets isolation from Htmiss
      if (deltaPhi(htmissvec, jets24[0]) < 0.5) vetoEvent;
      if (deltaPhi(htmissvec, jets24[1]) < 0.5) vetoEvent;
      if (deltaPhi(htmissvec, jets24[2]) < 0.3) vetoEvent;
      if (jets24.size() >= 4 && deltaPhi(htmissvec, jets24[3]) < 0.3) vetoEvent;

      // Count jet and b-jet multiplicities
      const size_t nj = jets24.size();
      size_t nbj = 0;
      for (const Jet& j : jets24)
        if (j.bTagged()) nbj += 1;


      ////////


      // Fill the aggregate signal regions first
      if (nj >= 3 && nbj == 0 && ht >  500*GeV && htmiss > 500*GeV) _h_srcountsagg[ 0]->fill(event.weight());
      if (nj >= 3 && nbj == 0 && ht > 1500*GeV && htmiss > 750*GeV) _h_srcountsagg[ 1]->fill(event.weight());
      if (nj >= 5 && nbj == 0 && ht >  500*GeV && htmiss > 500*GeV) _h_srcountsagg[ 2]->fill(event.weight());
      if (nj >= 5 && nbj == 0 && ht > 1500*GeV && htmiss > 750*GeV) _h_srcountsagg[ 3]->fill(event.weight());
      if (nj >= 9 && nbj == 0 && ht > 1500*GeV && htmiss > 750*GeV) _h_srcountsagg[ 4]->fill(event.weight());
      if (nj >= 3 && nbj >= 2 && ht >  500*GeV && htmiss > 500*GeV) _h_srcountsagg[ 5]->fill(event.weight());
      if (nj >= 3 && nbj >= 1 && ht >  750*GeV && htmiss > 750*GeV) _h_srcountsagg[ 6]->fill(event.weight());
      if (nj >= 5 && nbj >= 3 && ht >  500*GeV && htmiss > 500*GeV) _h_srcountsagg[ 7]->fill(event.weight());
      if (nj >= 5 && nbj >= 2 && ht > 1500*GeV && htmiss > 750*GeV) _h_srcountsagg[ 8]->fill(event.weight());
      if (nj >= 9 && nbj >= 3 && ht >  750*GeV && htmiss > 750*GeV) _h_srcountsagg[ 9]->fill(event.weight());
      if (nj >= 7 && nbj >= 1 && ht >  300*GeV && htmiss > 300*GeV) _h_srcountsagg[10]->fill(event.weight());
      if (nj >= 5 && nbj >= 1 && ht >  750*GeV && htmiss > 750*GeV) _h_srcountsagg[11]->fill(event.weight());


      // Nj bin and Nbj bins
      static const vector<double> njedges = {3., 5., 7., 9.};
      const size_t inj = binIndex(nj, njedges, true);
      static const vector<double> njbedges = {0., 1., 2., 3.};
      const size_t inbj = binIndex(nbj, njbedges, true);
      // HTmiss vs HT 2D bin
      int iht = 0;
      if (htmiss < 350*GeV) {
        iht = ht < 500 ? 1 : ht < 1000 ? 2 : 3;
      } else if (htmiss < 500*GeV && ht > 350*GeV) {
        iht = ht < 500 ? 4 : ht < 1000 ? 5 : 6;
      } else if (htmiss < 750*GeV && ht > 500*GeV) {
        iht = ht < 1000 ? 7 : 8;
      } else if (ht > 750*GeV) {
        iht = ht < 1500 ? 9 : 10;
      }
      if (iht == 0) vetoEvent;
      iht -= 1; //< change from the paper's indexing scheme to C++ zero-indexed
      // Total bin number
      const size_t ibin = 40*inj + 10*inbj + (size_t)iht;

      // Fill SR counter
      _h_srcounts[ibin]->fill(event.weight());

    }


    /// Normalise counters after the run
    void finalize() {

      const double sf = 12.9*crossSection()/femtobarn/sumOfWeights();
      scale(_h_srcounts, sf);
      scale(_h_srcountsagg, sf);

    }

    //@}


  private:

    /// @name Histograms
    //@{
    vector<CounterPtr> _h_srcounts, _h_srcountsagg;
    //@}


  };



  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(CMS_2016_PAS_SUS_16_14);


}
